package main

import (
	"bufio"
	_ "embed"
	"errors"
	"fmt"
	"math"
	"os"
	"os/exec"
	"path/filepath"
	"regexp"
	"sort"
	"strconv"
	"strings"

	"github.com/charmbracelet/glamour"
	"github.com/spf13/pflag"
)

type GPSEntry struct {
	TimeStr string
	TimeSec float64
	Lat     float64
	Lon     float64
	Height  float64 // The height in the subtitle file is relative to the takeoff position
}

var ffmpegBase = "ffmpeg"
var extraFfmpegArgs = []string{"-hide_banner", "-v", "warning"}

type FrameInfo struct {
	Path    string
	TimeSec float64
}

type Options = struct {
	inputVideo   string
	extractMode  string
	tagMode      string
	longHelp     bool
	useHeight    bool
	heightOffset float64
}

var options Options

func timeStrToSeconds(timeStr string) (float64, error) {
	parts := strings.Split(timeStr, ":")
	if len(parts) != 3 {
		return 0, fmt.Errorf("invalid time format: %s", timeStr)
	}

	secMs := strings.Split(parts[2], ",")
	if len(secMs) != 2 {
		return 0, fmt.Errorf("invalid seconds format: %s", parts[2])
	}

	hours, err := strconv.Atoi(parts[0])
	if err != nil {
		return 0, err
	}
	minutes, err := strconv.Atoi(parts[1])
	if err != nil {
		return 0, err
	}
	seconds, err := strconv.Atoi(secMs[0])
	if err != nil {
		return 0, err
	}
	millis, err := strconv.Atoi(secMs[1])
	if err != nil {
		return 0, err
	}

	totalSec := float64(hours*3600+minutes*60+seconds) + float64(millis)/1000.0
	return totalSec, nil
}

func extractSubtitles(videoPath, srtPath string) error {
	cmd := exec.Command(ffmpegBase, append(extraFfmpegArgs, "-y", "-i", videoPath, "-map", "0:s:0", srtPath)...)
	cmd.Stdout = os.Stdout
	cmd.Stderr = os.Stderr
	return cmd.Run()
}

func parseSRT(srtPath string) ([]GPSEntry, error) {
	file, err := os.Open(srtPath)
	if err != nil {
		return nil, err
	}
	defer file.Close()

	timeRe := regexp.MustCompile(`^(\d{2}:\d{2}:\d{2},\d{3})`)
	gpsRe := regexp.MustCompile(`GPS\s*\(\s*([-0-9.]+)\s*,\s*([-0-9.]+)`)
	heightRe := regexp.MustCompile(`H\s*([-\d.]+)m`)

	entries := []GPSEntry{}
	scanner := bufio.NewScanner(file)
	var currentTime string

	for scanner.Scan() {
		line := strings.TrimSpace(scanner.Text())

		if matches := timeRe.FindStringSubmatch(line); matches != nil {
			currentTime = matches[1]
		} else if gpsMatch := gpsRe.FindStringSubmatch(line); gpsMatch != nil && currentTime != "" {
			lon, err := strconv.ParseFloat(gpsMatch[1], 64)
			if err != nil {
				continue
			}
			lat, err := strconv.ParseFloat(gpsMatch[2], 64)
			if err != nil {
				continue
			}
			height := 0.0
			if hMatch := heightRe.FindStringSubmatch(line); hMatch != nil {
				h, err := strconv.ParseFloat(hMatch[1], 64)
				if err == nil {
					height = h
				}
			}
			timeSec, err := timeStrToSeconds(currentTime)
			if err != nil {
				continue
			}
			entries = append(entries, GPSEntry{
				TimeStr: currentTime,
				TimeSec: timeSec,
				Lat:     lat,
				Lon:     lon,
				Height:  height,
			})
			currentTime = ""
		}
	}

	if err := scanner.Err(); err != nil {
		return nil, err
	}

	return entries, nil
}

func getAvgFPS(videoPath string) (float64, error) {
	cmd := exec.Command("ffprobe", "-v", "0", "-select_streams", "v:0",
		"-show_entries", "stream=r_frame_rate", "-of", "csv=p=0", videoPath)
	output, err := cmd.Output()
	if err != nil {
		return 0, err
	}

	rateStr := strings.TrimSpace(string(output))
	parts := strings.Split(rateStr, "/")
	if len(parts) != 2 {
		return 0, fmt.Errorf("invalid frame rate format: %s", rateStr)
	}

	num, err := strconv.ParseFloat(parts[0], 64)
	if err != nil {
		return 0, err
	}
	den, err := strconv.ParseFloat(parts[1], 64)
	if err != nil {
		return 0, err
	}

	if den == 0 {
		return 0, errors.New("division by zero in frame rate")
	}
	return num / den, nil
}

func tagImage(imagePath string, lat, lon, height float64) error {
	latAbs := math.Abs(lat)
	lonAbs := math.Abs(lon)

	latRef := "N"
	if lat < 0 {
		latRef = "S"
	}
	lonRef := "E"
	if lon < 0 {
		lonRef = "W"
	}

	// Base arguments for latitude/longitude
	args := []string{
		"-overwrite_original",
		fmt.Sprintf("-GPSLatitude=%.9f", latAbs),
		fmt.Sprintf("-GPSLatitudeRef=%s", latRef),
		fmt.Sprintf("-GPSLongitude=%.9f", lonAbs),
		fmt.Sprintf("-GPSLongitudeRef=%s", lonRef),
	}

	// Add height tags if enabled
	if options.useHeight {
		adjustedHeight := height + options.heightOffset
		absAlt := math.Abs(adjustedHeight)
		altRef := "0" // Above sea level
		if adjustedHeight < 0 {
			altRef = "1" // Below sea level
		}
		args = append(args,
			fmt.Sprintf("-GPSAltitude=%.3f", absAlt),
			fmt.Sprintf("-GPSAltitudeRef=%s", altRef),
		)
	}

	args = append(args, imagePath)
	cmd := exec.Command("exiftool", args...)
	cmd.Stdout = os.Stdout
	cmd.Stderr = os.Stderr
	return cmd.Run()
}

func findClosestFrame(frames []FrameInfo, target float64) (FrameInfo, error) {
	closest := frames[0]
	minDiff := math.Abs(frames[0].TimeSec - target)

	for _, frame := range frames[1:] {
		diff := math.Abs(frame.TimeSec - target)
		if diff < minDiff {
			minDiff = diff
			closest = frame
		}
	}
	return closest, nil
}

// Helper function to find closest GPS entry
func findClosestGPS(entries []GPSEntry, target float64) (GPSEntry, error) {
	closest := entries[0]
	minDiff := math.Abs(entries[0].TimeSec - target)

	for _, entry := range entries[1:] {
		diff := math.Abs(entry.TimeSec - target)
		if diff < minDiff {
			minDiff = diff
			closest = entry
		}
	}
	return closest, nil
}

const earthRadius = 6371e3 // meters

// haversineDistance computes the great-circle distance between two GPS coordinates
// using the haversine formula.
func haversineDistance(lat1, lon1, lat2, lon2 float64) float64 {
	lat1R := lat1 * math.Pi / 180
	lon1R := lon1 * math.Pi / 180
	lat2R := lat2 * math.Pi / 180
	lon2R := lon2 * math.Pi / 180

	dLat := lat2R - lat1R
	dLon := lon2R - lon1R

	a := math.Sin(dLat/2)*math.Sin(dLat/2) +
		math.Cos(lat1R)*math.Cos(lat2R)*math.Sin(dLon/2)*math.Sin(dLon/2)
	return earthRadius * 2 * math.Asin(math.Sqrt(a))
}

// gpsCalculateDistances returns the cumulative haversine distance for each GPS entry.
// If useEle is true, height differences are incorporated into the distance.
func gpsCalculateDistances(entries []GPSEntry, useEle bool) []float64 {
	dists := make([]float64, len(entries))
	for i := 1; i < len(entries); i++ {
		d := haversineDistance(entries[i-1].Lat, entries[i-1].Lon, entries[i].Lat, entries[i].Lon)
		if useEle {
			dEle := entries[i].Height - entries[i-1].Height
			d = math.Sqrt(d*d + dEle*dEle)
		}
		dists[i] = dists[i-1] + d
	}
	return dists
}

// gpsRemoveDuplicates removes GPS entries whose horizontal position is identical
// to the previous entry (zero haversine distance)
func gpsRemoveDuplicates(entries []GPSEntry) []GPSEntry {
	if len(entries) <= 1 {
		return entries
	}
	dists := gpsCalculateDistances(entries, false)
	result := []GPSEntry{entries[0]}
	for i := 1; i < len(entries); i++ {
		if dists[i] > dists[i-1] {
			result = append(result, entries[i])
		}
	}
	return result
}

// computePCHIPSlopes computes tangent slopes for a PCHIP spline using the
// Fritsch-Carlson algorithm, which guarantees monotonicity within each interval.
func computePCHIPSlopes(x, y []float64) []float64 {
	n := len(x)
	slopes := make([]float64, n)
	if n < 2 {
		return slopes
	}

	delta := make([]float64, n-1)
	for i := range delta {
		delta[i] = (y[i+1] - y[i]) / (x[i+1] - x[i])
	}

	// One-sided endpoint slopes
	slopes[0] = delta[0]
	slopes[n-1] = delta[n-2]

	// Interior slopes via weighted harmonic mean (Fritsch-Carlson)
	for i := 1; i < n-1; i++ {
		if delta[i-1]*delta[i] <= 0 {
			slopes[i] = 0 // local extremum: preserve shape
		} else {
			h0 := x[i] - x[i-1]
			h1 := x[i+1] - x[i]
			w0 := 2*h1 + h0
			w1 := h1 + 2*h0
			slopes[i] = (w0 + w1) / (w0/delta[i-1] + w1/delta[i])
		}
	}

	// Clamp slopes to enforce monotonicity (|alpha|^2 + |beta|^2 <= 9)
	for i := 0; i < n-1; i++ {
		if delta[i] == 0 {
			slopes[i] = 0
			slopes[i+1] = 0
		} else {
			alpha := slopes[i] / delta[i]
			beta := slopes[i+1] / delta[i]
			tau := math.Sqrt(alpha*alpha + beta*beta)
			if tau > 3 {
				t := 3.0 / tau
				slopes[i] = t * alpha * delta[i]
				slopes[i+1] = t * beta * delta[i]
			}
		}
	}

	return slopes
}

// pchipEval evaluates a PCHIP spline (defined by knots x, values y, slopes m) at xq.
func pchipEval(x, y, m []float64, xq float64) float64 {
	n := len(x)
	if xq <= x[0] {
		return y[0]
	}
	if xq >= x[n-1] {
		return y[n-1]
	}

	i := sort.SearchFloat64s(x, xq) - 1

	h := x[i+1] - x[i]
	t := (xq - x[i]) / h

	// Cubic Hermite basis polynomials
	h00 := (1 + 2*t) * (1 - t) * (1 - t)
	h10 := t * (1 - t) * (1 - t)
	h01 := t * t * (3 - 2*t)
	h11 := t * t * (t - 1)

	return h00*y[i] + h10*h*m[i] + h01*y[i+1] + h11*h*m[i+1]
}

// interpolateGPS returns the interpolated GPS position at the given target time
// using piecewise cubic Hermite splines (PCHIP) parameterized by cumulative
// haversine distance
func interpolateGPS(entries []GPSEntry, target float64) (lat, lon, height float64, err error) {
	if len(entries) == 0 {
		return 0, 0, 0, errors.New("no GPS entries")
	}
	if target <= entries[0].TimeSec {
		return entries[0].Lat, entries[0].Lon, entries[0].Height, nil
	}
	if target >= entries[len(entries)-1].TimeSec {
		e := entries[len(entries)-1]
		return e.Lat, e.Lon, e.Height, nil
	}

	// Remove spatially duplicate points so cumulative distance is strictly increasing
	deduped := gpsRemoveDuplicates(entries)
	if len(deduped) == 1 {
		return deduped[0].Lat, deduped[0].Lon, deduped[0].Height, nil
	}

	n := len(deduped)
	cumDist := gpsCalculateDistances(deduped, true)

	lats := make([]float64, n)
	lons := make([]float64, n)
	heights := make([]float64, n)
	times := make([]float64, n)
	for i, e := range deduped {
		lats[i] = e.Lat
		lons[i] = e.Lon
		heights[i] = e.Height
		times[i] = e.TimeSec
	}

	// Build PCHIP splines for lat, lon, height, and time — all as functions of cumulative distance
	mLat := computePCHIPSlopes(cumDist, lats)
	mLon := computePCHIPSlopes(cumDist, lons)
	mHeight := computePCHIPSlopes(cumDist, heights)
	mTime := computePCHIPSlopes(cumDist, times)

	// Invert the time spline to find the cumulative distance d* at which time == target.
	// Time is monotonically increasing with distance after duplicate removal, so bisection converges.
	dLow, dHigh := cumDist[0], cumDist[n-1]
	for range 64 {
		dMid := (dLow + dHigh) / 2
		if pchipEval(cumDist, times, mTime, dMid) < target {
			dLow = dMid
		} else {
			dHigh = dMid
		}
		if dHigh-dLow < 1e-9 {
			break
		}
	}
	dTarget := (dLow + dHigh) / 2

	return pchipEval(cumDist, lats, mLat, dTarget),
		pchipEval(cumDist, lons, mLon, dTarget),
		pchipEval(cumDist, heights, mHeight, dTarget),
		nil
}

func extractFrames(videoPath, outDir, pattern string, args ...string) error {
	cmdArgs := extraFfmpegArgs
	cmdArgs = append(cmdArgs, "-y", "-i", videoPath)
	cmdArgs = append(cmdArgs, args...)
	cmdArgs = append(cmdArgs, filepath.Join(outDir, pattern))
	cmd := exec.Command(ffmpegBase, cmdArgs...)
	cmd.Stdout = os.Stdout
	cmd.Stderr = os.Stderr
	return cmd.Run()
}

func getFrameTimes(outDir string, fps float64) ([]FrameInfo, error) {
	files, err := filepath.Glob(filepath.Join(outDir, "frame_*.jpg"))
	if err != nil {
		return nil, err
	}

	frames := make([]FrameInfo, 0, len(files))
	for _, file := range files {
		base := filepath.Base(file)
		frameNumStr := strings.TrimPrefix(strings.TrimSuffix(base, ".jpg"), "frame_")
		frameNum, err := strconv.ParseInt(frameNumStr, 10, 64)
		if err != nil {
			continue
		}
		frames = append(frames, FrameInfo{
			Path:    file,
			TimeSec: float64(frameNum) / fps,
		})
	}

	// Sort by time
	sort.Slice(frames, func(i, j int) bool {
		return frames[i].TimeSec < frames[j].TimeSec
	})

	return frames, nil
}

func main() {
	pflag.StringVarP(&options.inputVideo, "input", "i", "", "Input video file")
	pflag.StringVar(&options.extractMode, "extract", "", "Extraction mode: fps=N, all, or direct")
	pflag.StringVar(&options.tagMode, "tag", "", "Tagging mode: direct, all, allip, none")
	pflag.BoolVarP(&options.longHelp, "long-help", "H", false, "Display long help")
	pflag.BoolVar(&options.useHeight, "use-height", false, "Enable height tagging in images")
	pflag.Float64Var(&options.heightOffset, "height-offset", 0.0, "Height offset to apply (meters)")
	pflag.Parse()

	inputVideo := options.inputVideo
	extractMode := options.extractMode
	tagMode := options.tagMode

	if options.longHelp {
		generateLongHelp()
		os.Exit(0)
	}

	if inputVideo == "" {
		fmt.Println("Error: Input video file is required")
		os.Exit(1)
	}
	if extractMode == "" || tagMode == "" {
		fmt.Println("Error: --extract and --tag flags are required")
		os.Exit(1)
	}

	if extractMode == "direct" && (tagMode != "direct" && tagMode != "none") {
		fmt.Println("Error: When using --extract direct, --tag must be direct or none")
		os.Exit(1)
	}

	deps := []string{"exiftool", "ffmpeg", "ffprobe"}
	for _, dep := range deps {
		if _, err := exec.LookPath(dep); err != nil && errors.Is(err, exec.ErrNotFound) {
			fmt.Printf("Error: dependency %s not found in path\n", dep)
			os.Exit(1)
		}
	}

	base := strings.TrimSuffix(inputVideo, filepath.Ext(inputVideo))
	srtPath := base + ".srt"
	outDir := base + "_jpegs"

	if err := os.MkdirAll(outDir, 0755); err != nil {
		fmt.Printf("Error creating output directory: %v\n", err)
		os.Exit(1)
	}

	fmt.Println("Extracting subtitles...")
	if err := extractSubtitles(inputVideo, srtPath); err != nil {
		fmt.Printf("Error extracting subtitles: %v\n", err)
		os.Exit(1)
	}

	fmt.Println("Parsing SRT entries...")
	entries, err := parseSRT(srtPath)
	if err != nil {
		fmt.Printf("Error parsing SRT: %v\n", err)
		os.Exit(1)
	}

	sort.Slice(entries, func(i, j int) bool {
		return entries[i].TimeSec < entries[j].TimeSec
	})

	var frames []FrameInfo
	var fps float64

	switch {
	case extractMode == "all":
		fmt.Println("Extracting all frames...")
		fps, err = getAvgFPS(inputVideo)
		if err != nil {
			fmt.Printf("Error getting FPS: %v\n", err)
			os.Exit(1)
		}
		if err := extractFrames(inputVideo, outDir, "frame_%08d.jpg", "-vsync", "0", "-start_number", "0", "-q:v", "2"); err != nil {
			fmt.Printf("Error extracting frames: %v\n", err)
			os.Exit(1)
		}
		frames, err = getFrameTimes(outDir, fps)
		if err != nil {
			fmt.Printf("Error getting frame times: %v\n", err)
			os.Exit(1)
		}

	case strings.HasPrefix(extractMode, "fps="):
		fpsStr := strings.TrimPrefix(extractMode, "fps=")
		fps, err = strconv.ParseFloat(fpsStr, 64)
		if err != nil {
			fmt.Printf("Invalid FPS value: %v\n", err)
			os.Exit(1)
		}
		fmt.Printf("Extracting at %.2f fps...\n", fps)
		if err := extractFrames(inputVideo, outDir, "frame_%08d.jpg", "-r", fmt.Sprintf("%.2f", fps), "-start_number", "0", "-q:v", "2"); err != nil {
			fmt.Printf("Error extracting frames: %v\n", err)
			os.Exit(1)
		}
		frames, err = getFrameTimes(outDir, fps)
		if err != nil {
			fmt.Printf("Error getting frame times: %v\n", err)
			os.Exit(1)
		}

	case extractMode == "direct":
		fmt.Println("Extracting GPS-linked frames...")
		for _, entry := range entries {
			timeFF := strings.Replace(entry.TimeStr, ",", ".", -1)
			frameName := fmt.Sprintf("direct_%s.jpg", strings.ReplaceAll(strings.ReplaceAll(entry.TimeStr, ":", "_"), ",", "_"))
			outPath := filepath.Join(outDir, frameName)

			cmd := exec.Command(ffmpegBase, append(extraFfmpegArgs, "-y", "-ss", timeFF, "-i", inputVideo,
				"-vframes", "1", "-q:v", "2", outPath)...)
			cmd.Stdout = os.Stdout
			cmd.Stderr = os.Stderr
			if err := cmd.Run(); err != nil {
				fmt.Printf("Error extracting frame at %s: %v\n", timeFF, err)
				continue
			}
			frames = append(frames, FrameInfo{
				Path:    outPath,
				TimeSec: entry.TimeSec,
			})
		}

	default:
		fmt.Printf("Unknown extraction mode: %s\n", extractMode)
		os.Exit(1)
	}

	// Tagging logic
	switch tagMode {
	case "direct":
		fmt.Println("Tagging direct frames...")
		for _, entry := range entries {
			closest, err := findClosestFrame(frames, entry.TimeSec)
			if err != nil {
				fmt.Printf("Error finding frame for time %.2f: %v\n", entry.TimeSec, err)
				continue
			}
			if err := tagImage(closest.Path, entry.Lat, entry.Lon, entry.Height); err != nil {
				fmt.Printf("Error tagging %s: %v\n", filepath.Base(closest.Path), err)
			} else {
				fmt.Printf("Tagged %s: lat=%.6f, lon=%.6f\n", filepath.Base(closest.Path), entry.Lat, entry.Lon)
			}
		}

	case "all":
		fmt.Println("Tagging all frames with closest GPS...")
		for _, frame := range frames {
			gps, err := findClosestGPS(entries, frame.TimeSec)
			if err != nil {
				fmt.Printf("Error finding GPS for frame %s: %v\n", filepath.Base(frame.Path), err)
				continue
			}
			if err := tagImage(frame.Path, gps.Lat, gps.Lon, gps.Height); err != nil {
				fmt.Printf("Error tagging %s: %v\n", filepath.Base(frame.Path), err)
			} else {
				fmt.Printf("Tagged %s: lat=%.6f, lon=%.6f\n", filepath.Base(frame.Path), gps.Lat, gps.Lon)
			}
		}

	case "allip":
		fmt.Println("Tagging all frames with interpolated GPS...")
		for _, frame := range frames {
			lat, lon, height, err := interpolateGPS(entries, frame.TimeSec)
			if err != nil {
				fmt.Printf("Error interpolating GPS for frame %s: %v\n", filepath.Base(frame.Path), err)
				continue
			}
			if err := tagImage(frame.Path, lat, lon, height); err != nil {
				fmt.Printf("Error tagging %s: %v\n", filepath.Base(frame.Path), err)
			} else {
				fmt.Printf("Tagged %s: lat=%.6f, lon=%.6f\n", filepath.Base(frame.Path), lat, lon)
			}
		}

	case "none":
		fmt.Println("Skipping tagging as requested")
	}

	fmt.Println("Operation completed successfully.")
}

//go:embed longhelp.md
var longHelp string

func generateLongHelp() {
	out, err := glamour.Render(longHelp, "dark")
	if err != nil {
		fmt.Printf("Error rendering markdown: %v\n", err)
		fmt.Println(longHelp)
		return
	}
	fmt.Println(out)
}
