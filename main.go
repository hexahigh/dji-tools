package main

import (
	"bufio"
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

	"github.com/spf13/pflag"
)

type GPSEntry struct {
	TimeStr string
	TimeSec float64
	Lat     float64
	Lon     float64
}

var ffmpegBase = "ffmpeg"
var extraFfmpegArgs = []string{"-hide_banner", "-v", "error"}

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
			timeSec, err := timeStrToSeconds(currentTime)
			if err != nil {
				continue
			}
			entries = append(entries, GPSEntry{
				TimeStr: currentTime,
				TimeSec: timeSec,
				Lat:     lat,
				Lon:     lon,
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

func sanitizeFilename(timeStr string) string {
	return strings.ReplaceAll(strings.ReplaceAll(timeStr, ":", "_"), ",", "_")
}

func tagImage(imagePath string, lat, lon float64) error {
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

	cmd := exec.Command("exiftool",
		"-overwrite_original",
		fmt.Sprintf("-GPSLatitude=%.9f", latAbs),
		fmt.Sprintf("-GPSLatitudeRef=%s", latRef),
		fmt.Sprintf("-GPSLongitude=%.9f", lonAbs),
		fmt.Sprintf("-GPSLongitudeRef=%s", lonRef),
		imagePath,
	)
	cmd.Stdout = os.Stdout
	cmd.Stderr = os.Stderr
	return cmd.Run()
}

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

func interpolateGPS(entries []GPSEntry, target float64) (lat, lon float64, err error) {
	if len(entries) == 0 {
		return 0, 0, errors.New("no GPS entries")
	}
	if target <= entries[0].TimeSec {
		return entries[0].Lat, entries[0].Lon, nil
	}
	if target >= entries[len(entries)-1].TimeSec {
		return entries[len(entries)-1].Lat, entries[len(entries)-1].Lon, nil
	}

	idx := sort.Search(len(entries), func(i int) bool {
		return entries[i].TimeSec >= target
	})
	if idx == 0 || idx == len(entries) {
		return entries[idx].Lat, entries[idx].Lon, nil
	}

	prev := entries[idx-1]
	next := entries[idx]

	total := next.TimeSec - prev.TimeSec
	weightPrev := (next.TimeSec - target) / total
	weightNext := (target - prev.TimeSec) / total

	lat = prev.Lat*weightPrev + next.Lat*weightNext
	lon = prev.Lon*weightPrev + next.Lon*weightNext
	return lat, lon, nil
}

func main() {
	var inputVideo, extractMode, tagMode string
	pflag.StringVarP(&inputVideo, "input", "i", "", "Input video file")
	pflag.StringVar(&extractMode, "extract", "", "Extraction mode: fps=N, all, or direct")
	pflag.StringVar(&tagMode, "tag", "", "Tagging mode: direct, all, allip, none")
	pflag.Parse()

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

	switch {
	case extractMode == "all":
		fmt.Println("Extracting all frames...")
		fps, err := getAvgFPS(inputVideo)
		if err != nil {
			fmt.Printf("Error getting FPS: %v\n", err)
			os.Exit(1)
		}
		cmd := exec.Command(ffmpegBase, append(extraFfmpegArgs, "-y", "-i", inputVideo, "-vsync", "0", "-start_number", "0",
			"-q:v", "2", filepath.Join(outDir, "frame_%08d.jpg"))...)
		cmd.Stdout = os.Stdout
		cmd.Stderr = os.Stderr
		if err := cmd.Run(); err != nil {
			fmt.Printf("Error extracting frames: %v\n", err)
			os.Exit(1)
		}
		if tagMode != "none" {
			files, err := filepath.Glob(filepath.Join(outDir, "frame_*.jpg"))
			if err != nil {
				fmt.Printf("Error listing frames: %v\n", err)
				os.Exit(1)
			}
			for _, file := range files {
				base := filepath.Base(file)
				frameNumStr := strings.TrimPrefix(strings.TrimSuffix(base, ".jpg"), "frame_")
				frameNum, err := strconv.ParseInt(frameNumStr, 10, 64)
				if err != nil {
					fmt.Printf("Error parsing frame number: %v\n", err)
					continue
				}
				timeSec := float64(frameNum) / fps
				switch tagMode {
				case "all":
					gps, err := findClosestGPS(entries, timeSec)
					if err != nil {
						fmt.Printf("Error finding GPS for frame %s: %v\n", base, err)
						continue
					}
					if err := tagImage(file, gps.Lat, gps.Lon); err != nil {
						fmt.Printf("Error tagging %s: %v\n", base, err)
					} else {
						fmt.Printf("Tagged %s: lat=%.6f, lon=%.6f\n", base, gps.Lat, gps.Lon)
					}
				case "allip":
					lat, lon, err := interpolateGPS(entries, timeSec)
					if err != nil {
						fmt.Printf("Error interpolating GPS for frame %s: %v\n", base, err)
						continue
					}
					if err := tagImage(file, lat, lon); err != nil {
						fmt.Printf("Error tagging %s: %v\n", base, err)
					} else {
						fmt.Printf("Tagged %s: lat=%.6f, lon=%.6f\n", base, lat, lon)
					}
				}
			}
		}

	case strings.HasPrefix(extractMode, "fps="):
		fpsStr := strings.TrimPrefix(extractMode, "fps=")
		fps, err := strconv.ParseFloat(fpsStr, 64)
		if err != nil {
			fmt.Printf("Invalid FPS value: %v\n", err)
			os.Exit(1)
		}
		fmt.Printf("Extracting at %.2f fps...\n", fps)
		cmd := exec.Command(ffmpegBase, append(extraFfmpegArgs, "-y", "-i", inputVideo, "-r", fmt.Sprintf("%.2f", fps),
			"-start_number", "0", "-q:v", "2", filepath.Join(outDir, "frame_%08d.jpg"))...)
		cmd.Stdout = os.Stdout
		cmd.Stderr = os.Stderr
		if err := cmd.Run(); err != nil {
			fmt.Printf("Error extracting frames: %v\n", err)
			os.Exit(1)
		}
		if tagMode != "none" {
			files, err := filepath.Glob(filepath.Join(outDir, "frame_*.jpg"))
			if err != nil {
				fmt.Printf("Error listing frames: %v\n", err)
				os.Exit(1)
			}
			for _, file := range files {
				base := filepath.Base(file)
				frameNumStr := strings.TrimPrefix(strings.TrimSuffix(base, ".jpg"), "frame_")
				frameNum, err := strconv.ParseInt(frameNumStr, 10, 64)
				if err != nil {
					fmt.Printf("Error parsing frame number: %v\n", err)
					continue
				}
				timeSec := float64(frameNum) / fps
				switch tagMode {
				case "all":
					gps, err := findClosestGPS(entries, timeSec)
					if err != nil {
						fmt.Printf("Error finding GPS for frame %s: %v\n", base, err)
						continue
					}
					if err := tagImage(file, gps.Lat, gps.Lon); err != nil {
						fmt.Printf("Error tagging %s: %v\n", base, err)
					} else {
						fmt.Printf("Tagged %s: lat=%.6f, lon=%.6f\n", base, gps.Lat, gps.Lon)
					}
				case "allip":
					lat, lon, err := interpolateGPS(entries, timeSec)
					if err != nil {
						fmt.Printf("Error interpolating GPS for frame %s: %v\n", base, err)
						continue
					}
					if err := tagImage(file, lat, lon); err != nil {
						fmt.Printf("Error tagging %s: %v\n", base, err)
					} else {
						fmt.Printf("Tagged %s: lat=%.6f, lon=%.6f\n", base, lat, lon)
					}
				}
			}
		}

	case extractMode == "direct":
		fmt.Println("Extracting GPS-linked frames...")
		for _, entry := range entries {
			timeFF := strings.Replace(entry.TimeStr, ",", ".", -1)
			frameName := sanitizeFilename(entry.TimeStr) + ".jpg"
			outPath := filepath.Join(outDir, frameName)
			cmd := exec.Command(ffmpegBase, append(extraFfmpegArgs, "-y", "-ss", timeFF, "-i", inputVideo,
				"-vframes", "1", "-q:v", "2", outPath)...)
			cmd.Stdout = os.Stdout
			cmd.Stderr = os.Stderr
			if err := cmd.Run(); err != nil {
				fmt.Printf("Error extracting frame at %s: %v\n", timeFF, err)
				continue
			}
			if tagMode == "direct" {
				if err := tagImage(outPath, entry.Lat, entry.Lon); err != nil {
					fmt.Printf("Error tagging %s: %v\n", frameName, err)
				} else {
					fmt.Printf("Tagged %s: lat=%.6f, lon=%.6f\n", frameName, entry.Lat, entry.Lon)
				}
			}
		}

	default:
		fmt.Printf("Unknown extraction mode: %s\n", extractMode)
		os.Exit(1)
	}

	fmt.Println("Operation completed successfully.")
}
