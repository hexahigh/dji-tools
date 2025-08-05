# Long Help

## Overview
This tool extracts frames from DJI drone videos and tags them with GPS metadata from the embedded subtitles.

## Flags

### Input (Required)
`--input, -i`
Specifies the input video file.
Example: `-i DJI_0001.MP4`

### Extract (Required)
`--extract`
Defines the frame extraction mode. Supports three modes:

1. **Fixed FPS Extraction**
   `--extract fps=N`
   Extracts frames at a specified frame rate (N frames per second).
   Example: `--extract fps=30` extracts 30 frames per second.

2. **All Frames Extraction**
   `--extract all`
   Extracts every frame from the video.
   *Note: May create a large number of files for long videos.*

3. **GPS-Linked Extraction**
   `--extract direct`
   Extracts only frames that have associated GPS data in the subtitle track.
   *Note: Requires subtitle track with GPS coordinates.*

### Tag (Required)
`--tag`
Determines how GPS metadata is applied to extracted frames:

1. **Direct Tagging**
   `--tag direct`
   Tags exactly one frame per GPS entry - the frame closest to the GPS timestamp.
   *Works with all extraction modes.*

2. **Closest Match Tagging**
   `--tag all`
   Tags all frames with the closest GPS data point in time.

3. **Interpolated Tagging**
   `--tag allip`
   Tags all frames with interpolated GPS coordinates between adjacent data points.

4. **No Tagging**
   `--tag none`
   Skips metadata tagging entirely (only extracts frames).


### Height Options
`--use-height`
Enables altitude/height tagging in images. Note that the altitudes stored in the subtitles are relative to the takeoff position.
Example: `--use-height`

`--height-offset`
Applies an offset to height values (in meters). Useful for adjusting for takeoff position.
Format: Floating point value
Example: `--height-offset=32.5` adds 32.5 meters to all heights

### Help
`--long-help, -H`
Displays this detailed help documentation.

## Constraints & Compatibility
- When using `--extract direct`, only these tagging modes are allowed:
  - `--tag direct` (tags extracted frames)
  - `--tag none` (just extracts frames)
- For `--extract fps` and `--extract all`, all tagging modes are supported

## File Naming Conventions
- **Direct extraction**: `direct_HH_MM_SS_sss.jpg`
- **FPS/All extraction**: `frame_00000001.jpg`
- Output directory: `[video_name]_jpegs/`

## Processing Workflow
1. Extract subtitles from video
2. Parse GPS coordinates from subtitle track
3. Extract frames according to selected mode
4. Tag frames with GPS metadata per selected method

## Example Commands

1. Extract GPS-linked frames and tag them directly:
```bash
program -i video.mp4 --extract direct --tag direct
```
2. Extract at 5 fps and tag closest GPS for all frames:
```bash
program -i video.mp4 --extract fps=5 --tag all
```
3. Extract all frames with interpolated GPS tags:
```bash
program -i video.mp4 --extract all --tag allip
```
4. Extract GPS-linked frames without tagging:
```bash
program -i video.mp4 --extract direct --tag none
```

## Dependencies
- `ffmpeg`: Video processing and frame extraction
- `ffprobe`: Video metadata analysis
- `exiftool`: GPS metadata embedding in JPEGs

## Output Structure
```
video_name_jpegs/
├── direct_12_34_56_789.jpg
├── direct_12_35_01_123.jpg
├── frame_00000001.jpg
├── frame_00000002.jpg
└── ...
```

## Troubleshooting
- Ensure video has embedded subtitle track with GPS data
- Verify all dependencies are in system PATH
- Check file permissions for output directory
- Use absolute paths for input files if encountering path issues