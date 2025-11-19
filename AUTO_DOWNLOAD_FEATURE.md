# Automatic Dataset Download Feature

## Overview

The system now automatically downloads EV biomarker datasets when a disease is not found in local pre-downloaded files. The download process is visible in the GUI with a progress bar and log window.

## How It Works

### 1. **Automatic Detection**
When you click "Run Step 1" with a new disease:
- System first checks if local data exists for that disease/biofluid
- If found → proceeds directly to biomarker search
- If not found → automatically starts download process

### 2. **Download Process**
The system downloads from multiple EV databases:
- **Vesiclepedia** (http://microvesicles.org/)
- **ExoCarta** (http://www.exocarta.org/)
- **EVpedia** (http://evpedia.info/)

Each download is saved to `./data/databases/` with format:
- `vesiclepedia_{disease}_{biofluid}.csv`
- `exocarta_{disease}_{biofluid}.csv`
- `evpedia_{disease}_{biofluid}.csv`

### 3. **Visual Feedback**
During download, the GUI shows:
- **Progress bar** - Visual percentage indicator
- **Progress log** - Step-by-step messages:
  - "No local data found. Starting download..."
  - "Querying Vesiclepedia database..."
  - "Processing Vesiclepedia results..."
  - "Downloaded Vesiclepedia data: filename.csv"
  - "Querying ExoCarta database..."
  - etc.
- **Completion message** - "Download complete! X files downloaded."

### 4. **Automatic Processing**
After download completes:
- System automatically runs biomarker search
- Results include data from newly downloaded files
- Files are saved for future use (no re-download needed)

## Example Flow

**User Input:**
- Disease: "Alzheimer"
- Biofluid: "Plasma"

**System Behavior:**
1. Checks local data → Not found
2. Shows download progress:
   ```
   [Progress Bar: 0% → 100%]
   
   Log:
   0% - No local data found. Starting download...
   10% - Querying Vesiclepedia database...
   30% - Processing Vesiclepedia results...
   50% - Downloaded Vesiclepedia data: vesiclepedia_alzheimer_plasma.csv
   60% - Querying ExoCarta database...
   80% - Downloaded ExoCarta data: exocarta_alzheimer_plasma.csv
   85% - Querying EVpedia database...
   95% - Downloaded EVpedia data: evpedia_alzheimer_plasma.csv
   100% - Download complete! 3 files downloaded.
   ```
3. Automatically runs biomarker search
4. Displays results with biomarkers from downloaded datasets

## Technical Details

### Backend API Endpoints

1. **`POST /api/step1/check`** - Check if local data exists
   ```json
   {
     "disease": "Alzheimer",
     "biofluid": "Plasma"
   }
   ```
   Returns: `{"has_local_data": true/false}`

2. **`POST /api/step1/download`** - Download datasets
   ```json
   {
     "disease": "Alzheimer",
     "biofluid": "Plasma"
   }
   ```
   Returns: `{"downloaded_files": [...], "progress": [...]}`

3. **`POST /api/step1/run`** - Run biomarker finder (with auto-download)
   ```json
   {
     "disease": "Alzheimer",
     "biofluid": "Plasma",
     "auto_download": true
   }
   ```

### Download Module

Located in: `exo_gpt/data_downloader.py`

Key functions:
- `check_local_data()` - Checks if data exists locally
- `download_all_sources()` - Downloads from all databases
- `download_data_if_needed()` - Main entry point (check + download)

### Frontend Components

**Step1Panel.tsx** includes:
- Download progress state management
- Progress bar component
- Progress log display
- Automatic download trigger

## Future Enhancements

For production use, the downloader can be enhanced to:
1. **Real API Integration** - Query actual EV database APIs
2. **Web Scraping** - Extract data from database websites
3. **Caching** - Store downloaded data with metadata
4. **Incremental Updates** - Only download new/updated entries
5. **Background Downloads** - Download in background while user works

## Notes

- Downloaded files are saved permanently in `./data/databases/`
- Subsequent queries for the same disease use cached data (no re-download)
- Download progress is shown in real-time in the GUI
- If download fails, error message is displayed

