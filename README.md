# UVVisAnalysis

A ROOT script for automated UV-Vis spectrophotometer data analysis. This tool reads UV-Vis data files, averages repeated measurements, performs baseline correction, and fits the results to a Rayleigh scattering model — all without user input.

## Features
- Reads UV-Vis data files with minimal setup.
- Plots raw data for quick visualization.
- Automatically averages repeated measurements.
- Detects baseline files and performs baseline correction.
- Fits the corrected data with a Rayleigh scattering model.

## Requirements
- **ROOT CERN** installed (developed and tested with ROOT 6.32.10).
- Tested on Ubuntu Ubuntu 24.04 LTS (the code has several system calls, maybe test carfully for other distributions)

## File Structure and Naming Convention
The script relies on specific file structure and naming rules for automation:

- **Data Files:** Must follow this format:
  - **First 2 lines**: Header information (ignored by the script).
  - **Subsequent lines**: Numerical data in the format `wavelength,absorption\n`.
  - **Last line**: Follows the same numerical data format (no special markers required).

- **File Naming:**
  - **Measurements:** `XXXX_N.txt` where:
    - `XXXX` is the identifier (any desired name).
    - `N` is the repetition number (e.g., 1, 2, 3...). Always include `_N` even if there's only one measurement (_1 in this case).
  - **Baseline Files:**
    - Use tags `baseline_before` and `baseline_after` in the filename for baseline detection.
    - Baseline files may also include `_N` for repeated measurements. The baseline subtraction occurs **after averaging**.

## Usage
1. **Compile the Code:**
   ```bash
   root -l Analysis.C+
   ```
2. **Run the Analysis:**
   ```bash
   root -l 'Analysis.C'
   ```
3. **Output Files:**
   - `RawSpectrum.png` — Raw data visualization.
   - `Averages.png` — Averaged data with baseline curves.
   - `BaselineCorrected.png` — Final corrected data with Rayleigh fit.

## Key Functions Overview
- **`FancyPlot()`** — Enhances ROOT's visual settings for clearer plots.
- **`GetFileNames()`** — Collects `.txt` filenames in the working directory.
- **`SanitizeFile()`** — Cleans data by removing problematic lines.
- **`ReadFile_wavelength()` / `ReadFile_absorption()`** — Reads wavelength and absorption data from files.
- **`findDuplicates()`** — Identifies and averages repeated measurements.
- **`FitRayleigh()`** — Fits corrected data with a Rayleigh scattering model.

## Known Limitations
- Assumes file structure and naming conventions are strictly followed.
- Baseline correction assumes a smooth baseline behavior — unexpected fluctuations may reduce accuracy.

## Future Improvements
- Implement enhanced error handling for file I/O issues.
- Add user-defined parameter control for fitting ranges and smoothing options.

## Contributing
Contributions are welcome! Feel free to submit issues or pull requests to improve functionality, compatibility, or documentation.

## License
This project is licensed under the MIT License.

