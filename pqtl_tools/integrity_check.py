#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Check that each subdirectory under BASE_DIR contains the following six files
inside its 'plot' folder:
    - all_Manhattan.png
    - all_qqplot.png
    - female_Manhattan.png
    - female_qqplot.png
    - male_Manhattan.png
    - male_qqplot.png

If any files are missing, generate a CSV report named missing_files.csv
listing the subdirectory and the missing file(s).
"""

import csv
from pathlib import Path

# Root directory containing pqtl project subfolders
BASE_DIR = Path("/gpfs/chencao/ysbioinfor/project/proteohubProject/web_out/pqtl")

# List of expected files in each 'plot' directory
EXPECTED_FILES = [
    "all_Manhattan.png",
    "all_qqplot.png",
    "female_Manhattan.png",
    "female_qqplot.png",
    "male_Manhattan.png",
    "male_qqplot.png",
]

# Name of the output CSV file
OUTPUT_CSV = "missing_files.csv"

def main():
    missing_records = []

    # Iterate over each item in the base directory
    for subdir in BASE_DIR.iterdir():
        if not subdir.is_dir():
            # Skip if it's not a directory
            continue

        plot_dir = subdir / "plot"
        if not plot_dir.is_dir():
            # If the 'plot' directory is missing, mark all expected files as missing
            for fname in EXPECTED_FILES:
                missing_records.append({
                    "folder": subdir.name,
                    "missing_file": f"plot/ (directory missing), expected {fname}"
                })
            continue

        # Check each expected file
        for fname in EXPECTED_FILES:
            file_path = plot_dir / fname
            if not file_path.is_file():
                missing_records.append({
                    "folder": subdir.name,
                    "missing_file": fname
                })

    # Write the missing file records to CSV
    with open(OUTPUT_CSV, mode="w", newline="", encoding="utf-8") as csvfile:
        fieldnames = ["folder", "missing_file"]
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        for record in missing_records:
            writer.writerow(record)

    print(f"Check completed: found {len(missing_records)} missing file records.")
    print(f"See details in: {OUTPUT_CSV}")

if __name__ == "__main__":
    main()