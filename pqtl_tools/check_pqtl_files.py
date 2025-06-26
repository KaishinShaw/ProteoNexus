#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Check the completeness of files under
    /data/part1/proteo_a001/web_out/pqtl/*/{plot,output}/
and write any missing files to missing_files.csv.
"""

import csv
import sys
from pathlib import Path


def find_missing_files(base_dir: Path) -> list[dict]:
    """
    Traverse each subdirectory of base_dir and check for expected files
    in 'plot' and 'output' subdirectories. Return a list of dicts
    recording any missing files.
    """
    plot_expected = [
        "all_Manhattan.png",
        "all_qqplot.png",
        "female_Manhattan.png",
        "female_qqplot.png",
        "male_Manhattan.png",
        "male_qqplot.png",
    ]
    output_expected = [
        "sig_summ_all2.assoc.tsv",
        "sig_summ_female2.assoc.tsv",
        "sig_summ_male2.assoc.tsv",
        "summ_all2.assoc.txt.gz",
        "summ_female2.assoc.txt.gz",
        "summ_male2.assoc.txt.gz",
    ]

    missing = []

    if not base_dir.is_dir():
        raise NotADirectoryError(f"{base_dir} does not exist or is not a directory")

    for study_dir in base_dir.iterdir():
        if not study_dir.is_dir():
            continue

        # Check the 'plot' subdirectory
        plot_dir = study_dir / "plot"
        for filename in plot_expected:
            if not (plot_dir / filename).is_file():
                missing.append({
                    "directory": str(plot_dir),
                    "missing_file": filename,
                })

        # Check the 'output' subdirectory
        output_dir = study_dir / "output"
        for filename in output_expected:
            if not (output_dir / filename).is_file():
                missing.append({
                    "directory": str(output_dir),
                    "missing_file": filename,
                })

    return missing


def write_missing_to_csv(missing: list[dict], csv_path: Path) -> None:
    """
    Write the list of missing-file records to a CSV file.
    """
    try:
        with csv_path.open("w", newline="", encoding="utf-8") as csvfile:
            fieldnames = ["directory", "missing_file"]
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
            writer.writeheader()
            for record in missing:
                writer.writerow(record)
    except Exception as e:
        print(f"ERROR: Failed to write CSV file at {csv_path}: {e}", file=sys.stderr)
        sys.exit(1)


def main():
    base_dir = Path("/data/part1/proteo_a001/web_out/pqtl")
    output_csv = Path("missing_files.csv")

    try:
        missing_files = find_missing_files(base_dir)
    except NotADirectoryError as e:
        print(f"ERROR: {e}", file=sys.stderr)
        sys.exit(1)

    write_missing_to_csv(missing_files, output_csv)
    print(f"Check complete: found {len(missing_files)} missing files.")
    print(f"Results written to: {output_csv.resolve()}")


if __name__ == "__main__":
    main()