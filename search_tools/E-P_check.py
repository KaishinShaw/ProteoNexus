#!/usr/bin/env python3
"""
directory_file_completeness_checker.py

This script scans each first-level subdirectory under a given base directory,
verifies the presence of required output files for three categories
("all", "female", "male") whenever an effect size file exists,
and reports any missing files as a pandas DataFrame.
"""

import sys
from pathlib import Path
from typing import List, Dict

import pandas as pd


def find_missing_files(base_dir: Path) -> List[Dict[str, str]]:
    """
    Traverse first-level subdirectories of base_dir. For each subdirectory
    that contains an 'effect_size_{category}.tsv' file, check for the
    presence of associated enrichment and plot files. Record any missing files.

    Parameters:
        base_dir (Path): The root directory containing subject subdirectories.

    Returns:
        List[Dict[str, str]]: A list of records with keys 'directory',
                               'category', and 'missing_files'.
    """
    records: List[Dict[str, str]] = []

    for subdir in sorted(base_dir.iterdir()):
        if not subdir.is_dir():
            continue

        for category in ("all", "female", "male"):
            effect_file = subdir / f"effect_size_{category}.tsv"
            if not effect_file.exists():
                continue

            required_files = [
                subdir / f"enrichment_GO_{category}.tsv",
                subdir / f"enrichment_KEGG_{category}.tsv",
                subdir / f"enrichment_plot_{category}.png",
                subdir / f"volcano_plot_{category}.png",
            ]
            missing = [p.name for p in required_files if not p.exists()]

            if missing:
                records.append({
                    "directory": str(subdir),
                    "category": category,
                    "missing_files": ", ".join(missing),
                })

    return records


def main():
    base_dir = Path("./E-P")
    if not base_dir.is_dir():
        sys.exit(f"Error: Base directory '{base_dir}' does not exist or is not a directory.")

    missing_records = find_missing_files(base_dir)

    if missing_records:
        df = pd.DataFrame(missing_records)
        # Display the table in interactive environments:
        try:
            from IPython.display import display  # type: ignore
            display(df)
        except ImportError:
            print(df.to_string(index=False))
    else:
        print("All required files are present in every subdirectory.")


if __name__ == "__main__":
    main()