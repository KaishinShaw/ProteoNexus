#!/usr/bin/env python3
"""
Validate presence and exclusivity of expected files for each rsID directory
based on a manifest CSV, and report any missing or extra directories/files.

Expected files under /.../web_out/G-P/{rsID}/:
  - effect_size_{sex}.tsv
  - volcano_plot_{sex}.png
"""

import csv
import sys
from pathlib import Path
from typing import Dict, List, Set, Tuple


def load_manifest(manifest_path: Path) -> List[Dict[str, str]]:
    """
    Read the manifest CSV and return a list of entries with keys 'rsid' and 'sex'.
    """
    if not manifest_path.is_file():
        raise FileNotFoundError(f"Manifest file not found: {manifest_path}")

    entries: List[Dict[str, str]] = []
    with manifest_path.open(newline='', encoding='utf-8') as f:
        reader = csv.DictReader(f, delimiter=',')
        for row in reader:
            rsid = row.get('rsID') or row.get('rsid')
            sex = row.get('sex')
            if not rsid or not sex:
                # Skip rows missing required information
                continue
            entries.append({'rsid': rsid.strip(), 'sex': sex.strip()})
    return entries


def scan_directories(base_dir: Path) -> Set[str]:
    """
    Return the set of subdirectory names under base_dir.
    """
    if not base_dir.exists() or not base_dir.is_dir():
        return set()
    return {p.name for p in base_dir.iterdir() if p.is_dir()}


def check_files_for_rsid(base_dir: Path, rsid: str, sex: str) -> Tuple[Set[str], Set[str]]:
    """
    For a given rsID directory and sex, return two sets:
      - missing_files: expected but not found
      - extra_files: found but not expected
    """
    expected = {
        f"effect_size_{sex}.tsv",
        f"volcano_plot_{sex}.png",
    }
    dir_path = base_dir / rsid
    if not dir_path.is_dir():
        # Directory itself is missing; treat all expected as missing
        return expected.copy(), set()

    actual_files = {p.name for p in dir_path.iterdir() if p.is_file()}
    missing = expected - actual_files
    extra = actual_files - expected
    return missing, extra


def main():
    manifest_path = Path("/gpfs/chencao/ysbioinfor/project/proteohubProject/web_out/SEARCH/gpd_manifest.csv")
    base_dir = Path("/gpfs/chencao/ysbioinfor/project/proteohubProject/web_out/G-P")

    try:
        manifest = load_manifest(manifest_path)
    except FileNotFoundError as e:
        print(f"ERROR: {e}", file=sys.stderr)
        sys.exit(1)

    if not manifest:
        print(f"ERROR: No valid entries found in manifest {manifest_path}", file=sys.stderr)
        sys.exit(1)

    # Collect expected vs actual directories
    expected_dirs = {entry['rsid'] for entry in manifest}
    actual_dirs = scan_directories(base_dir)

    missing_dirs = expected_dirs - actual_dirs
    extra_dirs = actual_dirs - expected_dirs

    # Check files inside each expected directory
    missing_files: Dict[str, Set[str]] = {}
    extra_files: Dict[str, Set[str]] = {}

    for entry in manifest:
        rsid = entry['rsid']
        sex = entry['sex']
        missing, extra = check_files_for_rsid(base_dir, rsid, sex)
        if missing:
            missing_files[rsid] = missing
        if extra:
            extra_files[rsid] = extra

    # Reporting
    def print_section(title: str, items: List[str]):
        print(f"\n=== {title} ===")
        if items:
            for line in items:
                print(line)
        else:
            print("None")

    # Format directory reports
    print_section(
        "Missing Directories",
        sorted(missing_dirs)
    )
    print_section(
        "Extra Directories",
        sorted(extra_dirs)
    )

    # Format file reports per rsID
    missing_files_list = [
        f"{rsid}: missing {', '.join(sorted(files))}"
        for rsid, files in sorted(missing_files.items())
    ]
    extra_files_list = [
        f"{rsid}: extra {', '.join(sorted(files))}"
        for rsid, files in sorted(extra_files.items())
    ]

    print_section("Missing Files per rsID", missing_files_list)
    print_section("Extra Files per rsID", extra_files_list)


if __name__ == "__main__":
    main()