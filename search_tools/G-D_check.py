#!/usr/bin/env python3
"""
Check presence and exclusivity of expected files for each rsID directory
based on a manifest CSV, and report any missing or extra directories/files.
"""

import csv
from pathlib import Path
import sys


def load_manifest(manifest_path):
    """
    Load the manifest CSV and return a list of dicts with keys: rsID, sex.
    """
    entries = []
    with manifest_path.open(newline='') as f:
        reader = csv.DictReader(f)
        for row in reader:
            rsid = row.get('rsID') or row.get('rsid')
            sex = row.get('sex')
            if not rsid or not sex:
                continue
            entries.append({'rsid': rsid, 'sex': sex})
    return entries


def scan_directories(base_dir):
    """
    Return a set of directory names under base_dir.
    """
    if not base_dir.exists() or not base_dir.is_dir():
        return set()
    return {p.name for p in base_dir.iterdir() if p.is_dir()}


def check_files_for_rsid(base_dir, rsid, sex):
    """
    For a given rsID directory and sex, return two sets:
    - missing_files: expected but not found
    - extra_files: found but not expected
    """
    expected = {
        f"effect_size_{sex}.tsv",
        f"top_trait_{sex}.png"
    }
    dir_path = base_dir / rsid
    if not dir_path.exists() or not dir_path.is_dir():
        # We handle missing directories at a higher level
        return expected.copy(), set()

    actual = {p.name for p in dir_path.iterdir() if p.is_file()}
    missing = expected - actual
    extra = actual - expected
    return missing, extra


def main():
    manifest_path = Path("/gpfs/chencao/ysbioinfor/project/proteohubProject/web_out/SEARCH/gpd_manifest.csv")
    base_dir = Path("/gpfs/chencao/ysbioinfor/project/proteohubProject/web_out/G-D")

    # Load manifest
    manifest = load_manifest(manifest_path)
    if not manifest:
        print(f"Error: No valid entries found in manifest {manifest_path}", file=sys.stderr)
        sys.exit(1)

    expected_dirs = {entry['rsid'] for entry in manifest}
    actual_dirs = scan_directories(base_dir)

    # Determine missing and extra directories
    missing_dirs = expected_dirs - actual_dirs
    extra_dirs = actual_dirs - expected_dirs

    # Check files within each expected directory
    missing_files = {}
    extra_files = {}
    for entry in manifest:
        rsid = entry['rsid']
        sex = entry['sex']
        if rsid not in actual_dirs:
            continue
        missing, extra = check_files_for_rsid(base_dir, rsid, sex)
        if missing:
            missing_files[rsid] = missing
        if extra:
            extra_files[rsid] = extra

    # Report
    print("\n=== Missing Directories ===")
    if missing_dirs:
        for d in sorted(missing_dirs):
            print(d)
    else:
        print("None")

    print("\n=== Extra Directories ===")
    if extra_dirs:
        for d in sorted(extra_dirs):
            print(d)
    else:
        print("None")

    print("\n=== Missing Files per rsID ===")
    if missing_files:
        for rsid, files in sorted(missing_files.items()):
            print(f"{rsid}: {', '.join(sorted(files))}")
    else:
        print("None")

    print("\n=== Extra Files per rsID ===")
    if extra_files:
        for rsid, files in sorted(extra_files.items()):
            print(f"{rsid}: {', '.join(sorted(files))}")
    else:
        print("None")


if __name__ == "__main__":
    main()