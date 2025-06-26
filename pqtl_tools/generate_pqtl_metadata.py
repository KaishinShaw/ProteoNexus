#!/usr/bin/env python3
"""
generate_pqtl_metadata.py

Scan the pQTL share directory tree and produce a JSON metadata file
describing which proteins (first‐level subdirectories) go into each
letter‐based tarball, and what output/plot files each protein contains.

Usage:
    python generate_pqtl_metadata.py \
        --root-dir /gpfs/chencao/ysbioinfor/project/proteohubProject/pqtl_share \
        [--output-file dataset_description.json]
"""

import os
import json
import argparse
from datetime import datetime
from collections import defaultdict

def scan_pqtl_share(root_dir):
    """
    Scan first‐level subdirectories under root_dir, group them by
    initial letter (a–z), and collect their output/plot files.
    Returns a dict: { letter: [ protein_info, ... ], ... }
    """
    groups = defaultdict(list)
    for entry in os.listdir(root_dir):
        path = os.path.join(root_dir, entry)
        if not os.path.isdir(path):
            continue
        if not entry:
            continue
        first = entry[0].lower()
        if first < 'a' or first > 'z':
            continue
        # Collect files under output/ and plot/
        info = {
            "name": entry,
            "output_files": [],
            "plot_files": []
        }
        out_dir = os.path.join(path, "output")
        if os.path.isdir(out_dir):
            info["output_files"] = sorted(os.listdir(out_dir))
        plot_dir = os.path.join(path, "plot")
        if os.path.isdir(plot_dir):
            info["plot_files"] = sorted(os.listdir(plot_dir))
        groups[first].append(info)
    # Sort each group's protein list by name
    for letter in groups:
        groups[letter].sort(key=lambda x: x["name"].lower())
    return groups

def build_metadata(groups, date_stamp):
    """
    Given the grouped proteins and date, build the final metadata dict.
    """
    metadata = {
        "date_generated": date_stamp,
        "tarballs": {}
    }
    for letter, proteins in sorted(groups.items()):
        tarball_name = f"ProteoNexus_pQTL_protein_{letter}_{date_stamp}.tar.gz"
        metadata["tarballs"][tarball_name] = {
            "letter": letter,
            "proteins": proteins
        }
    return metadata

def main():
    parser = argparse.ArgumentParser(
        description="Generate JSON metadata for pQTL protein tarballs."
    )
    parser.add_argument(
        "--root-dir",
        default="/gpfs/chencao/ysbioinfor/project/proteohubProject/pqtl_share",
        help="Root directory containing protein subdirectories."
    )
    parser.add_argument(
        "--output-file",
        default="dataset_description.json",
        help="Path to write the JSON metadata."
    )
    args = parser.parse_args()

    root = args.root_dir
    if not os.path.isdir(root):
        print(f"Error: root directory '{root}' does not exist.", file=sys.stderr)
        exit(1)

    # YYYYMMDD stamp
    today = datetime.now().strftime("%Y%m%d")

    groups = scan_pqtl_share(root)
    metadata = build_metadata(groups, today)

    # Write JSON
    out_path = args.output_file
    with open(out_path, "w", encoding="utf-8") as fh:
        json.dump(metadata, fh, indent=2, ensure_ascii=False)

    print(f"Metadata written to {out_path}")

if __name__ == "__main__":
    main()