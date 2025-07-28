#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
merge_sig_summ_female2.py

This script searches for all `sig_summ_female2.assoc.tsv` files under the
directory structure:

    /gpfs/chencao/ysbioinfor/project/proteohubProject/web_out/pqtl/*/output/

For each file found, it:
  1. Infers the protein name from the parent directory (the '*' wildcard).
  2. Reads the TSV into a pandas DataFrame.
  3. Inserts a column (by default named "protein_name") as the first column,
     populating it with the inferred protein name.
  4. Collects all annotated DataFrames and concatenates them vertically.

Finally, the script writes out the merged DataFrame as both:
 - A CSV file.
 - A pickle file.

Usage:
    python merge_sig_summ_female2.py
    python merge_sig_summ_female2.py \
        --input-dir /path/to/pqtl \
        --output-csv merged_female2.csv \
        --output-pkl merged_female2.pkl
"""

import os
import sys
import glob
import logging
import argparse
from pathlib import Path

import pandas as pd


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Merge multiple sig_summ_female2.assoc.tsv files into one DataFrame."
    )
    parser.add_argument(
        "--input-dir", "-i",
        default="/gpfs/chencao/ysbioinfor/project/proteohubProject/web_out/pqtl",
        help="Base directory containing per-protein subdirectories."
    )
    parser.add_argument(
        "--output-csv", "-c",
        default="merged_sig_summ_female2.csv",
        help="Path to write the merged CSV file."
    )
    parser.add_argument(
        "--output-pkl", "-p",
        default="merged_sig_summ_female2.pkl",
        help="Path to write the merged pickle file."
    )
    parser.add_argument(
        "--column-name",
        default="protein_name",
        help="Name of the new column for storing the protein name."
    )
    return parser.parse_args()


def find_assoc_files(input_dir: str) -> list:
    """
    Locate all sig_summ_female2.assoc.tsv files under the given directory.
    """
    pattern = os.path.join(input_dir, "*", "output", "sig_summ_female2.assoc.tsv")
    return sorted(glob.glob(pattern))


def load_and_annotate(file_path: str, column_name: str) -> pd.DataFrame:
    """
    Read a TSV file into a DataFrame and prepend a column with the protein name.
    """
    protein_name = Path(file_path).parent.parent.name
    df = pd.read_csv(file_path, sep="\t", header=0)
    df.insert(0, column_name, protein_name)
    return df


def main():
    args = parse_arguments()

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s"
    )

    logging.info(f"Searching for 'sig_summ_female2.assoc.tsv' under: {args.input_dir}")
    files = find_assoc_files(args.input_dir)
    if not files:
        logging.error("No files found. Exiting.")
        sys.exit(1)

    dfs = []
    for fp in files:
        try:
            dfs.append(load_and_annotate(fp, args.column_name))
            logging.debug(f"Loaded {fp}")
        except Exception as e:
            logging.exception(f"Failed to process {fp}: {e}")
            sys.exit(1)

    merged = pd.concat(dfs, ignore_index=True)
    logging.info(f"Writing merged CSV to {args.output_csv}")
    merged.to_csv(args.output_csv, index=False)
    logging.info(f"Writing merged pickle to {args.output_pkl}")
    merged.to_pickle(args.output_pkl)
    logging.info("Done.")


if __name__ == "__main__":
    main()