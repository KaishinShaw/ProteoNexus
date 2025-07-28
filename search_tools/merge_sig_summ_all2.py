
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
merge_sig_summ_all2.py

This script searches for all `sig_summ_all2.assoc.tsv` files under the
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
    python merge_sig_summ_all2.py
    python merge_sig_summ_all2.py \
        --input-dir /gpfs/chencao/ysbioinfor/project/proteohubProject/web_out/pqtl \
        --output-csv merged_all_proteins.csv \
        --output-pkl merged_all_proteins.pkl

Dependencies:
    - Python 3.6+
    - pandas
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
        description="Merge multiple sig_summ_all2.assoc.tsv files into one DataFrame."
    )
    parser.add_argument(
        "--input-dir", "-i",
        default="/gpfs/chencao/ysbioinfor/project/proteohubProject/web_out/pqtl",
        help="Base directory containing per-protein subdirectories under 'pqtl/*/output/'."
    )
    parser.add_argument(
        "--output-csv", "-c",
        default="merged_sig_summ_all2.csv",
        help="Path to write the merged CSV file."
    )
    parser.add_argument(
        "--output-pkl", "-p",
        default="merged_sig_summ_all2.pkl",
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
    Locate all sig_summ_all2.assoc.tsv files under the given directory.

    Parameters
    ----------
    input_dir : str
        Base directory containing subdirectories of the form '*/output/'.

    Returns
    -------
    list of str
        Full paths to all matching TSV files.
    """
    pattern = os.path.join(input_dir, "*", "output", "sig_summ_all2.assoc.tsv")
    files = glob.glob(pattern)
    return sorted(files)


def load_and_annotate(file_path: str, column_name: str) -> pd.DataFrame:
    """
    Read a TSV file into a DataFrame and prepend a column with the protein name.

    Parameters
    ----------
    file_path : str
        Path to the sig_summ_all2.assoc.tsv file.
    column_name : str
        Column name to use for the protein identifier.

    Returns
    -------
    pd.DataFrame
        The annotated DataFrame.
    """
    # Infer the protein name from the parent directory of 'output'
    protein_name = Path(file_path).parent.parent.name

    # Read the TSV into a DataFrame
    df = pd.read_csv(file_path, sep="\t", header=0)

    # Insert the protein column at position 0
    df.insert(0, column_name, protein_name)
    return df


def main():
    args = parse_arguments()

    # Configure logging
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S"
    )

    logging.info(f"Searching for '*.assoc.tsv' files under: {args.input_dir}")
    assoc_files = find_assoc_files(args.input_dir)

    if not assoc_files:
        logging.error("No 'sig_summ_all2.assoc.tsv' files found. Exiting.")
        sys.exit(1)

    logging.info(f"Found {len(assoc_files)} files. Beginning to load and annotate...")

    dataframes = []
    for fp in assoc_files:
        try:
            df = load_and_annotate(fp, column_name=args.column_name)
            dataframes.append(df)
            logging.debug(f"Loaded and annotated: {fp} (rows: {len(df)})")
        except Exception as exc:
            logging.exception(f"Failed to process file '{fp}': {exc}")
            sys.exit(1)

    logging.info("Concatenating all DataFrames...")
    merged_df = pd.concat(dataframes, ignore_index=True)

    logging.info(f"Writing merged CSV to: {args.output_csv}")
    merged_df.to_csv(args.output_csv, index=False)

    logging.info(f"Writing merged pickle to: {args.output_pkl}")
    merged_df.to_pickle(args.output_pkl)

    logging.info("Done.")


if __name__ == "__main__":
    main()