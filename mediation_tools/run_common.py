#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This script reads two CSV files (trait and protein), identifies the common entries
in their second columns, filters each DataFrame to keep only those common entries,
and writes the results to CSV files in the /common directory. 

The trait output file is named:
    .../{trait_filename}_{protein_filename}_common.csv

The protein output file is named:
    .../{protein_filename}_{trait_filename}_common.csv
"""

import os
import pandas as pd

def main(trait_path: str, protein_path: str, output_dir: str = "/gpfs/chencao/ysbioinfor/project/proteohubProject/snp_FDR/common/all") -> None:
    # Ensure output directory exists
    os.makedirs(output_dir, exist_ok=True)

    # Extract base filenames without extensions
    trait_filename = os.path.splitext(os.path.basename(trait_path))[0]
    protein_filename = os.path.splitext(os.path.basename(protein_path))[0]

    # Construct output file paths
    trait_outfile = os.path.join(
        output_dir,
        f"{trait_filename}_{protein_filename}_trait.csv"
    )
    protein_outfile = os.path.join(
        output_dir,
        f"{protein_filename}_{trait_filename}_protein.csv"
    )

    # 1. Read CSV files (assumes header on first row)
    df_trait = pd.read_csv(trait_path, dtype=str, encoding="utf-8")
    df_protein = pd.read_csv(protein_path, dtype=str, encoding="utf-8")

    # 2. Identify the common entries in the first column (index 0)
    col_name = df_trait.columns[0]
    trait_values = set(df_trait[col_name])
    protein_values = set(df_protein[col_name])
    common_values = trait_values & protein_values

    # 3. Filter each DataFrame by the intersection
    df_trait_common = df_trait[df_trait[col_name].isin(common_values)]
    df_protein_common = df_protein[df_protein[col_name].isin(common_values)]

    # 4. Write filtered results to CSV (no row indices, include headers)
    df_trait_common.to_csv(trait_outfile, index=False, encoding="utf-8")
    df_protein_common.to_csv(protein_outfile, index=False, encoding="utf-8")

    # 5. Summary output
    print(f"Found {len(common_values)} common entries in column '{col_name}'.")
    print(f"Trait common file: {trait_outfile} ({len(df_trait_common)} rows including header).")
    print(f"Protein common file: {protein_outfile} ({len(df_protein_common)} rows including header).")


if __name__ == "__main__":
    # Example usage: adjust paths as needed
    trait_csv_path   = "/gpfs/chencao/ysbioinfor/project/proteohubProject/snp_FDR/trait/all/BH_005.csv"
    protein_csv_path = "/gpfs/chencao/ysbioinfor/project/proteohubProject/snp_FDR/protein/all/bonferroni_001.csv"
    main(trait_csv_path, protein_csv_path)