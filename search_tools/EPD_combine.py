#!/usr/bin/env python3
"""
Concatenate all TSV files under /gpfs/chencao/ysbioinfor/project/proteohubProject/web_out/E-P-D/,
and prepend the following four columns (in this order):

 1. modules : fixed value "E-P-D"
 2. PHE     : the immediate parent directory name ({PHE number})
 3. PHD     : the leading part of the filename (everything before the last underscore)
 4. gender  : the trailing part of the filename (male / female / all)

Resulting in a single pandas.DataFrame. To save to disk, uncomment the last two lines.
"""

from pathlib import Path
import pandas as pd

# Configuration
ROOT_DIR    = Path("/gpfs/chencao/ysbioinfor/project/proteohubProject/web_out/E-P-D/")
TSV_PATTERN = "**/*.tsv"                # recursive glob pattern
NEW_COLS    = ["modules", "PHE", "PHD", "gender"]

data_frames = []

# Traverse all .tsv files under ROOT_DIR
for tsv_file in ROOT_DIR.glob(TSV_PATTERN):
    if not tsv_file.is_file():
        continue  # skip directories, hidden files, etc.

    # Extract metadata from path and filename
    phm_code = tsv_file.parent.name       # {PHM number}
    stem     = tsv_file.stem              # filename without .tsv
    parts    = stem.split("_")            # e.g. ["PHD123", "male"]

    if len(parts) < 2:
        raise ValueError(f"Filename does not conform to expected pattern: {tsv_file}")

    gender = parts[-1]                    # male / female / all
    phd    = "_".join(parts[:-1])         # allows underscores in the PHD code

    # Read the TSV into a DataFrame (as all strings to avoid type coercion)
    df = pd.read_csv(tsv_file, sep="\t", dtype=str)

    # Prepend our new columns
    df.insert(0, "modules", "E-P-D")
    df.insert(1, "PHE",      phm_code)
    df.insert(2, "PHD",      phd)
    df.insert(3, "gender",   gender)

    data_frames.append(df)

# Ensure we found at least one file
if not data_frames:
    raise RuntimeError(f"No TSV files found under: {ROOT_DIR}")

# Concatenate all the partial DataFrames
combined_df = pd.concat(data_frames, ignore_index=True)

# Reorder columns so that NEW_COLS appear on the far left,
# in case any TSV already had columns with those names.
combined_df = combined_df[NEW_COLS + [c for c in combined_df.columns if c not in NEW_COLS]]

# Example output
print(combined_df.head())
print(f"\nTotal rows: {len(combined_df):,}")

output_path = ROOT_DIR / "E-P-D_merged.tsv"
combined_df.to_csv('E-P-D_merged.tsv', sep="\t", index=False)