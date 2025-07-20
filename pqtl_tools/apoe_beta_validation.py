#!/usr/bin/env python3
"""
--------------------------------------------------------------------
Utility script for comparing discovery GWAS results with an existing
summary-statistics file.

1.  Concatenate all **extension-less** text files in the current
    working directory into a single DataFrame named `discovery`.

2.  Load `summ_all.assoc.txt`, keeping rows whose Wald p value
    (`p_wald`) is below a user-defined threshold (default: 0.05).

3.  Intersect the two tables on
        CHROM / GENPOS / ALLELE0 / ALLELE1
    and compute the Pearson correlation coefficient between
    `discovery.BETA` and `assoc.beta` for the overlapping SNPs.

Author: <Kaishin>
--------------------------------------------------------------------
"""

from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd
from scipy.stats import pearsonr

# ---------------------------------------------------------------------
# 1. Load and merge discovery files
# ---------------------------------------------------------------------
def load_discovery() -> pd.DataFrame:
    """
    Read every file in the current directory whose filename has
    *no* extension (e.g. ``sample`` rather than ``sample.txt``) and
    concatenate them row-wise.

    Returns
    -------
    pandas.DataFrame
        Tidy discovery table with standardised column names and
        numeric columns coerced to floats.
    """
    extensionless = [f for f in Path(".").iterdir()
                     if f.is_file() and f.suffix == ""]

    if not extensionless:
        raise FileNotFoundError(
            "No extension-less files were found in the current directory."
        )

    frames: list[pd.DataFrame] = []
    for fp in extensionless:
        df = pd.read_csv(
            fp,
            sep=r"\s+",
            header=None,   # raw files lack a header row
            skiprows=1,    # begin reading from the 2nd line
            engine="c"
        )
        frames.append(df)

    discovery = pd.concat(frames, ignore_index=True)

    discovery.columns = [
        "CHROM", "GENPOS", "ID",
        "ALLELE0", "ALLELE1",
        "A1FREQ", "INFO", "N", "TEST",
        "BETA", "SE", "CHISQ", "LOG10P", "EXTRA"
    ]

    # Force numeric columns; non-parsable values become NaN
    numeric_cols = ["GENPOS", "A1FREQ", "INFO", "N",
                    "BETA", "SE", "CHISQ", "LOG10P"]
    discovery[numeric_cols] = discovery[numeric_cols].apply(
        pd.to_numeric, errors="coerce"
    )

    return discovery


# ---------------------------------------------------------------------
# 2. Load the *assoc* file
# ---------------------------------------------------------------------
def load_assoc(path: str | Path, p_threshold: float = 1e-4) -> pd.DataFrame:
    """
    Load an association results file and retain rows with
    ``p_wald < p_threshold``.

    Parameters
    ----------
    path : str or pathlib.Path
        File system path to *summ_all.assoc.txt*.
    p_threshold : float, default 1e-4
        Wald-test p-value cutoff.

    Returns
    -------
    pandas.DataFrame
        Filtered association table.
    """
    assoc_cols = [
        "chr", "rs", "ps", "n_mis", "n_obs",
        "allele1", "allele0", "af",
        "beta", "se", "p_wald", "pip_susie"
    ]

    assoc = pd.read_csv(
        path,
        sep=r"\s+",
        header=0,        # first line already contains a header
        names=assoc_cols,
        engine="c"
    )

    assoc = assoc[assoc["p_wald"] < p_threshold].copy()
    assoc["beta"] = pd.to_numeric(assoc["beta"], errors="coerce")

    return assoc


# ---------------------------------------------------------------------
# 3 & 4. Merge tables and compute Pearson correlation
# ---------------------------------------------------------------------
def compute_correlation(
    discovery: pd.DataFrame,
    assoc: pd.DataFrame,
) -> None:
    """
    Merge `discovery` and `assoc` on genomic coordinates and alleles,
    then report the Pearson correlation between effect sizes.

    Raises
    ------
    ValueError
        If there is no overlap or fewer than two valid SNPs.
    """
    merged = discovery.merge(
        assoc,
        left_on=["CHROM", "GENPOS", "ALLELE0", "ALLELE1"],
        right_on=["chr",   "ps",     "allele0", "allele1"],
        how="inner",
        suffixes=("_disc", "_assoc")
    )

    if merged.empty:
        raise ValueError(
            "No overlap between discovery and assoc on the specified keys."
        )

    subset = merged[["BETA", "beta"]].dropna()
    if subset.shape[0] < 2:
        raise ValueError(
            "Not enough overlapping SNPs with valid effect sizes to "
            "compute a Pearson correlation."
        )

    r, p_val = pearsonr(subset["BETA"], subset["beta"])

    print(f"Number of overlapping SNPs : {subset.shape[0]}")
    print(f"Pearson r                  : {r: .5f}")
    print(f"p-value                    : {p_val: .3e}")


# ---------------------------------------------------------------------
# Main execution
# ---------------------------------------------------------------------
def main() -> None:
    discovery = load_discovery()

    assoc_path = (
        "/gpfs/chencao/ysbioinfor/project/proteohubProject/"
        "pqtl_share/apoe/output/summ_all.assoc.txt"
    )
    assoc = load_assoc(assoc_path)

    compute_correlation(discovery, assoc)


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:  # pylint: disable=broad-except
        sys.exit(f"Error: {exc}")