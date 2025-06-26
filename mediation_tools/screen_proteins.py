#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Parallel screening of Bonferroni-adjusted q-values and SuSiE PIP values.

A SNP–protein pair is reported as a “hit” when
    (q_value < bonferroni_threshold)  AND  (pip > pip_threshold)

Output:
    • One CSV file per (bonferroni_threshold, pip_threshold) pair
    • Columns: snp_name | protein | q_value | pip
"""

from __future__ import annotations

from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Iterable, List, Tuple

import pandas as pd

# ──────────────────────────────  configuration  ──────────────────────────────
BONFERRONI_PKL: Path = Path("q_snps_bonferroni.pkl")
PIP_PKL:        Path = Path("snps_with_pip_susie.pkl")

# Parameter grid
BONFERRONI_THRESHOLDS: tuple[float, ...] = (0.05, 0.01)
PIP_THRESHOLDS:        tuple[float, ...] = (0.90, 0.85, 0.80)

MAX_THREADS: int = 12     # parallel worker threads
# ──────────────────────────────────────────────────────────────────────────────


def load_and_align() -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Read both pickle files, harmonise dtypes, sort by (chr, rs),
    and return data frames with identical indices and column names.
    """
    bonf = pd.read_pickle(BONFERRONI_PKL)
    pip  = pd.read_pickle(PIP_PKL)

    # Basic header validation
    if bonf.columns.tolist()[:2] != ["chr", "rs"] \
       or pip.columns.tolist()[:2]  != ["chr", "rs"]:
        raise ValueError("The first two columns of both files must be "
                         "named 'chr' and 'rs'.")
    if bonf.columns.tolist()[2:] != pip.columns.tolist()[2:]:
        raise ValueError("Protein (column 3 onward) names differ.")

    # Harmonise dtypes
    bonf["chr"] = bonf["chr"].astype("int64")
    pip["chr"]  = pip["chr"].astype("int64")
    bonf["rs"]  = bonf["rs"].astype("string")
    pip["rs"]   = pip["rs"].astype("string")

    # Align on the intersection of (chr, rs)
    bonf = bonf.set_index(["chr", "rs"]).sort_index()
    pip  = pip.set_index(["chr", "rs"]).sort_index()

    common_idx = bonf.index.intersection(pip.index)
    if common_idx.empty:
        raise ValueError("No overlapping (chr, rs) pairs between the two "
                         "input files.")

    bonf = bonf.loc[common_idx].reset_index()
    pip  = pip.loc[common_idx].reset_index()
    return bonf, pip


def _collect_hits_for_protein(
    protein: str,
    bonf: pd.DataFrame,
    pip: pd.DataFrame,
    bonf_thr: float,
    pip_thr: float,
) -> List[Tuple[str, str, float, float]]:
    """
    Gather (rs, protein, q_value, pip) records for one protein column.
    Designed to be executed inside a ThreadPoolExecutor.
    """
    mask = (bonf[protein] < bonf_thr) & (pip[protein] > pip_thr)
    if not mask.any():
        return []

    rs          = bonf.loc[mask, "rs"]
    q_values    = bonf.loc[mask, protein]
    pip_values  = pip.loc[mask, protein]

    return [
        (snp, protein, float(qv), float(pv))
        for snp, qv, pv in zip(rs, q_values, pip_values, strict=True)
    ]


def screen_hits(
    bonf: pd.DataFrame,
    pip: pd.DataFrame,
    bonf_thr: float,
    pip_thr: float,
    *,
    n_threads: int = MAX_THREADS,
) -> pd.DataFrame:
    """
    Identify all SNP–protein hits for the given threshold pair.
    """
    protein_cols: Iterable[str] = bonf.columns[2:]   # skip chr, rs

    records: List[Tuple[str, str, float, float]] = []
    with ThreadPoolExecutor(max_workers=n_threads) as pool:
        futures = {
            pool.submit(_collect_hits_for_protein,
                        protein, bonf, pip, bonf_thr, pip_thr): protein
            for protein in protein_cols
        }
        for fut in as_completed(futures):
            records.extend(fut.result())

    return pd.DataFrame(records,
                        columns=["snp_name", "protein", "q_value", "pip"])


def export_csv(df: pd.DataFrame, bonf_thr: float, pip_thr: float) -> Path:
    """
    Write a result frame to disk and return its path.
    """
    fname = f"bonferroni_{bonf_thr:.3g}-pip_{pip_thr:.3g}.csv"
    path  = Path.cwd() / fname
    df.to_csv(path, index=False)
    return path


def main() -> None:
    # 1. Load both data sets once
    bonf_df, pip_df = load_and_align()

    # 2. Iterate over the 2 × 3 parameter grid
    for bonf_thr in BONFERRONI_THRESHOLDS:
        for pip_thr in PIP_THRESHOLDS:
            hits_df = screen_hits(
                bonf_df, pip_df,
                bonf_thr=bonf_thr,
                pip_thr=pip_thr,
                n_threads=MAX_THREADS,
            )
            csv_path = export_csv(hits_df, bonf_thr, pip_thr)

            print(
                f"[{bonf_thr:.3g}, {pip_thr:.3g}] "
                f"→ {len(hits_df):,} hits → '{csv_path.name}'"
            )


if __name__ == "__main__":
    main()