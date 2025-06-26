#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Efficient row-wise multiple-testing correction for very large SNP
p-value matrices.

Changes compared with the previous version
──────────────────────────────────────────
1. Saving the corrected matrix to HDF5 is *optional*  
   – pass a valid `output_hdf5` path or `None`.

2. Exporting significant associations to CSV is *optional*  
   – pass a valid `output_csv` path or `None`.

3. NEW: (optional) export the *entire* corrected matrix as a pandas
   DataFrame and pickle it.  
   • Pass a path via `output_pkl` (or `None` to skip).  
   • The pickle has exactly the same shape/column names as the
     original input; first two columns are `"chr"` (int32) and `"rs"`
     (str); the remaining columns contain the corrected q-values as
     float32.

Everything else – algorithms, memory footprint, multithreading – is
unchanged.
"""

from __future__ import annotations

import math
import pickle
import logging
from pathlib import Path
from typing import List, Tuple, Literal, Callable

import h5py
import numpy as np
import pandas as pd
from concurrent.futures import ThreadPoolExecutor, as_completed

# ──────────────────────────────────────────────────────────────────────────────
# logging
# ──────────────────────────────────────────────────────────────────────────────
logging.basicConfig(level=logging.INFO,
                    format="%(asctime)s │ %(levelname)s │ %(message)s")
logger = logging.getLogger(__name__)

# ──────────────────────────────────────────────────────────────────────────────
# I/O helpers
# ──────────────────────────────────────────────────────────────────────────────
def load_snp_dataframe(pkl_path: str | Path) -> pd.DataFrame:
    """
    Load a pickle containing a *pandas* DataFrame with at least the two
    mandatory columns ``"chr"`` and ``"rs"``.
    """
    pkl_path = Path(pkl_path)
    logger.info("Loading DataFrame from %s", pkl_path)
    with pkl_path.open("rb") as fh:
        df = pickle.load(fh)

    if not isinstance(df, pd.DataFrame):
        raise TypeError(f"Pickle must hold a DataFrame, found {type(df)}")
    if not {"chr", "rs"}.issubset(df.columns):
        raise ValueError("DataFrame must contain columns 'chr' and 'rs'")

    logger.info("Loaded DataFrame shape = %s", df.shape)
    return df


def extract_pvalue_matrix(df: pd.DataFrame
                          ) -> Tuple[pd.Series, pd.Series, np.ndarray, List[str]]:
    """
    Separate meta-data from numerical p-value matrix.
    """
    chr_col: pd.Series = df["chr"].astype(np.int32)
    rs_col: pd.Series = df["rs"].astype(str)
    pval_columns: List[str] = [c for c in df.columns if c not in ("chr", "rs")]

    p_values = df[pval_columns].to_numpy(dtype=np.float32, copy=False)
    n_missing = np.isnan(p_values).sum()
    logger.info("Missing values: %d (%.2f%%)",
                n_missing, 100 * n_missing / p_values.size)
    return chr_col, rs_col, p_values, pval_columns

# ──────────────────────────────────────────────────────────────────────────────
# Multiple-testing correction kernels
# ──────────────────────────────────────────────────────────────────────────────
def _bh_rowwise(p: np.ndarray) -> np.ndarray:
    m, n = p.shape
    valid = ~np.isnan(p)
    k = np.arange(1, n + 1, dtype=np.float32)
    cnt = valid.sum(axis=1, dtype=np.int32)

    order = np.argsort(p, axis=1, kind="mergesort")
    p_sorted = np.take_along_axis(p, order, axis=1)

    with np.errstate(divide="ignore", invalid="ignore"):
        q_tmp = (cnt[:, None] / k) * p_sorted

    q_tmp = np.where(np.isnan(q_tmp), np.inf, q_tmp)
    q_sorted = np.minimum.accumulate(q_tmp[:, ::-1], axis=1)[:, ::-1]
    np.clip(q_sorted, 0, 1, out=q_sorted)

    inv = np.argsort(order, axis=1)
    q = np.take_along_axis(q_sorted, inv, axis=1)
    q[~valid] = np.nan
    return q.astype(np.float32, copy=False)


def _bonferroni_rowwise(p: np.ndarray) -> np.ndarray:
    valid = ~np.isnan(p)
    m_valid = valid.sum(axis=1, keepdims=True, dtype=np.float32)
    with np.errstate(invalid="ignore"):
        q = p * m_valid
    np.clip(q, 0, 1, out=q)
    return q.astype(np.float32, copy=False)

# ──────────────────────────────────────────────────────────────────────────────
# Multithreaded driver
# ──────────────────────────────────────────────────────────────────────────────
def apply_correction_parallel(
        p_values: np.ndarray,
        method: Literal["bh", "bonferroni"] = "bh",
        n_threads: int = 10,
        chunk_rows: int | None = None) -> np.ndarray:

    if method not in {"bh", "bonferroni"}:
        raise ValueError("method must be 'bh' or 'bonferroni'")

    kernel: Callable[[np.ndarray], np.ndarray] = (
        _bh_rowwise if method == "bh" else _bonferroni_rowwise
    )

    m, n = p_values.shape
    if chunk_rows is None:
        # target chunk ~512 kB of float32
        chunk_rows = max(1, int(512_000 // (n * 4)))
    chunk_rows = min(chunk_rows, m)
    n_chunks = math.ceil(m / chunk_rows)
    logger.info("Parallel %s correction: %d × %d matrix → %d chunks × %d rows "
                "on %d threads", method.upper(), m, n, n_chunks, chunk_rows,
                n_threads)

    q_values = np.empty_like(p_values, dtype=np.float32)

    futures = {}
    with ThreadPoolExecutor(max_workers=n_threads) as exe:
        for idx in range(n_chunks):
            lo = idx * chunk_rows
            hi = min(m, lo + chunk_rows)
            view = p_values[lo:hi]          # shared read-only view
            fut = exe.submit(kernel, view)  # NumPy releases the GIL
            futures[fut] = (lo, hi)

        for fut in as_completed(futures):
            lo, hi = futures[fut]
            q_values[lo:hi] = fut.result()  # zero-copy assignment

    return q_values

# ──────────────────────────────────────────────────────────────────────────────
# Result export helpers
# ──────────────────────────────────────────────────────────────────────────────

def save_results_to_hdf5(
        output_path: str | Path,
        chr_col: pd.Series,
        rs_col: pd.Series,
        q_values: np.ndarray,
        pval_columns: List[str],
        method: str) -> None:
    """
    Store corrected q-values together with rich metadata in a compressed
    HDF5 file.
    """
    output_path = Path(output_path)
    logger.info("Writing HDF5 to %s", output_path)

    with h5py.File(output_path, "w") as hf:
        grp = hf.create_group("snp_data")

        # numeric chromosome vector
        grp.create_dataset("chr",
                           data=chr_col.to_numpy(np.int32),
                           compression="gzip", compression_opts=9)

        # rs IDs – use variable-length UTF-8 + supply *object* array
        str_t = h5py.string_dtype(encoding="utf-8")
        grp.create_dataset("rs",
                           data=rs_col.to_numpy(dtype=object),
                           dtype=str_t,
                           compression="gzip", compression_opts=9)

        # q-value matrix
        chunk_r = min(1000, q_values.shape[0])      # ~1 MB chunks
        chunk_c = min(100,  q_values.shape[1])
        grp.create_dataset("q_values",
                           data=q_values,
                           dtype="f4",
                           chunks=(chunk_r, chunk_c),
                           compression="gzip", compression_opts=9,
                           shuffle=True)

        # column names
        grp.create_dataset("column_names/pvalue_columns",
                           data=np.asarray(pval_columns, dtype=object),
                           dtype=str_t)

        # metadata
        hf.attrs.update({
            "description": "Row-wise multiple-testing correction for SNPs",
            "method": method,
            "created_by": "o3 SNP pipeline",
            "creation_date": pd.Timestamp.now().isoformat()
        })

    logger.info("HDF5 saved.")


def export_significant_results(
        rs_col: pd.Series,
        q_values: np.ndarray,
        pval_columns: List[str],
        threshold: float,
        output_csv: str | Path) -> None:
    """
    Flat file with every (SNP, test) pair that passes `threshold`.
    """
    logger.info("Exporting q < %.3g to %s", threshold, output_csv)
    sig_mask = q_values < threshold
    rows, cols = np.where(sig_mask)
    if rows.size == 0:
        logger.info("No significant associations found.")
        return

    df_out = pd.DataFrame({
        "snp_name": rs_col.iloc[rows].values,
        "pvalue_column": np.asarray(pval_columns)[cols],
        "q_value": q_values[rows, cols]
    }).sort_values("q_value")
    df_out.to_csv(output_csv, index=False)
    logger.info("Significant associations: %d (unique SNPs = %d, tests = %d)",
                len(df_out), df_out["snp_name"].nunique(),
                df_out["pvalue_column"].nunique())


def build_corrected_dataframe(chr_col: pd.Series,
                              rs_col: pd.Series,
                              q_values: np.ndarray,
                              pval_columns: List[str]) -> pd.DataFrame:
    """
    Assemble a DataFrame identical in shape/column order to the original
    input, but containing *corrected* q-values.
    """
    df_q = pd.DataFrame(q_values, columns=pval_columns, dtype=np.float32)
    df_q.insert(0, "rs",  rs_col.astype(str))
    df_q.insert(0, "chr", chr_col.astype(np.int32))
    return df_q

# ──────────────────────────────────────────────────────────────────────────────
# Main entry
# ──────────────────────────────────────────────────────────────────────────────
def run_pipeline(input_pkl: str | Path,
                 correction_method: Literal["bh", "bonferroni"] = "bh",
                 significance_threshold: float | None = None,
                 n_threads: int = 10,
                 *,
                 output_hdf5: str | Path | None = None,
                 output_csv: str | Path | None = None,
                 output_pkl: str | Path | None = None) -> None:
    """
    End-to-end driver.

    Parameters
    ----------
    input_pkl : str | Path
        Pickled DataFrame with p-values and columns ``"chr"``, ``"rs"``.
    correction_method : {"bh", "bonferroni"}
        Row-wise multiple-testing procedure.
    significance_threshold : float | None
        Threshold for calling associations significant.  If None, defaults
        to 0.05 for BH and 0.01 for Bonferroni.
    n_threads : int
        Number of worker threads.
    output_hdf5 / output_csv / output_pkl : path | None
        Pass a path to enable the corresponding export; use ``None`` to skip.
    """
    default_thresh = 0.05 if correction_method == "bh" else 0.01
    threshold = default_thresh if significance_threshold is None \
                else significance_threshold

    df = load_snp_dataframe(input_pkl)
    chr_col, rs_col, p_vals, pval_cols = extract_pvalue_matrix(df)

    q_vals = apply_correction_parallel(p_vals,
                                       method=correction_method,
                                       n_threads=n_threads)

    # --- Optional outputs ----------------------------------------------------
    if output_hdf5 is not None:
        save_results_to_hdf5(output_hdf5, chr_col, rs_col,
                             q_vals, pval_cols, correction_method)

    if output_csv is not None:
        export_significant_results(rs_col, q_vals, pval_cols,
                                   threshold=threshold,
                                   output_csv=output_csv)

    if output_pkl is not None:
        df_q = build_corrected_dataframe(chr_col, rs_col, q_vals, pval_cols)
        logger.info("Saving corrected matrix to pickle %s", output_pkl)
        with Path(output_pkl).open("wb") as fh:
            pickle.dump(df_q, fh, protocol=pickle.HIGHEST_PROTOCOL)
        logger.info("Pickle saved (shape = %s).", df_q.shape)

    logger.info("\n%s\nFinished %s correction – %d SNPs × %d tests\n"
                "Threshold = %.3g (threads = %d)\n%s",
                "═"*60, correction_method.upper(),
                p_vals.shape[0], p_vals.shape[1], threshold, n_threads,
                "═"*60)

# ──────────────────────────────────────────────────────────────────────────────
# CLI example
# ──────────────────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    run_pipeline(
        input_pkl="p_snps_cleaned.pkl",
        correction_method="bonferroni",          # "bh" or "bonferroni"
        significance_threshold=None,     # None → default per method
        n_threads=12,                    # ≤ physical cores
        # Optional outputs:
        output_hdf5="q_snps_bonferroni.h5",                # e.g. "q_snps_corrected.h5"
        output_csv=None,         # or None
        output_pkl="q_snps_bonferroni.pkl"  # or None
    )