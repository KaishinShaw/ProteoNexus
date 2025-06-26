#!/usr/bin/env python3
"""
SNP–Phenotype Association Analysis Pipeline

This script reads SNP information from a PLINK .bim file and aggregates
SuSiE posterior inclusion probabilities (pip_susie) from compressed
association results for multiple phenotypes into a single DataFrame.
It then cleans the resulting matrix, saves it to a pickle file, and
generates summary statistics.
"""

import logging
import gzip
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
import warnings

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S"
)
logger = logging.getLogger(__name__)

# Suppress pandas performance warnings
warnings.filterwarnings("ignore", category=pd.errors.PerformanceWarning)


class SNPPhenotypeAnalyzer:
    """
    Aggregates pip_susie values for multiple phenotypes into a single matrix.
    """

    def __init__(self, base_dir: str, bim_file: str):
        """
        Initializes the analyzer with the base directory containing phenotype
        subdirectories and the path to the PLINK .bim file.

        Args:
            base_dir: Path to the directory containing phenotype subdirectories.
            bim_file: Path to the PLINK .bim file.
        """
        self.base_dir = Path(base_dir)
        self.bim_file = Path(bim_file)

    def load_snp_data(self) -> pd.DataFrame:
        """
        Loads SNP identifiers from the BIM file.

        Returns:
            DataFrame with columns ['chr', 'rs'].
        """
        logger.info(f"Loading BIM file: {self.bim_file}")
        try:
            snp_df = pd.read_csv(
                self.bim_file,
                sep="\t",
                header=None,
                names=["chr", "rs", "pos_cm", "pos_bp", "a1", "a2"],
                usecols=[0, 1]
            )
            snp_df.columns = ["chr", "rs"]
            logger.info(f"Loaded {len(snp_df):,} SNPs")
            return snp_df
        except Exception as e:
            logger.error(f"Failed to read BIM file: {e}")
            raise

    def list_phenotypes(self) -> List[str]:
        """
        Lists all phenotype subdirectories in the base directory.

        Returns:
            Sorted list of phenotype names.
        """
        try:
            phenotypes = [d.name for d in self.base_dir.iterdir() if d.is_dir()]
            phenotypes.sort()
            logger.info(f"Found {len(phenotypes)} phenotypes")
            return phenotypes
        except Exception as e:
            logger.error(f"Failed to list phenotypes: {e}")
            raise

    def initialize_matrix(self,
                          snp_df: pd.DataFrame,
                          phenotypes: List[str]) -> pd.DataFrame:
        """
        Creates an empty DataFrame indexed by SNPs, with one column per phenotype.

        Args:
            snp_df: DataFrame containing SNP 'chr' and 'rs'.
            phenotypes: List of phenotype column names.

        Returns:
            DataFrame with columns ['chr', 'rs'] + phenotypes, all phenotype
            values initialized to NaN.
        """
        logger.info("Initializing results matrix...")
        empty = pd.DataFrame(index=snp_df.index,
                             columns=phenotypes,
                             dtype=np.float64)
        matrix = pd.concat([snp_df, empty], axis=1)
        logger.info(f"Initialized matrix of shape {matrix.shape}")
        return matrix

    def load_pip_susie(self, phenotype: str) -> Optional[Dict[str, float]]:
        """
        Loads 'rs' -> 'pip_susie' mappings for a single phenotype.

        Args:
            phenotype: Name of the phenotype subdirectory.

        Returns:
            Dictionary mapping rs IDs to pip_susie values, or None on failure.
        """
        assoc_path = self.base_dir / phenotype / "output" / "summ_all2.assoc.txt.gz"
        if not assoc_path.exists():
            logger.warning(f"Missing file for phenotype '{phenotype}': {assoc_path}")
            return None

        try:
            with gzip.open(assoc_path, "rt") as f:
                header = f.readline().strip().split("\t")
                rs_idx = 1
                pip_idx = header.index("pip_susie") if "pip_susie" in header else -1

                if pip_idx < 0:
                    logger.error(f"'pip_susie' column not found in {phenotype} header")
                    return None

                # rewind and read only the needed columns
                f.seek(0)
                df = pd.read_csv(
                    f,
                    sep="\t",
                    usecols=[rs_idx, pip_idx],
                    dtype={rs_idx: str, pip_idx: np.float64}
                )

            df.columns = ["rs", "pip_susie"]
            mapping = dict(zip(df["rs"], df["pip_susie"]))
            logger.info(f"Loaded {len(mapping):,} pip_susie values for '{phenotype}'")
            return mapping

        except Exception as e:
            logger.error(f"Error reading associations for '{phenotype}': {e}")
            return None

    def aggregate(self,
                  matrix: pd.DataFrame,
                  phenotypes: List[str]) -> pd.DataFrame:
        """
        Fills the matrix with pip_susie values for all phenotypes.

        Args:
            matrix: DataFrame initialized by initialize_matrix().
            phenotypes: List of phenotype names.

        Returns:
            DataFrame with pip_susie values filled in.
        """
        logger.info("Aggregating pip_susie across all phenotypes...")
        successes, failures = [], []

        for idx, pheno in enumerate(phenotypes, start=1):
            logger.info(f"[{idx}/{len(phenotypes)}] Processing '{pheno}'")
            mapping = self.load_pip_susie(pheno)
            if mapping:
                matrix[pheno] = matrix["rs"].map(mapping)
                successes.append(pheno)
            else:
                failures.append(pheno)

        logger.info(f"Successfully processed {len(successes)} phenotypes")
        if failures:
            logger.warning(f"Failed to process {len(failures)} phenotypes: {failures}")

        return matrix

    def clean(self,
              matrix: pd.DataFrame) -> Tuple[pd.DataFrame, List[str]]:
        """
        Drops phenotype columns that are entirely NaN.

        Args:
            matrix: DataFrame after aggregation.

        Returns:
            A tuple of (cleaned DataFrame, list of dropped columns).
        """
        logger.info("Cleaning matrix: dropping all-NaN columns...")
        all_nan = matrix.isna().all()
        to_drop = matrix.columns[all_nan].tolist()
        cleaned = matrix.loc[:, ~all_nan]

        if to_drop:
            logger.info(f"Dropped {len(to_drop)} columns: {to_drop}")
        else:
            logger.info("No all-NaN columns to drop")

        return cleaned, to_drop

    def summarize(self,
                  cleaned: pd.DataFrame,
                  phenotypes: List[str]) -> None:
        """
        Logs summary statistics: dimensions, number of SNPs, and coverage per phenotype.

        Args:
            cleaned: Cleaned DataFrame.
            phenotypes: Original list of phenotype names.
        """
        logger.info(f"Final matrix shape: {cleaned.shape}")
        snp_count = len(cleaned)
        retained = [c for c in cleaned.columns if c in phenotypes]
        logger.info(f"Number of retained phenotypes: {len(retained)}")

        stats = []
        for pheno in retained:
            non_na = cleaned[pheno].notna().sum()
            coverage = non_na / snp_count * 100
            stats.append((pheno, non_na, coverage))

        stats.sort(key=lambda x: x[2], reverse=True)

        logger.info("Top 5 phenotypes by coverage:")
        for pheno, count, pct in stats[:5]:
            logger.info(f"  {pheno}: {count:,} SNPs ({pct:.2f}%)")

        if len(stats) > 10:
            logger.info("...")
            logger.info("Bottom 5 phenotypes by coverage:")
            for pheno, count, pct in stats[-5:]:
                logger.info(f"  {pheno}: {count:,} SNPs ({pct:.2f}%)")

    def run(self,
            output_pickle: str = "snps_with_pip_susie.pkl") -> pd.DataFrame:
        """
        Executes the full pipeline: load SNPs, list phenotypes, initialize matrix,
        aggregate pip_susie, clean, save, and summarize.

        Args:
            output_pickle: Filename for the pickled output.

        Returns:
            The cleaned DataFrame containing pip_susie values.
        """
        snp_df = self.load_snp_data()
        phenotypes = self.list_phenotypes()
        matrix = self.initialize_matrix(snp_df, phenotypes)
        matrix = self.aggregate(matrix, phenotypes)
        cleaned, dropped = self.clean(matrix)

        logger.info(f"Saving cleaned matrix to '{output_pickle}'")
        cleaned.to_pickle(output_pickle)
        logger.info("Save completed successfully")

        self.summarize(cleaned, phenotypes)
        return cleaned


def main():
    BASE_DIR = "/gpfs/chencao/ysbioinfor/project/proteohubProject/web_out/G-P/"
    BIM_FILE = "/gpfs/chencao/ysbioinfor/Datasets/ukb/geno/EUR_protein/hm3/all/merge.bim"
    OUTPUT_FILE = "snps_with_pip_susie.pkl"

    analyzer = SNPPhenotypeAnalyzer(BASE_DIR, BIM_FILE)
    analyzer.run(OUTPUT_FILE)
    logger.info("Pipeline execution completed.")


if __name__ == "__main__":
    main()