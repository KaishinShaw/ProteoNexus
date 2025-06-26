#!/usr/bin/env python3
"""
SNP-Phenotype Association Analysis Pipeline

This script processes SNP data and phenotype associations to create a comprehensive
p-value matrix for genetic analysis.
"""

import pandas as pd
import numpy as np
import gzip
import os
import logging
from typing import Dict, List, Optional, Tuple
from pathlib import Path
import warnings

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)
logger = logging.getLogger(__name__)

# Suppress specific pandas warnings
warnings.filterwarnings('ignore', category=pd.errors.PerformanceWarning)


class SNPPhenotypeAnalyzer:
    """Analyzes SNP-phenotype associations and creates p-value matrices."""
    
    def __init__(self, base_dir: str, bim_file: str):
        """
        Initialize the analyzer with data paths.
        
        Args:
            base_dir: Base directory containing phenotype subdirectories
            bim_file: Path to the BIM file containing SNP information
        """
        self.base_dir = Path(base_dir)
        self.bim_file = Path(bim_file)
        self.p_snps: Optional[pd.DataFrame] = None
        
    def load_snp_data(self) -> pd.DataFrame:
        """
        Load SNP data from BIM file.
        
        Returns:
            DataFrame containing chromosome and rs columns
        """
        logger.info(f"Reading BIM file: {self.bim_file}")
        
        try:
            # Read BIM file (PLINK format)
            snp_data = pd.read_csv(
                self.bim_file,
                sep='\t',
                header=None,
                names=['chr', 'rs', 'pos_cm', 'pos_bp', 'a1', 'a2'],
                usecols=[0, 1]  # Only read first two columns
            )
            snp_data.columns = ['chr', 'rs']
            
            logger.info(f"Loaded {len(snp_data):,} SNPs")
            return snp_data
            
        except Exception as e:
            logger.error(f"Failed to load BIM file: {e}")
            raise
    
    def get_phenotype_columns(self) -> List[str]:
        """
        Get sorted list of phenotype directories.
        
        Returns:
            Sorted list of phenotype names
        """
        try:
            phenotypes = [
                d.name for d in self.base_dir.iterdir()
                if d.is_dir()
            ]
            phenotypes.sort()
            logger.info(f"Found {len(phenotypes)} phenotype directories")
            return phenotypes
            
        except Exception as e:
            logger.error(f"Failed to list phenotype directories: {e}")
            raise
    
    def initialize_dataframe(self, snp_data: pd.DataFrame, phenotypes: List[str]) -> pd.DataFrame:
        """
        Initialize DataFrame with SNP data and empty phenotype columns.
        
        Args:
            snp_data: DataFrame with SNP information
            phenotypes: List of phenotype column names
            
        Returns:
            DataFrame with SNP data and initialized phenotype columns
        """
        logger.info("Initializing phenotype columns...")
        
        # Create empty DataFrame for phenotype columns to avoid fragmentation
        phenotype_df = pd.DataFrame(
            index=snp_data.index,
            columns=phenotypes,
            dtype=np.float64
        )
        
        # Concatenate SNP data with phenotype columns
        p_snps = pd.concat([snp_data, phenotype_df], axis=1)
        
        logger.info(f"Initialized DataFrame with shape: {p_snps.shape}")
        return p_snps
    
    def load_phenotype_associations(self, phenotype: str) -> Optional[Dict[str, float]]:
        """
        Load p-values for a specific phenotype.
        
        Args:
            phenotype: Phenotype name
            
        Returns:
            Dictionary mapping rs IDs to p-values, or None if loading fails
        """
        file_path = self.base_dir / phenotype / "output" / "summ_all.assoc.txt.gz"
        
        if not file_path.exists():
            logger.warning(f"File not found for {phenotype}: {file_path}")
            return None
        
        try:
            # Read compressed association file
            with gzip.open(file_path, 'rt') as f:
                # Read header to determine column positions
                header = f.readline().strip().split('\t')
                
                # Find rs and p_wald column indices
                try:
                    rs_idx = 1  # Usually second column
                    p_wald_idx = header.index('p_wald') if 'p_wald' in header else -1
                except ValueError:
                    p_wald_idx = -1  # Default to last column
                
                # Reset file pointer
                f.seek(0)
                
                # Read data efficiently
                summ = pd.read_csv(
                    f,
                    sep='\t',
                    usecols=[rs_idx, p_wald_idx],
                    dtype={rs_idx: str, p_wald_idx: np.float64}
                )
            
            # Ensure proper column names
            summ.columns = ['rs', 'p_wald']
            
            # Create mapping dictionary
            rs_to_p_wald = dict(zip(summ['rs'], summ['p_wald']))
            
            logger.info(f"Loaded {len(rs_to_p_wald):,} associations for {phenotype}")
            return rs_to_p_wald
            
        except Exception as e:
            logger.error(f"Error processing {phenotype}: {e}")
            return None
    
    def process_all_phenotypes(self, p_snps: pd.DataFrame, phenotypes: List[str]) -> pd.DataFrame:
        """
        Process all phenotype associations and fill p-values.
        
        Args:
            p_snps: DataFrame with SNP data and empty phenotype columns
            phenotypes: List of phenotype names
            
        Returns:
            DataFrame with filled p-values
        """
        logger.info("Processing phenotype associations...")
        
        successful_phenotypes = []
        failed_phenotypes = []
        
        for i, phenotype in enumerate(phenotypes, 1):
            logger.info(f"Processing {phenotype} ({i}/{len(phenotypes)})...")
            
            rs_to_p_wald = self.load_phenotype_associations(phenotype)
            
            if rs_to_p_wald:
                # Vectorized mapping for efficiency
                p_snps[phenotype] = p_snps['rs'].map(rs_to_p_wald)
                successful_phenotypes.append(phenotype)
            else:
                failed_phenotypes.append(phenotype)
        
        # Log summary
        logger.info(f"Successfully processed: {len(successful_phenotypes)} phenotypes")
        if failed_phenotypes:
            logger.warning(f"Failed to process: {len(failed_phenotypes)} phenotypes")
        
        return p_snps
    
    def clean_dataframe(self, p_snps: pd.DataFrame) -> Tuple[pd.DataFrame, List[str]]:
        """
        Remove columns with all NA values.
        
        Args:
            p_snps: DataFrame with p-values
            
        Returns:
            Tuple of (cleaned DataFrame, list of removed columns)
        """
        logger.info("Cleaning dataset by removing columns with all NA values...")
        
        # Identify columns with all NA values
        na_mask = p_snps.isna().all()
        na_columns = p_snps.columns[na_mask].tolist()
        
        # Drop columns with all NA values
        p_snps_cleaned = p_snps.loc[:, ~na_mask]
        
        if na_columns:
            logger.info(f"Removed {len(na_columns)} columns with all NA values")
        else:
            logger.info("No columns with all NA values found")
        
        return p_snps_cleaned, na_columns
    
    def generate_summary_statistics(self, p_snps_cleaned: pd.DataFrame, phenotypes: List[str]) -> None:
        """
        Generate and log summary statistics.
        
        Args:
            p_snps_cleaned: Cleaned DataFrame
            phenotypes: Original list of phenotypes
        """
        logger.info("\nGenerating summary statistics...")
        
        # Basic statistics
        logger.info(f"Dataset dimensions: {p_snps_cleaned.shape}")
        logger.info(f"Number of SNPs: {len(p_snps_cleaned):,}")
        
        # Phenotype statistics
        retained_phenotypes = [col for col in p_snps_cleaned.columns if col in phenotypes]
        logger.info(f"Number of phenotypes retained: {len(retained_phenotypes)}")
        
        # Coverage statistics
        logger.info("\nPhenotype coverage statistics:")
        coverage_stats = []
        
        for phenotype in retained_phenotypes:
            non_na_count = p_snps_cleaned[phenotype].notna().sum()
            coverage_pct = (non_na_count / len(p_snps_cleaned)) * 100
            coverage_stats.append({
                'phenotype': phenotype,
                'snp_count': non_na_count,
                'coverage_pct': coverage_pct
            })
        
        # Sort by coverage
        coverage_stats.sort(key=lambda x: x['coverage_pct'], reverse=True)
        
        # Display top and bottom 5
        logger.info("Top 5 phenotypes by coverage:")
        for stat in coverage_stats[:5]:
            logger.info(f"  {stat['phenotype']}: {stat['snp_count']:,} SNPs ({stat['coverage_pct']:.2f}%)")
        
        if len(coverage_stats) > 10:
            logger.info("...")
            logger.info("Bottom 5 phenotypes by coverage:")
            for stat in coverage_stats[-5:]:
                logger.info(f"  {stat['phenotype']}: {stat['snp_count']:,} SNPs ({stat['coverage_pct']:.2f}%)")
    
    def run(self, output_pickle: str = 'p_snps_cleaned.pkl') -> pd.DataFrame:
        """
        Run the complete analysis pipeline.
        
        Args:
            output_pickle: Path for output pickle file
            
        Returns:
            Cleaned DataFrame with p-values
        """
        # Load SNP data
        snp_data = self.load_snp_data()
        
        # Get phenotype list
        phenotypes = self.get_phenotype_columns()
        
        # Initialize DataFrame
        p_snps = self.initialize_dataframe(snp_data, phenotypes)
        
        # Process phenotypes
        p_snps = self.process_all_phenotypes(p_snps, phenotypes)
        
        # Clean DataFrame
        p_snps_cleaned, removed_columns = self.clean_dataframe(p_snps)
        
        # Save results
        logger.info(f"Saving cleaned dataset to {output_pickle}...")
        p_snps_cleaned.to_pickle(output_pickle)
        logger.info("Pickle file saved successfully!")
        
        # Generate summary
        self.generate_summary_statistics(p_snps_cleaned, phenotypes)
        
        return p_snps_cleaned


def main():
    """Main execution function."""
    # Configuration
    BASE_DIR = "/gpfs/chencao/ysbioinfor/project/proteohubProject/web_out/G-P/"
    BIM_FILE = "/gpfs/chencao/ysbioinfor/Datasets/ukb/geno/EUR_protein/hm3/all/merge.bim"
    OUTPUT_FILE = "p_snps_cleaned.pkl"
    
    # Create analyzer and run pipeline
    analyzer = SNPPhenotypeAnalyzer(BASE_DIR, BIM_FILE)
    p_snps_cleaned = analyzer.run(OUTPUT_FILE)
    
    logger.info("Analysis completed successfully!")
    return p_snps_cleaned


if __name__ == "__main__":
    main()
