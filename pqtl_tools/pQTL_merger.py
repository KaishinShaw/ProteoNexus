#!/usr/bin/env python3
"""
Multi-threaded merger for proteomics QTL analysis results.

This script merges significant SNP-protein associations from CSV files with
detailed association statistics from compressed text files, utilizing parallel
processing for optimal performance.
"""

import os
import gzip
import pandas as pd
import numpy as np
from pathlib import Path
from typing import List, Dict, Tuple, Optional
import multiprocessing as mp
from functools import partial
import logging
from datetime import datetime
import warnings

# Suppress pandas warnings for cleaner output
warnings.filterwarnings('ignore')

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)
logger = logging.getLogger(__name__)


class ProteoQTLMerger:
    """
    A class to handle the merging of proteomics QTL analysis results.
    
    Attributes:
        base_dir (Path): Base directory for SEARCH files
        pqtl_dir (Path): Base directory for pQTL output files
        n_cores (int): Number of CPU cores to use for parallel processing
    """
    
    def __init__(self, base_dir: str, pqtl_dir: str, n_cores: int = 10):
        """
        Initialize the ProteoQTLMerger.
        
        Args:
            base_dir: Path to the directory containing CSV files
            pqtl_dir: Path to the directory containing protein-specific outputs
            n_cores: Number of CPU cores to use (default: 10)
        """
        self.base_dir = Path(base_dir)
        self.pqtl_dir = Path(pqtl_dir)
        self.n_cores = min(n_cores, mp.cpu_count())
        
        # Define file mappings
        self.file_mappings = {
            'male': {
                'csv': 'male_bonferroni_0.05-pip_0.8.csv',
                'txt': 'summ_male2.assoc.txt.gz',
                'output': 'merged_sig_bon_male2.csv'
            },
            'female': {
                'csv': 'female_bonferroni_0.05-pip_0.8.csv',
                'txt': 'summ_female2.assoc.txt.gz',
                'output': 'merged_sig_bon_female2.csv'
            },
            'all': {
                'csv': 'all_bonferroni_0.05-pip_0.8.csv',
                'txt': 'summ_all2.assoc.txt.gz',
                'output': 'merged_sig_bon_all2.csv'
            }
        }
        
        logger.info(f"Initialized ProteoQTLMerger with {self.n_cores} cores")
    
    def read_compressed_file(self, filepath: Path) -> Optional[pd.DataFrame]:
        """
        Read a compressed tab-delimited file.
        
        Args:
            filepath: Path to the compressed file
            
        Returns:
            DataFrame containing the file contents, or None if error
        """
        try:
            with gzip.open(filepath, 'rt') as f:
                df = pd.read_csv(f, sep='\t', low_memory=False)
            return df
        except Exception as e:
            logger.error(f"Error reading {filepath}: {str(e)}")
            return None
    
    def process_protein_group(self, args: Tuple[str, pd.DataFrame, str]) -> pd.DataFrame:
        """
        Process a single protein group by merging CSV data with association statistics.
        
        Args:
            args: Tuple containing (protein_name, protein_df, txt_filename)
            
        Returns:
            Merged DataFrame for the protein
        """
        protein_name, protein_df, txt_filename = args
        
        # Construct path to the association file
        assoc_file_path = self.pqtl_dir / protein_name / 'output' / txt_filename
        
        if not assoc_file_path.exists():
            logger.warning(f"Association file not found: {assoc_file_path}")
            return protein_df
        
        # Read the association file
        assoc_df = self.read_compressed_file(assoc_file_path)
        
        if assoc_df is None or assoc_df.empty:
            logger.warning(f"Empty or invalid association file for protein: {protein_name}")
            return protein_df
        
        # Merge on SNP names
        merged_df = protein_df.merge(
            assoc_df,
            left_on='snp_name',
            right_on='rs',
            how='left'
        )
        
        # Remove duplicate 'rs' column as it's redundant with 'snp_name'
        if 'rs' in merged_df.columns:
            merged_df = merged_df.drop(columns=['rs'])
        
        return merged_df
    
    def process_dataset(self, dataset_type: str) -> None:
        """
        Process a complete dataset (male, female, or all).
        
        Args:
            dataset_type: Type of dataset to process ('male', 'female', or 'all')
        """
        logger.info(f"Processing {dataset_type} dataset...")
        
        # Get file paths
        csv_file = self.base_dir / self.file_mappings[dataset_type]['csv']
        txt_filename = self.file_mappings[dataset_type]['txt']
        output_file = self.base_dir / self.file_mappings[dataset_type]['output']
        
        # Read the CSV file
        if not csv_file.exists():
            logger.error(f"CSV file not found: {csv_file}")
            return
        
        logger.info(f"Reading CSV file: {csv_file}")
        csv_df = pd.read_csv(csv_file)
        
        # Get unique proteins
        unique_proteins = csv_df['protein'].unique()
        logger.info(f"Found {len(unique_proteins)} unique proteins")
        
        # Prepare arguments for parallel processing
        process_args = []
        for protein in unique_proteins:
            protein_df = csv_df[csv_df['protein'] == protein].copy()
            process_args.append((protein, protein_df, txt_filename))
        
        # Process proteins in parallel
        logger.info(f"Starting parallel processing with {self.n_cores} cores...")
        with mp.Pool(processes=self.n_cores) as pool:
            results = pool.map(self.process_protein_group, process_args)
        
        # Combine results
        logger.info("Combining results...")
        merged_df = pd.concat(results, ignore_index=True)
        
        # Sort by protein and SNP name for consistency
        merged_df = merged_df.sort_values(['protein', 'snp_name'])
        
        # Save the merged results
        logger.info(f"Saving merged results to: {output_file}")
        merged_df.to_csv(output_file, index=False)
        
        # Log summary statistics
        total_rows = len(merged_df)
        matched_rows = merged_df['chr'].notna().sum() if 'chr' in merged_df.columns else 0
        logger.info(f"Dataset {dataset_type} complete: {total_rows} total rows, "
                   f"{matched_rows} rows with association data")
    
    def run(self) -> None:
        """
        Execute the complete merging process for all datasets.
        """
        start_time = datetime.now()
        logger.info("Starting ProteoQTL merger...")
        
        # Process each dataset
        for dataset_type in ['male', 'female', 'all']:
            try:
                self.process_dataset(dataset_type)
            except Exception as e:
                logger.error(f"Error processing {dataset_type} dataset: {str(e)}")
                continue
        
        # Calculate and log execution time
        end_time = datetime.now()
        duration = end_time - start_time
        logger.info(f"Processing complete. Total execution time: {duration}")


def validate_environment() -> bool:
    """
    Validate that required directories exist.
    
    Returns:
        True if environment is valid, False otherwise
    """
    base_dir = Path('/gpfs/chencao/ysbioinfor/project/proteohubProject/web_out/SEARCH/')
    pqtl_dir = Path('/gpfs/chencao/ysbioinfor/project/proteohubProject/web_out/pqtl/')
    
    if not base_dir.exists():
        logger.error(f"Base directory does not exist: {base_dir}")
        return False
    
    if not pqtl_dir.exists():
        logger.error(f"pQTL directory does not exist: {pqtl_dir}")
        return False
    
    return True


def main():
    """
    Main execution function.
    """
    # Validate environment
    if not validate_environment():
        logger.error("Environment validation failed. Exiting.")
        return 1
    
    # Define paths
    base_dir = '/gpfs/chencao/ysbioinfor/project/proteohubProject/web_out/SEARCH/'
    pqtl_dir = '/gpfs/chencao/ysbioinfor/project/proteohubProject/web_out/pqtl/'
    
    # Initialize and run the merger
    merger = ProteoQTLMerger(
        base_dir=base_dir,
        pqtl_dir=pqtl_dir,
        n_cores=24
    )
    
    try:
        merger.run()
        logger.info("All processing completed successfully!")
        return 0
    except Exception as e:
        logger.error(f"Fatal error during processing: {str(e)}")
        return 1


if __name__ == "__main__":
    exit(main())
