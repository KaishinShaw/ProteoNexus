#!/usr/bin/env python3
"""
Script to verify the completeness and modification times of pQTL analysis files.
Checks all protein folders for required output and plot files.
"""

import os
import glob
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Tuple
import sys

# Define required files for each subdirectory
REQUIRED_FILES = {
    'output': [
        'cs_all2.rds',
        'cs_female2.rds',
        'cs_male2.rds',
        'sig_summ_all2.assoc.tsv',
        'sig_summ_female2.assoc.tsv',
        'sig_summ_male2.assoc.tsv',
        'summ_all2.assoc.txt.gz',
        'summ_female2.assoc.txt.gz',
        'summ_male2.assoc.txt.gz',
        'summ_all.assoc.txt.gz',
        'summ_female.assoc.txt.gz',
        'summ_male.assoc.txt.gz'
    ],
    'plot': [
        'all_Manhattan.png',
        'all_qqplot.png',
        'female_Manhattan.png',
        'female_qqplot.png',
        'male_Manhattan.png',
        'male_qqplot.png'
    ]
}

# Base directory path
BASE_DIR = '/gpfs/chencao/ysbioinfor/project/proteohubProject/web_out/pqtl'

# Cutoff date: September 1, 2025
CUTOFF_DATE = datetime(2025, 9, 1)


class FileValidator:
    """Class to validate protein analysis files"""
    
    def __init__(self, base_dir: str):
        """
        Initialize the validator with base directory
        
        Args:
            base_dir: Base directory containing protein folders
        """
        self.base_dir = Path(base_dir)
        self.missing_files = {}
        self.outdated_files = {}
        
    def get_protein_folders(self) -> List[str]:
        """
        Get all protein folders in the base directory
        
        Returns:
            List of protein folder names
        """
        try:
            protein_folders = [d.name for d in self.base_dir.iterdir() 
                             if d.is_dir()]
            return sorted(protein_folders)
        except Exception as e:
            print(f"Error accessing base directory: {e}")
            return []
    
    def check_file_modification_time(self, file_path: Path) -> bool:
        """
        Check if file modification time is after the cutoff date
        
        Args:
            file_path: Path to the file
            
        Returns:
            True if file is modified after cutoff date, False otherwise
        """
        try:
            mtime = datetime.fromtimestamp(file_path.stat().st_mtime)
            return mtime >= CUTOFF_DATE
        except Exception:
            return False
    
    def validate_protein_folder(self, protein_name: str) -> Tuple[List[str], List[str]]:
        """
        Validate files in a single protein folder
        
        Args:
            protein_name: Name of the protein (folder name)
            
        Returns:
            Tuple of (missing_files, outdated_files)
        """
        protein_path = self.base_dir / protein_name
        missing = []
        outdated = []
        
        for subdir, required_files in REQUIRED_FILES.items():
            subdir_path = protein_path / subdir
            
            for file_name in required_files:
                file_path = subdir_path / file_name
                
                if not file_path.exists():
                    missing.append(f"{subdir}/{file_name}")
                elif not self.check_file_modification_time(file_path):
                    mtime = datetime.fromtimestamp(file_path.stat().st_mtime)
                    outdated.append(f"{subdir}/{file_name} (modified: {mtime.strftime('%Y-%m-%d %H:%M:%S')})")
        
        return missing, outdated
    
    def validate_all_proteins(self):
        """
        Validate all protein folders and collect issues
        """
        protein_folders = self.get_protein_folders()
        
        if not protein_folders:
            print("No protein folders found in the base directory.")
            return
        
        print(f"Found {len(protein_folders)} protein folders to check...")
        print("-" * 80)
        
        for protein_name in protein_folders:
            missing, outdated = self.validate_protein_folder(protein_name)
            
            if missing:
                self.missing_files[protein_name] = missing
            if outdated:
                self.outdated_files[protein_name] = outdated
    
    def generate_report(self):
        """
        Generate and print the validation report
        """
        print("\n" + "=" * 80)
        print("VALIDATION REPORT")
        print("=" * 80)
        
        # Report missing files
        if self.missing_files:
            print("\n## MISSING FILES ##")
            print("-" * 40)
            for protein_name, files in sorted(self.missing_files.items()):
                print(f"\nProtein: {protein_name}")
                print(f"Directory: {self.base_dir}/{protein_name}/")
                print("Missing files:")
                for file_name in files:
                    print(f"  - {file_name}")
        else:
            print("\n## MISSING FILES ##")
            print("No missing files found.")
        
        # Report outdated files
        if self.outdated_files:
            print("\n## FILES WITH OUTDATED MODIFICATION TIME ##")
            print(f"(Files modified before {CUTOFF_DATE.strftime('%Y-%m-%d')})")
            print("-" * 40)
            for protein_name, files in sorted(self.outdated_files.items()):
                print(f"\nProtein: {protein_name}")
                print(f"Directory: {self.base_dir}/{protein_name}/")
                print("Outdated files:")
                for file_info in files:
                    print(f"  - {file_info}")
        else:
            print("\n## FILES WITH OUTDATED MODIFICATION TIME ##")
            print("All files have valid modification times.")
        
        # Summary statistics
        print("\n" + "=" * 80)
        print("SUMMARY STATISTICS")
        print("=" * 80)
        total_proteins = len(self.get_protein_folders())
        proteins_with_issues = len(set(list(self.missing_files.keys()) + 
                                      list(self.outdated_files.keys())))
        proteins_complete = total_proteins - proteins_with_issues
        
        print(f"Total protein folders checked: {total_proteins}")
        print(f"Protein folders with all files present and up-to-date: {proteins_complete}")
        print(f"Protein folders with issues: {proteins_with_issues}")
        
        if self.missing_files:
            total_missing = sum(len(files) for files in self.missing_files.values())
            print(f"Total missing files: {total_missing}")
        
        if self.outdated_files:
            total_outdated = sum(len(files) for files in self.outdated_files.values())
            print(f"Total outdated files: {total_outdated}")
    
    def export_issues_to_file(self, output_file: str = "validation_issues.txt"):
        """
        Export validation issues to a text file
        
        Args:
            output_file: Name of the output file
        """
        with open(output_file, 'w') as f:
            f.write(f"pQTL File Validation Report\n")
            f.write(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
            f.write(f"Base directory: {self.base_dir}\n")
            f.write(f"Cutoff date: {CUTOFF_DATE.strftime('%Y-%m-%d')}\n")
            f.write("=" * 80 + "\n\n")
            
            if self.missing_files:
                f.write("MISSING FILES\n")
                f.write("-" * 40 + "\n")
                for protein_name, files in sorted(self.missing_files.items()):
                    f.write(f"\nProtein: {protein_name}\n")
                    for file_name in files:
                        f.write(f"  {self.base_dir}/{protein_name}/{file_name}\n")
            
            if self.outdated_files:
                f.write("\n\nOUTDATED FILES\n")
                f.write("-" * 40 + "\n")
                for protein_name, files in sorted(self.outdated_files.items()):
                    f.write(f"\nProtein: {protein_name}\n")
                    for file_info in files:
                        f.write(f"  {self.base_dir}/{protein_name}/{file_info}\n")
        
        print(f"\nDetailed report exported to: {output_file}")


def main():
    """Main execution function"""
    print("pQTL File Validation Tool")
    print("=" * 80)
    print(f"Base directory: {BASE_DIR}")
    print(f"Cutoff date: {CUTOFF_DATE.strftime('%Y-%m-%d')}")
    print("=" * 80)
    
    # Check if base directory exists
    if not os.path.exists(BASE_DIR):
        print(f"Error: Base directory does not exist: {BASE_DIR}")
        sys.exit(1)
    
    # Create validator and run validation
    validator = FileValidator(BASE_DIR)
    validator.validate_all_proteins()
    validator.generate_report()
    
    # Export issues to file if any issues found
    if validator.missing_files or validator.outdated_files:
        validator.export_issues_to_file()
        print("\nValidation completed with issues found.")
        sys.exit(1)
    else:
        print("\nValidation completed successfully! All files are present and up-to-date.")
        sys.exit(0)


if __name__ == "__main__":
    main()