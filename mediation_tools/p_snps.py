import pandas as pd
import numpy as np
import gzip
import os

# Task 1: Read tab-separated file and create all_snps dataframe
print("Reading merge.bim file...")
all_snps = pd.read_csv('/gpfs/chencao/ysbioinfor/Datasets/ukb/geno/EUR_protein/hm3/all/merge.bim', 
                       sep='\t', 
                       header=None)

# Task 2: Keep only first two columns and rename them
all_snps = all_snps.iloc[:, :2]
all_snps.columns = ['chr', 'rs']

# Task 3: Read the file again (as specified in instructions)
# This seems redundant but following the instructions exactly
all_snps = pd.read_csv('/gpfs/chencao/ysbioinfor/Datasets/ukb/geno/EUR_protein/hm3/all/merge.bim', 
                       sep='\t', 
                       header=None)
all_snps = all_snps.iloc[:, :2]
all_snps.columns = ['chr', 'rs']

# Task 4: Create empty columns with specified names
phenotype_columns = [
    'PHD01001', 'PHD01004', 'PHD01007', 'PHD01010', 'PHD01013', 'PHD03001', 'PHD04002', 'PHD04005',
    'PHD05003', 'PHD06003', 'PHD08002', 'PHD08005', 'PHD09002', 'PHD10003', 'PHD10006', 'PHD10009',
    'PHD11002', 'PHD12001', 'PHD12004', 'PHD01002', 'PHD01005', 'PHD01008', 'PHD01011', 'PHD01014',
    'PHD03002', 'PHD04003', 'PHD05001', 'PHD06001', 'PHD07001', 'PHD08003', 'PHD08006', 'PHD10001',
    'PHD10004', 'PHD10007', 'PHD10010', 'PHD11003', 'PHD12002', 'PHD12005', 'PHD01003', 'PHD01006',
    'PHD01009', 'PHD01012', 'PHD02001', 'PHD04001', 'PHD04004', 'PHD05002', 'PHD06002', 'PHD08001',
    'PHD08004', 'PHD09001', 'PHD10002', 'PHD10005', 'PHD10008', 'PHD11001', 'PHD11004', 'PHD12003',
    'PHD12006'
]

# Initialize empty columns with NA (float)
for col in phenotype_columns:
    all_snps[col] = np.nan

# Task 5: Rename dataframe to p_snps
p_snps = all_snps

# Task 6: Iterate through each phenotype and fill p_wald values
print("Processing phenotype associations...")
for phenotype in phenotype_columns:
    file_path = f'/gpfs/chencao/ysbioinfor/project/proteohubProject/web_out/GWAS-D/{phenotype}/output/summ_all.assoc.txt.gz'
    
    # Check if file exists
    if not os.path.exists(file_path):
        print(f"Warning: File not found for {phenotype}: {file_path}")
        continue
    
    try:
        # Read compressed file
        with gzip.open(file_path, 'rt') as f:
            summ = pd.read_csv(f, sep='\t')
        
        # Keep only the second column (rs) and last column (p_wald)
        # Assuming the second column is at index 1 and last column is p_wald
        summ = summ.iloc[:, [1, -1]]
        summ.columns = ['rs', 'p_wald']
        
        # Create a mapping dictionary for efficient lookup
        rs_to_p_wald = dict(zip(summ['rs'], summ['p_wald']))
        
        # Fill p_wald values in p_snps based on matching rs values
        p_snps[phenotype] = p_snps['rs'].map(rs_to_p_wald)
        
        print(f"Successfully processed {phenotype}")
        
    except Exception as e:
        print(f"Error processing {phenotype}: {str(e)}")

# Task 7: Save p_snps as CSV file
output_path = 'p_snps_results.csv'
print(f"Saving results to {output_path}...")
p_snps.to_csv(output_path, index=False)
print("CSV file saved successfully!")

# Task 8: Drop columns with all NA values and save as pickle
print("\nCleaning dataset by removing columns with all NA values...")

# Identify columns with all NA values
na_columns = p_snps.columns[p_snps.isna().all()].tolist()

# Create cleaned dataframe by dropping columns with all NA values
p_snps_cleaned = p_snps.dropna(axis=1, how='all')

# Display cleaning summary
if na_columns:
    print(f"Removed {len(na_columns)} columns with all NA values:")
    for col in na_columns:
        if col in phenotype_columns:  # Only list phenotype columns
            print(f"  - {col}")
else:
    print("No columns with all NA values found.")

# Save cleaned dataframe as pickle file
pickle_output_path = 'p_snps_cleaned.pkl'
print(f"\nSaving cleaned dataset to {pickle_output_path}...")
p_snps_cleaned.to_pickle(pickle_output_path)
print("Pickle file saved successfully!")

# Display final summary statistics
print("\nFinal Summary Statistics:")
print(f"Original dataset dimensions: {p_snps.shape}")
print(f"Cleaned dataset dimensions: {p_snps_cleaned.shape}")
print(f"Number of SNPs: {len(p_snps_cleaned)}")
print(f"Number of phenotypes retained: {len([col for col in p_snps_cleaned.columns if col in phenotype_columns])}")

# Optional: Display retained phenotypes
retained_phenotypes = [col for col in p_snps_cleaned.columns if col in phenotype_columns]
print(f"\nRetained phenotypes ({len(retained_phenotypes)}):")
for phenotype in retained_phenotypes:
    non_na_count = p_snps_cleaned[phenotype].notna().sum()
    coverage_pct = (non_na_count / len(p_snps_cleaned)) * 100
    print(f"  {phenotype}: {non_na_count} SNPs ({coverage_pct:.2f}% coverage)")
