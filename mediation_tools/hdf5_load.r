# Example: Loading and using the HDF5 file
with h5py.File('q_snps_fdr_corrected.h5', 'r') as hf:
    # Access data
    chr_data = hf['snp_data/chr'][:]
    rs_data = hf['snp_data/rs'][:]
    q_values = hf['snp_data/q_values'][:]
    
    # Access metadata
    col_names = hf['column_names/pvalue_columns'][:]
    stats = dict(hf['statistics'].attrs)
    
    # Reconstruct DataFrame if needed
    df_results = pd.DataFrame(q_values, columns=col_names)
    df_results.insert(0, 'chr', chr_data)
    df_results.insert(1, 'rs', rs_data)