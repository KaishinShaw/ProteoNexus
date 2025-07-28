# ------------------------------------------------------------------------------
# Script:    annotate_gene_positions.R
# Purpose:   Read filtered summary CSVs, match on geneid1 to gencode annotation,
#            and append 'start' and 'end' genomic coordinates to each table.
# Author:    Hydraulik
# Date:      2025-07-28
# ------------------------------------------------------------------------------

# Ensure strings are not coerced to factors
options(stringsAsFactors = FALSE)

# ------------------------------------------------------------------------------
# 1. Configuration
# ------------------------------------------------------------------------------
# Base directory containing the input CSV files
base_dir <- "C:/Users/shaok/Desktop/ProteoNexus_Citation/GRABBING"

# List of input CSV filenames
csv_files <- c(
    "merged_sig_summ_all2_with_geneid1_filtered.csv",
    "merged_sig_summ_female2_with_geneid1_filtered.csv",
    "merged_sig_summ_male2_with_geneid1_filtered.csv"
)

# ------------------------------------------------------------------------------
# 2. Verify that 'gencode' exists in the environment and has required columns
# ------------------------------------------------------------------------------
if (!exists("gencode")) {
    stop("Cannot find object 'gencode' in the environment. Please load it before running this script.")
}

required_cols <- c("geneid1", "start", "end")
if (!all(required_cols %in% colnames(gencode))) {
    stop(
        "The 'gencode' object must contain the following columns:\n",
        paste(required_cols, collapse = ", ")
    )
}

# ------------------------------------------------------------------------------
# 3. Build lookup tables for fast, case‐insensitive matching on geneid1
# ------------------------------------------------------------------------------
# Normalize gencode$geneid1 to lower case
gencode_lc <- tolower(gencode$geneid1)

# Create named vectors mapping geneid1 -> start / end
start_lookup <- setNames(gencode$start, gencode_lc)
end_lookup   <- setNames(gencode$end,   gencode_lc)

# ------------------------------------------------------------------------------
# 4. Process each CSV: read, annotate, and write out annotated version
# ------------------------------------------------------------------------------
for (fname in csv_files) {
    # Construct full path and read the CSV
    infile <- file.path(base_dir, fname)
    message("Reading file: ", infile)
    df <- read.csv(infile, stringsAsFactors = FALSE)
    
    # Ensure the 'geneid1' column is present
    if (!"geneid1" %in% colnames(df)) {
        stop("Input file ", fname, " does not contain a 'geneid1' column.")
    }
    
    # Perform case‐insensitive lookup of start/end positions
    query_keys   <- tolower(df$geneid1)
    df$start     <- start_lookup[query_keys]
    df$end       <- end_lookup[query_keys]
    
    # Summary of annotation success
    n_total   <- nrow(df)
    n_matched <- sum(!is.na(df$start) & !is.na(df$end))
    message("  Annotated ", n_matched, " out of ", n_total, " rows.")
    
    # Write the annotated table to a new CSV
    out_name <- sub("\\.csv$", "_annotated.csv", fname)
    outfile  <- file.path(base_dir, out_name)
    message("Writing annotated data to: ", outfile)
    write.csv(df, outfile, row.names = FALSE)
}

message("All files have been processed successfully.")