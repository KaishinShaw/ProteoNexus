#!/usr/bin/env Rscript

# -----------------------------------------------------------------------------
# Script:    update_and_filter_geneid1.R
# Purpose:   Read two CSVs, overwrite geneid1 for specific protein_name values,
#            then remove rows where geneid1 is NA or empty.
# Usage:     Place this script in the directory with the 3 CSVs and run:
#            Rscript update_and_filter_geneid1.R
# -----------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
})

# Mapping from protein_name (lowercase) to ENSG IDs
ensg_mapping <- c(
  bap18     = "ENSG00000258315",
  cert      = "ENSG00000113163",
  ddx58     = "ENSG00000107201",
  dusp13    = "ENSG00000079393",
  ebi3_il27 = "ENSG00000105246",
  fam172a   = "ENSG00000113391",
  gba       = "ENSG00000177628",
  hla_a     = "ENSG00000206503",
  hla_e     = "ENSG00000204592",
  leg1      = "ENSG00000100097",
  ment      = "ENSG00000143443",
  mylpf     = "ENSG00000180209",
  skiv2l    = "ENSG00000204351",
  tdgf1     = "ENSG00000241186",
  wars      = "ENSG00000140105"
)

# CSV files to process
files <- c(
  "merged_sig_summ_all2_with_geneid1.csv",
  "merged_sig_summ_female2_with_geneid1.csv",
  "merged_sig_summ_male2_with_geneid1.csv"
)

for (fname in files) {
  if (!file.exists(fname)) {
    warning("File not found, skipping: ", fname)
    next
  }
  
  # Read input, coercing all columns to character
  df <- read_csv(fname, col_types = cols(.default = "c"))
  
  df_updated <- df %>%
    # create lowercase version of protein_name for lookup
    mutate(protein_lc = str_to_lower(protein_name)) %>%
    # look up in mapping; if found, overwrite geneid1
    mutate(
      geneid1 = coalesce(
        ensg_mapping[protein_lc],
        geneid1
      )
    ) %>%
    # drop helper column
    select(-protein_lc) %>%
    # remove rows where geneid1 is NA or empty string
    filter(!is.na(geneid1) & geneid1 != "")
  
  # Construct output filename
  out_fname <- sub("\\.csv$", "_filtered.csv", fname)
  
  # Write the filtered data
  write_csv(df_updated, out_fname)
  message("Wrote filtered file: ", out_fname)
}