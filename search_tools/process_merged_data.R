#!/usr/bin/env Rscript

# process_proteohub_merged_data.R
# Author: Hydraulik
# Date:   2025-06-23
# Purpose:
#   1. Read merged TSV and supporting CSV metadata files.
#   2. Enrich each merged table with exposure, measurement, and disease annotations.
#   3. Reorder columns so that newly added metadata follows the key identifier columns.
#   4. Export the processed tables to TSV for archival.

# Load required packages
if (!requireNamespace("tidyverse", quietly = TRUE)) {
  install.packages("tidyverse")
}
library(tidyverse)

# 1) Read the primary merged tables (TSV)
exp_prot_disease <- read_tsv("exp_prot_disease_merged.tsv", col_types = cols())
gene_prot_disease <- read_tsv("gene_prot_disease_merged.tsv", col_types = cols())
meas_prot_disease <- read_tsv("meas_prot_disease_merged.tsv", col_types = cols())

# 2) Read the metadata tables (CSV)
protein_hub_vars <- read_csv("protein_hub_variables.csv", col_types = cols())
exposures        <- read_csv("exposures.csv",             col_types = cols())
measurements     <- read_csv("measurements.csv",          col_types = cols())

# 3) Standardize and rename columns in metadata tables
exposures <- exposures %>%
  rename(
    Exposure_ID            = PHID,
    Exposure_Category      = Category,
    Exposure_Reported_Trait = Reported.Trait
  )

measurements <- measurements %>%
  rename(
    Measurement_ID             = PFID,
    Measurement_Category       = Category,
    Measurement_Reported_Trait = Reported.Trait
  )

protein_hub_vars <- protein_hub_vars %>%
  rename(
    Disease_ID               = PHID,
    Disease_ICD10_Category   = ICD10.category,
    Disease_Reported_Trait   = Reported.Trait,
    Disease_ICD10_Code       = ICD10
  )

# 4) Enrich exp_prot_disease
exp_prot_disease <- exp_prot_disease %>%
  # Join exposure metadata on PHE
  left_join(
    exposures %>% select(Exposure_ID, Exposure_Category, Exposure_Reported_Trait),
    by = c("PHE" = "Exposure_ID")
  ) %>%
  # Join disease metadata on PHD
  left_join(
    protein_hub_vars %>% select(Disease_ID, Disease_ICD10_Category, Disease_Reported_Trait, Disease_ICD10_Code),
    by = c("PHD" = "Disease_ID")
  ) %>%
  # Reorder: place exposure fields after PHE and disease fields after PHD
  relocate(Exposure_Category, Exposure_Reported_Trait, .after = PHE) %>%
  relocate(Disease_ICD10_Category, Disease_Reported_Trait, Disease_ICD10_Code, .after = PHD)

# 5) Enrich gene_prot_disease
gene_prot_disease <- gene_prot_disease %>%
  # Join disease metadata on PHD
  left_join(
    protein_hub_vars %>% select(Disease_ID, Disease_ICD10_Category, Disease_Reported_Trait, Disease_ICD10_Code),
    by = c("PHD" = "Disease_ID")
  ) %>%
  # Reorder: place disease fields after PHD
  relocate(Disease_ICD10_Category, Disease_Reported_Trait, Disease_ICD10_Code, .after = PHD)

# 6) Enrich meas_prot_disease
meas_prot_disease <- meas_prot_disease %>%
  # Join measurement metadata on PHM
  left_join(
    measurements %>% select(Measurement_ID, Measurement_Category, Measurement_Reported_Trait),
    by = c("PHM" = "Measurement_ID")
  ) %>%
  # Join disease metadata on PHD
  left_join(
    protein_hub_vars %>% select(Disease_ID, Disease_ICD10_Category, Disease_Reported_Trait, Disease_ICD10_Code),
    by = c("PHD" = "Disease_ID")
  ) %>%
  # Reorder: place measurement fields after PHM and disease fields after PHD
  relocate(Measurement_Category, Measurement_Reported_Trait, .after = PHM) %>%
  relocate(Disease_ICD10_Category, Disease_Reported_Trait, Disease_ICD10_Code, .after = PHD)

# 7) Inspect the resulting column order
cat("Columns in exp_prot_disease:\n"); print(colnames(exp_prot_disease))
cat("\nColumns in gene_prot_disease:\n"); print(colnames(gene_prot_disease))
cat("\nColumns in meas_prot_disease:\n"); print(colnames(meas_prot_disease))

# 8) Write processed tables back to TSV (for archival)
write_tsv(exp_prot_disease,       "exp_prot_disease_processed.tsv")
write_tsv(gene_prot_disease,      "gene_prot_disease_processed.tsv")
write_tsv(meas_prot_disease,      "meas_prot_disease_processed.tsv")