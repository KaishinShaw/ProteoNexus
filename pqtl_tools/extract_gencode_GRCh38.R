###############################################################################
# Script:   extract_gencode_gene_table.R
# Purpose:  Import the GENCODE v48 “basic” GTF file and export a tidy
#           data.frame that contains one row per gene with commonly used
#           metadata fields.
# Author:   Hydraulik
# Date:     2025-07-28
# Requires: Bioconductor package “rtracklayer”
###############################################################################

library(rtracklayer)

## ---------------------------------------------------------------------------
## 1. Read GENCODE annotation
## ---------------------------------------------------------------------------
# Import the compressed GTF file.  The result is a GRanges object.
gtf <- import("gencode.v48.basic.annotation.gtf.gz")

## ---------------------------------------------------------------------------
## 2. Keep only gene-level features
## ---------------------------------------------------------------------------
genes <- gtf[gtf$type == "gene"]

## ---------------------------------------------------------------------------
## 3. Helper: return NA when the supplied vector is NULL
## ---------------------------------------------------------------------------
get_or_na <- function(x) {
  if (is.null(x)) NA else x
}

## ---------------------------------------------------------------------------
## 4. Build a tidy data.frame
## ---------------------------------------------------------------------------
gene_tbl <- data.frame(
  chr                 = as.character(seqnames(genes)),          # chromosome with “chr” prefix
  database            = genes$source,                           # annotation source
  start               = start(genes),                           # 1-based start coordinate
  end                 = end(genes),                             # 1-based end   coordinate
  description         = get_or_na(genes$gene_description),      # free-text description
  gene_id             = genes$gene_id,                          # Ensembl gene ID (with version)
  gene_type           = genes$gene_type,                        # Ensembl biotype
  gene_symbol         = genes$gene_name,                        # HGNC symbol
  gene_id_stripped    = sub("\\..*$", "", genes$gene_id),       # Ensembl ID without version
  chr_numeric         = sub("^chr", "", as.character(seqnames(genes))), # chromosome without “chr”
  stringsAsFactors    = FALSE
)

## ---------------------------------------------------------------------------
## 5. Quick look at the result
## ---------------------------------------------------------------------------
head(gene_tbl)