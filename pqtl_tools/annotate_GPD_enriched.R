# -------------------------------------------------------------------------
# Script to read GPD_enriched.tsv and merge.bim, then annotate with ps, start, end
# -------------------------------------------------------------------------

# 1. Define file paths
base_dir   <- "C:/Users/shaok/Desktop/ProteoNexus_Citation/GRABBING"
gpd_file   <- file.path(base_dir, "GPD_enriched.tsv")
bim_file   <- file.path(base_dir, "merge.bim")
rdata_file <- file.path(base_dir, "gene_tbl37.RData")

# 2. Read input files
#    - GPD_enriched.tsv has a header, tab‐separated
#    - merge.bim has no header; the columns will be V1, V2, V3, V4, … by default
GPD_enriched <- read.delim(gpd_file,
                           header = TRUE,
                           stringsAsFactors = FALSE)

bim <- read.table(bim_file,
                  header = FALSE,
                  stringsAsFactors = FALSE)

# 3. Prepend columns "start", "end", "ps" (initialized to NA)
GPD_enriched <- cbind(
    start = NA_integer_,
    end   = NA_integer_,
    ps    = NA_real_,
    GPD_enriched
)

# 4. Annotate GPD_enriched$ps by matching rsID → bim$V2, pulling bim$V4
bim_idx <- match(GPD_enriched$rsID, bim$V2)
# bim_idx is NA where no match; bim$V4[bim_idx] will yield NA there
GPD_enriched$ps <- bim$V4[bim_idx]

# 5. Load gene_tbl37 (provides data.frame gene_tbl37 with columns gene_symbol, start, end)
load(rdata_file)
#    (Assumes gene_tbl37 has at least: gene_symbol, start, end)

# 6. Annotate start/end by matching protein → gene_tbl37$gene_symbol (case‐insensitive)
gene_idx <- match(
    tolower(GPD_enriched$protein),
    tolower(gene_tbl37$gene_symbol)
)
GPD_enriched$start <- gene_tbl37$start[gene_idx]
GPD_enriched$end   <- gene_tbl37$end[gene_idx]

# 7. Inspect the result
head(GPD_enriched)

# 8. (Optional) write out annotated table
write.table(
    GPD_enriched,
    file = file.path(base_dir, "GPD_enriched_annotated.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)