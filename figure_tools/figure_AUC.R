# figure_AUC.R
###############################################################################
# Section 1.  pheatmap workflow
###############################################################################

# ---- Install and load the required packages ---------------------------------
if (!requireNamespace("pheatmap", quietly = TRUE)) install.packages("pheatmap")
if (!requireNamespace("RColorBrewer", quietly = TRUE)) install.packages("RColorBrewer")

library(pheatmap)
library(RColorBrewer)

# ---- Load the data ----------------------------------------------------------
# Assumes a CSV file named `models.csv` in the working directory.
# The first column contains the row names: "all", "female", "male".
data <- read.csv(
    "models.csv",
    header      = TRUE,
    row.names   = 1,
    check.names = FALSE,
    na.strings  = c("", "NA")
)

# ---- Define a colour palette ------------------------------------------------
# Sequential, colour-blind-friendly palette (Yellow-Green-Blue).
my_palette <- colorRampPalette(brewer.pal(9, "YlGnBu"))(100)

# ---- Draw the heat-map and export as a high-resolution PDF ------------------
pdf("heatmap.pdf", width = 14, height = 4)   # Adjust dimensions as required
pheatmap(
    mat               = data,
    color             = my_palette,  # Colour palette
    cluster_rows      = FALSE,       # Disable row clustering
    cluster_cols      = FALSE,       # Disable column clustering
    show_rownames     = TRUE,        # Show row names on the left
    show_colnames     = FALSE,       # Hide column names
    na_col            = "grey90",    # Light grey for missing values
    border_color      = NA,          # No cell borders
    fontsize_row      = 10,          # Row-name font size
    cellwidth         = 12,          # Cell width  (px when printed)
    cellheight        = 12,          # Cell height (px when printed)
    legend            = TRUE,        # Display legend
    treeheight_row    = 0,           # Suppress row dendrogram
    treeheight_col    = 0            # Suppress column dendrogram
)
dev.off()

###############################################################################
# Section 2.  ComplexHeatmap workflow (more advanced & journal-grade)
###############################################################################

# ---- Install and load the required packages ---------------------------------
pkgs <- c("ComplexHeatmap", "circlize", "scico", "grid")
to_install <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(to_install)) install.packages(to_install)

library(ComplexHeatmap)
library(circlize)
library(scico)
library(grid)

# ---- Convert to a matrix ----------------------------------------------------
# Re-use the ‘data’ object created above.
mat <- as.matrix(data)

# ---- 1. Select an appropriate, SCI-compliant colour palette -----------------
has_pos <- any(mat > 0, na.rm = TRUE)
has_neg <- any(mat < 0, na.rm = TRUE)

if (has_pos && has_neg) {
    # Diverging palette: “vik” (blue–white–red, colour-blind-friendly).
    max_abs      <- max(abs(mat), na.rm = TRUE)
    palette      <- scico(255, palette = "vik")
    col_fun      <- colorRamp2(seq(-max_abs, max_abs, length.out = 255), palette)
    legend_at    <- round(seq(-max_abs, max_abs, length.out = 5), 2)
    legend_title <- "Value"
} else {
    # Sequential palette: “lajolla” (dark blue → light yellow, high contrast).
    palette      <- scico(255, palette = "lajolla", direction = -1)
    col_fun      <- colorRamp2(
        seq(min(mat, na.rm = TRUE), max(mat, na.rm = TRUE), length.out = 255),
        palette
    )
    legend_at    <- round(seq(min(mat, na.rm = TRUE),
                              max(mat, na.rm = TRUE), length.out = 5), 2)
    legend_title <- "AUC"
}

# ---- 2. Build the Heatmap object -------------------------------------------
ht <- Heatmap(
    mat,
    name                  = legend_title,
    col                   = col_fun,
    cluster_rows          = FALSE,
    cluster_columns       = FALSE,
    row_names_side        = "left",
    show_row_names        = TRUE,
    show_column_names     = FALSE,
    na_col                = "#EEEEEE",          # Light grey for missing values
    border                = NA,                 # No grid borders
    column_title          = "AUC Heatmap",
    column_title_side     = "top",
    column_title_gp       = gpar(fontsize = 16, fontface = "bold"),
    heatmap_legend_param  = list(
        title_gp            = gpar(fontsize = 12, fontface = "bold"),
        labels_gp           = gpar(fontsize = 10),
        at                  = legend_at,
        legend_height       = unit(40, "mm"),
        legend_graph_width  = unit(4,  "mm")
    )
)

# ---- 3. Export to a high-resolution PDF -------------------------------------
pdf("AUC_Heatmap.pdf",
    width  = 8,
    height = 3,
    useDingbats = FALSE)   # Prevent font-embedding issues in LaTeX workflows
draw(ht, heatmap_legend_side = "right")
dev.off()