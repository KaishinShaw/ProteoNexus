# figure_AUC_SCI_publication.R
################################################################################
# Publication-Quality AUC Heatmap Generation for Scientific Journals
# Author: Hydraulik
# Date: 2025-01-30
# Description: Generate high-quality heatmaps for AUC values without cell 
#              annotations, suitable for SCI journal submission
################################################################################

# --- Package Management -------------------------------------------------------
# Function to install and load packages
load_packages <- function(packages, bioc = FALSE) {
    for (pkg in packages) {
        if (!requireNamespace(pkg, quietly = TRUE)) {
            if (bioc) {
                if (!requireNamespace("BiocManager", quietly = TRUE)) {
                    install.packages("BiocManager")
                }
                BiocManager::install(pkg, update = FALSE, ask = FALSE)
            } else {
                install.packages(pkg, repos = "https://cran.r-project.org")
            }
        }
        library(pkg, character.only = TRUE)
    }
}

# Load required packages
cran_packages <- c("pheatmap", "RColorBrewer", "viridisLite", "Cairo", 
                   "grid", "scales")
bioc_packages <- c("ComplexHeatmap", "circlize")

load_packages(cran_packages)
load_packages(bioc_packages, bioc = TRUE)

# --- Font Setup ---------------------------------------------------------------
# Define font family with fallback
font_family <- "sans"  # Use default sans-serif font for compatibility

# --- Data Loading and Validation ----------------------------------------------
# Check if data file exists
data_file <- "models.csv"
if (!file.exists(data_file)) {
    stop(paste("Error: Data file", data_file, "not found in working directory"))
}

# Load data with error handling
tryCatch({
    data <- read.csv(
        data_file,
        header      = TRUE,
        row.names   = 1,
        check.names = FALSE,
        na.strings  = c("", "NA", "N/A", "n/a"),
        stringsAsFactors = FALSE
    )
}, error = function(e) {
    stop(paste("Error reading data file:", e$message))
})

# Validate data
if (nrow(data) == 0 || ncol(data) == 0) {
    stop("Error: Data file is empty")
}

# --- Data Preprocessing -------------------------------------------------------
# Define row mapping for publication (in the order you want them displayed)
row_map <- c(
    all    = "Sex-combined",
    male   = "Male", 
    female = "Female"
)

# Check if expected rows exist
expected_rows <- names(row_map)
missing_rows <- expected_rows[!expected_rows %in% rownames(data)]
if (length(missing_rows) > 0) {
    warning(paste("Missing rows:", paste(missing_rows, collapse = ", ")))
    # Use only available rows
    row_map <- row_map[names(row_map) %in% rownames(data)]
}

# Reorder and rename rows
if (length(row_map) > 0) {
    data_ordered <- data[names(row_map), , drop = FALSE]
    rownames(data_ordered) <- row_map
} else {
    stop("Error: No expected rows found in data")
}

# Convert to matrix
mat <- as.matrix(data_ordered)

# Check for all NA matrix
if (all(is.na(mat))) {
    stop("Error: All values in the matrix are NA")
}

# --- Define Professional Color Palettes ---------------------------------------
# Scientific color palettes suitable for SCI publications

# 1. Sequential palette for positive values (Nature-style blue gradient)
nature_blue <- colorRampPalette(c("#F7FBFF", "#DEEBF7", "#C6DBEF", 
                                  "#9ECAE1", "#6BAED6", "#4292C6", 
                                  "#2171B5", "#08519C", "#08306B"))(100)

# 2. Viridis palette (colorblind-friendly)
viridis_palette <- viridis(100, option = "D", begin = 0.05, end = 0.95)

# 3. Custom scientific palette (Cell/Science style)
cell_palette <- colorRampPalette(c("#FFF5F0", "#FEE0D2", "#FCBBA1", 
                                   "#FC9272", "#FB6A4A", "#EF3B2C", 
                                   "#CB181D", "#A50F15", "#67000D"))(100)

################################################################################
# SECTION 1: Simple Publication-Quality Heatmap with pheatmap
################################################################################

# Calculate breaks for legend
data_range <- range(mat, na.rm = TRUE)
breaks_seq <- seq(data_range[1], data_range[2], length.out = 101)
legend_breaks <- pretty(data_range, n = 5)

# --- Generate pheatmap with professional settings -----------------------------
# Calculate dimensions for optimal layout
n_cols <- ncol(mat)
n_rows <- nrow(mat)

# Calculate cell dimensions for good proportions
cellheight_base <- 20  # Base height for each cell
cellwidth_calc <- 35   # Width for each cell

# Define output dimensions (inches)
fig_width <- max(8, (n_cols * cellwidth_calc / 72) + 2.5)  # Add space for legend
fig_height <- max(3, (n_rows * cellheight_base / 72) + 1.5)  # Add space for margins

# Create high-resolution PDF
cairo_pdf(
    filename  = "Figure_1_AUC_Heatmap_pheatmap.pdf",
    width     = fig_width,
    height    = fig_height,
    family    = font_family,
    pointsize = 10
)

# Create heatmap with scientific styling
pheatmap(
    mat               = mat,
    color             = nature_blue,  # Professional blue gradient
    breaks            = breaks_seq,
    cluster_rows      = FALSE,        # No clustering to maintain order
    cluster_cols      = FALSE,
    show_rownames     = TRUE,         # Show row names on left
    show_colnames     = FALSE,        # Hide column names for cleaner look
    fontsize_row      = 12,           # Font size for row labels
    cellwidth         = cellwidth_calc,
    cellheight        = cellheight_base,
    legend            = TRUE,
    legend_breaks     = legend_breaks,
    legend_labels     = sprintf("%.2f", legend_breaks),
    border_color      = "white",      # White borders for clarity
    na_col            = "#E5E5E5",    # Light gray for NA values
    main              = "",            # No title (add in figure caption)
    gaps_row          = NULL,
    gaps_col          = NULL,
    display_numbers   = FALSE,         # No numbers displayed
    annotation_legend = TRUE,
    annotation_names_row = TRUE,
    annotation_names_col = FALSE
)

# Add legend title manually
grid.text("AUC", x = unit(0.85, "npc"), y = unit(0.95, "npc"), 
          gp = gpar(fontsize = 11, fontface = "bold"))

dev.off()

################################################################################
# SECTION 2: Advanced Publication-Quality Heatmap with ComplexHeatmap
################################################################################

# --- Determine appropriate color scheme ---------------------------------------
has_negative <- any(mat < 0, na.rm = TRUE)
has_positive <- any(mat > 0, na.rm = TRUE)

if (has_negative && has_positive) {
    # Use diverging palette centered at zero
    max_abs <- max(abs(data_range))
    col_fun <- colorRamp2(
        breaks = c(-max_abs, -max_abs/2, 0, max_abs/2, max_abs),
        colors = c("#053061", "#4393C3", "#F7F7F7", "#D6604D", "#67001F")
    )
    legend_title <- "AUC"
    legend_at <- seq(-max_abs, max_abs, length.out = 5)
} else {
    # Use sequential palette
    col_fun <- colorRamp2(
        breaks = seq(data_range[1], data_range[2], length.out = 100),
        colors = viridis_palette
    )
    legend_title <- "AUC"
    legend_at <- pretty(data_range, n = 5)
}

# --- Create ComplexHeatmap with publication settings --------------------------
# Calculate dimensions for optimal layout
heatmap_height <- unit(n_rows * 0.8, "cm")
heatmap_width <- unit(n_cols * 1.2, "cm")

# Create heatmap without cell values
ht <- Heatmap(
    matrix              = mat,
    name                = legend_title,
    col                 = col_fun,
    
    # Clustering settings
    cluster_rows        = FALSE,
    cluster_columns     = FALSE,
    show_row_dend       = FALSE,
    show_column_dend    = FALSE,
    
    # Row/column names
    show_row_names      = TRUE,       # Show row names
    show_column_names   = FALSE,      # Hide column names
    row_names_side      = "left",     # Row names on the left side
    
    # Text formatting
    row_names_gp        = gpar(fontsize = 12, fontfamily = font_family),
    
    # Cell formatting
    rect_gp             = gpar(col = "white", lwd = 0.5),
    na_col              = "#E5E5E5",
    
    # Legend settings - optimized for layout
    heatmap_legend_param = list(
        title           = legend_title,
        title_position  = "topleft",
        at              = legend_at,
        labels          = sprintf("%.2f", legend_at),
        title_gp        = gpar(fontsize = 11, fontface = "bold", 
                               fontfamily = font_family),
        labels_gp       = gpar(fontsize = 10, fontfamily = font_family),
        legend_height   = unit(min(40, n_rows * 10), "mm"),  # Scale with heatmap
        grid_width      = unit(4, "mm"),
        border          = "black"
    ),
    
    # Size settings
    width               = heatmap_width,
    height              = heatmap_height
)

# --- Export to multiple formats -----------------------------------------------
# Calculate figure dimensions
fig_width  <- max(8, (n_cols * 1.2 / 2.54) + 2)  # Convert cm to inches, add margin
fig_height <- max(3, (n_rows * 0.8 / 2.54) + 1)

# 1. PDF (vector format, best for publications)
cairo_pdf(
    filename  = "Figure_2_AUC_Heatmap_ComplexHeatmap.pdf",
    width     = fig_width,
    height    = fig_height,
    family    = font_family,
    pointsize = 10
)
draw(ht, 
     heatmap_legend_side = "right",
     padding = unit(c(5, 10, 5, 5), "mm"),
     merge_legend = TRUE)
dev.off()

# 2. TIFF (600 DPI for print quality)
tiff(
    filename    = "Figure_2_AUC_Heatmap_ComplexHeatmap.tiff",
    width       = fig_width,
    height      = fig_height,
    units       = "in",
    res         = 600,  # High resolution for print
    compression = "lzw"
)
draw(ht, 
     heatmap_legend_side = "right",
     padding = unit(c(5, 10, 5, 5), "mm"),
     merge_legend = TRUE)
dev.off()

# 3. EPS (vector format for certain journals)
setEPS()
postscript(
    file       = "Figure_2_AUC_Heatmap_ComplexHeatmap.eps",
    width      = fig_width,
    height     = fig_height,
    paper      = "special",
    horizontal = FALSE,
    onefile    = FALSE,
    family     = font_family
)
draw(ht, 
     heatmap_legend_side = "right",
     padding = unit(c(5, 10, 5, 5), "mm"),
     merge_legend = TRUE)
dev.off()

# 4. PNG (for presentations/web, 300 DPI)
png(
    filename = "Figure_2_AUC_Heatmap_ComplexHeatmap.png",
    width    = fig_width,
    height   = fig_height,
    units    = "in",
    res      = 300
)
draw(ht, 
     heatmap_legend_side = "right",
     padding = unit(c(5, 10, 5, 5), "mm"),
     merge_legend = TRUE)
dev.off()

# --- Create a compact version with 5:1 aspect ratio ---------------------------
# Recalculate dimensions for 5:1 ratio
compact_height <- 3  # cm
compact_width <- compact_height * 5  # 5:1 ratio

ht_compact <- Heatmap(
    matrix              = mat,
    name                = legend_title,
    col                 = col_fun,
    cluster_rows        = FALSE,
    cluster_columns     = FALSE,
    show_row_names      = TRUE,
    show_column_names   = FALSE,
    row_names_side      = "left",
    row_names_gp        = gpar(fontsize = 10, fontfamily = font_family),
    rect_gp             = gpar(col = "white", lwd = 0.5),
    na_col              = "#E5E5E5",
    heatmap_legend_param = list(
        title           = legend_title,
        title_position  = "topleft",
        at              = legend_at,
        labels          = sprintf("%.2f", legend_at),
        title_gp        = gpar(fontsize = 10, fontface = "bold"),
        labels_gp       = gpar(fontsize = 9),
        legend_height   = unit(20, "mm"),  # Smaller for compact layout
        grid_width      = unit(3, "mm"),
        border          = "black"
    ),
    width               = unit(compact_width, "cm"),
    height              = unit(compact_height, "cm")
)

# Export compact version
cairo_pdf(
    filename  = "Figure_3_AUC_Heatmap_Compact_5to1.pdf",
    width     = 10,
    height    = 2,
    family    = font_family,
    pointsize = 10
)
draw(ht_compact, 
     heatmap_legend_side = "right",
     padding = unit(c(2, 5, 2, 2), "mm"),
     merge_legend = TRUE)
dev.off()

# --- Summary message ----------------------------------------------------------
cat("\n================== Heatmap Generation Complete ==================\n")
cat("Generated files (all without cell values):\n")
cat("\npheatmap versions:\n")
cat("- Figure_1_AUC_Heatmap_pheatmap.pdf\n")
cat("\nComplexHeatmap versions:\n")
cat("- Figure_2_AUC_Heatmap_ComplexHeatmap.pdf\n")
cat("- Figure_2_AUC_Heatmap_ComplexHeatmap.tiff (600 DPI)\n")
cat("- Figure_2_AUC_Heatmap_ComplexHeatmap.eps\n")
cat("- Figure_2_AUC_Heatmap_ComplexHeatmap.png (300 DPI)\n")
cat("- Figure_3_AUC_Heatmap_Compact_5to1.pdf (5:1 aspect ratio)\n")
cat("\nFeatures:\n")
cat("- Row labels (Sex-combined, Male, Female) on the left side\n")
cat("- Optimized legend size and positioning\n")
cat("- No cell value annotations for cleaner appearance\n")
cat("- Professional color schemes: Nature blue / Viridis\n")
cat("================================================================\n")

library(ComplexHeatmap)
library(ggplotify)
library(cowplot)
library(grid)
library(gtable)
library(ggplot2)
library(ggplotify)

p_heatmap <- as.ggplot(~ draw(ht_compact,
                              heatmap_legend_side = "right",
                              padding = unit(c(2,5,2,2), "mm"),
                              merge_legend = TRUE))