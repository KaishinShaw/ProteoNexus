# figure_AUC_SCI_publication.R
################################################################################
# Publication-Quality AUC Heatmap Generation for Scientific Journals
# Author: Hydraulik
# Date: 2025-01-30
# Description: Generate high-quality heatmaps for AUC values with professional
#              scientific color schemes suitable for SCI journal submission
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

# 2. Diverging palette for positive/negative values (Publication standard)
sci_diverging <- colorRampPalette(c("#053061", "#2166AC", "#4393C3", 
                                    "#92C5DE", "#D1E5F0", "#F7F7F7", 
                                    "#FDDBC7", "#F4A582", "#D6604D", 
                                    "#B2182B", "#67001F"))(100)

# 3. Viridis palette (colorblind-friendly)
viridis_palette <- viridis(100, option = "D", begin = 0.05, end = 0.95)

# 4. Custom scientific palette (Cell/Science style)
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
# Calculate dimensions for 5:1 ratio
n_cols <- ncol(mat)
n_rows <- nrow(mat)

# For pheatmap, we need to calculate cellwidth and cellheight
# to achieve 5:1 ratio for the heatmap part
cellheight_base <- 15  # Base height for each cell
cellwidth_calc <- (cellheight_base * n_rows * 5) / n_cols  # Calculate width for 5:1 ratio

# Define output dimensions (inches) - adjusted for 5:1 ratio
fig_width <- 10   # Total figure width including legend
fig_height <- 2.5  # Compact height for 5:1 ratio

# Create high-resolution PDF (5:1 ratio version)
cairo_pdf(
    filename  = "Figure_1_AUC_Heatmap_pheatmap_5to1.pdf",
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
    cluster_rows      = FALSE,
    cluster_cols      = FALSE,
    show_rownames     = TRUE,
    show_colnames     = FALSE,        # Hide column names
    fontsize_row      = 10,           # Adjusted for compact layout
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
    display_numbers   = FALSE          # No numbers for this version
)

dev.off()

# Alternative version with values displayed (normal proportions)
cairo_pdf(
    filename  = "Figure_1_AUC_Heatmap_pheatmap_with_values.pdf",
    width     = 8,
    height    = 3,
    family    = font_family,
    pointsize = 10
)

pheatmap(
    mat               = mat,
    color             = nature_blue,
    breaks            = breaks_seq,
    cluster_rows      = FALSE,
    cluster_cols      = FALSE,
    show_rownames     = TRUE,
    show_colnames     = FALSE,        # Hide column names
    fontsize_row      = 12,
    fontsize_number   = 8,            # Font size for cell values
    cellwidth         = 25,
    cellheight        = 22,
    legend            = TRUE,
    legend_breaks     = legend_breaks,
    legend_labels     = sprintf("%.2f", legend_breaks),
    border_color      = "white",
    na_col            = "#E5E5E5",
    main              = "",
    display_numbers   = TRUE,         # Show values
    number_format     = "%.3f",       # Format for displayed numbers
    number_color      = "black"       # Color for numbers
)

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
    legend_title <- "AUC Value"
    legend_at <- seq(-max_abs, max_abs, length.out = 5)
} else {
    # Use sequential palette
    col_fun <- colorRamp2(
        breaks = seq(data_range[1], data_range[2], length.out = 100),
        colors = viridis_palette
    )
    legend_title <- "AUC Score"
    legend_at <- pretty(data_range, n = 5)
}

# --- Create ComplexHeatmap with publication settings --------------------------
# Calculate dimensions for 5:1 ratio
heatmap_height <- 1.5  # Total height in cm for 3 rows
heatmap_width <- heatmap_height * 5  # Width = 5 * height for 5:1 ratio

# Version without cell values (5:1 ratio)
ht_5to1 <- Heatmap(
    matrix              = mat,
    name                = legend_title,
    col                 = col_fun,
    
    # Clustering settings
    cluster_rows        = FALSE,
    cluster_columns     = FALSE,
    show_row_dend       = FALSE,
    show_column_dend    = FALSE,
    
    # Row/column names
    show_row_names      = TRUE,
    show_column_names   = FALSE,      # Hide column names
    row_names_side      = "left",
    
    # Text formatting
    row_names_gp        = gpar(fontsize = 10, fontfamily = font_family),
    
    # Cell formatting
    rect_gp             = gpar(col = "white", lwd = 0.5),
    na_col              = "#E5E5E5",
    
    # Legend settings
    heatmap_legend_param = list(
        title           = legend_title,
        title_position  = "topleft",
        at              = legend_at,
        labels          = sprintf("%.2f", legend_at),
        title_gp        = gpar(fontsize = 10, fontface = "bold", 
                               fontfamily = font_family),
        labels_gp       = gpar(fontsize = 9, fontfamily = font_family),
        legend_height   = unit(25, "mm"),  # Smaller legend for compact layout
        grid_width      = unit(3, "mm"),
        border          = "black"
    ),
    
    # Size settings for 5:1 ratio
    width               = unit(heatmap_width, "cm"),
    height              = unit(heatmap_height, "cm")
)

# Version with cell values (normal proportions)
ht_with_values <- Heatmap(
    matrix              = mat,
    name                = legend_title,
    col                 = col_fun,
    
    # Clustering settings
    cluster_rows        = FALSE,
    cluster_columns     = FALSE,
    show_row_dend       = FALSE,
    show_column_dend    = FALSE,
    
    # Row/column names
    show_row_names      = TRUE,
    show_column_names   = FALSE,      # Hide column names
    row_names_side      = "left",
    
    # Text formatting
    row_names_gp        = gpar(fontsize = 12, fontfamily = font_family),
    
    # Cell formatting
    rect_gp             = gpar(col = "white", lwd = 0.5),
    na_col              = "#E5E5E5",
    
    # Cell text - show values
    cell_fun            = function(j, i, x, y, width, height, fill) {
        if (!is.na(mat[i, j])) {
            # Determine text color based on background
            val <- mat[i, j]
            bg_lightness <- col2rgb(fill)[1,1] / 255
            text_col <- ifelse(bg_lightness > 0.5, "black", "white")
            
            grid.text(sprintf("%.3f", val), x, y, 
                      gp = gpar(fontsize = 8, col = text_col, fontfamily = font_family))
        }
    },
    
    # Legend settings
    heatmap_legend_param = list(
        title           = legend_title,
        title_position  = "topleft",
        at              = legend_at,
        labels          = sprintf("%.2f", legend_at),
        title_gp        = gpar(fontsize = 10, fontface = "bold", 
                               fontfamily = font_family),
        labels_gp       = gpar(fontsize = 9, fontfamily = font_family),
        legend_height   = unit(30, "mm"),
        grid_width      = unit(4, "mm"),
        border          = "black"
    ),
    
    # Size settings (normal proportions)
    width               = unit(ncol(mat) * 0.8, "cm"),
    height              = unit(nrow(mat) * 0.8, "cm")
)

# --- Export to multiple formats -----------------------------------------------
# Calculate figure dimensions for 5:1 ratio version
fig_width_5to1  <- 10  # Total width including legend
fig_height_5to1 <- 2   # Compact height

# 1. PDF without values (5:1 ratio)
cairo_pdf(
    filename  = "Figure_2_AUC_Heatmap_ComplexHeatmap_5to1.pdf",
    width     = fig_width_5to1,
    height    = fig_height_5to1,
    family    = font_family,
    pointsize = 10
)
draw(ht_5to1, 
     heatmap_legend_side = "right",
     padding = unit(c(2, 5, 2, 2), "mm"),
     merge_legend = TRUE)
dev.off()

# 2. PDF with values (normal proportions)
fig_width_normal  <- max(7, ncol(mat) * 0.3 + 2.5)
fig_height_normal <- max(3, nrow(mat) * 0.4 + 1.5)

cairo_pdf(
    filename  = "Figure_2_AUC_Heatmap_ComplexHeatmap_with_values.pdf",
    width     = fig_width_normal,
    height    = fig_height_normal,
    family    = font_family,
    pointsize = 10
)
draw(ht_with_values, 
     heatmap_legend_side = "right",
     padding = unit(c(5, 10, 5, 5), "mm"),
     merge_legend = TRUE)
dev.off()

# 3. TIFF (600 DPI for print quality - 5:1 ratio)
tiff(
    filename    = "Figure_2_AUC_Heatmap_ComplexHeatmap_5to1.tiff",
    width       = fig_width_5to1,
    height      = fig_height_5to1,
    units       = "in",
    res         = 600,  # High resolution for print
    compression = "lzw"
)
draw(ht_5to1, 
     heatmap_legend_side = "right",
     padding = unit(c(2, 5, 2, 2), "mm"),
     merge_legend = TRUE)
dev.off()

# 4. EPS (vector format - 5:1 ratio)
setEPS()
postscript(
    file       = "Figure_2_AUC_Heatmap_ComplexHeatmap_5to1.eps",
    width      = fig_width_5to1,
    height     = fig_height_5to1,
    paper      = "special",
    horizontal = FALSE,
    onefile    = FALSE,
    family     = font_family
)
draw(ht_5to1, 
     heatmap_legend_side = "right",
     padding = unit(c(2, 5, 2, 2), "mm"),
     merge_legend = TRUE)
dev.off()

# 5. PNG (for presentations/web, 300 DPI - 5:1 ratio)
png(
    filename = "Figure_2_AUC_Heatmap_ComplexHeatmap_5to1.png",
    width    = fig_width_5to1,
    height   = fig_height_5to1,
    units    = "in",
    res      = 300
)
draw(ht_5to1, 
     heatmap_legend_side = "right",
     padding = unit(c(2, 5, 2, 2), "mm"),
     merge_legend = TRUE)
dev.off()

# --- Summary message ----------------------------------------------------------
cat("\n================== Heatmap Generation Complete ==================\n")
cat("Generated files:\n")
cat("\n5:1 Ratio versions (without values):\n")
cat("- Figure_1_AUC_Heatmap_pheatmap_5to1.pdf\n")
cat("- Figure_2_AUC_Heatmap_ComplexHeatmap_5to1.pdf\n")
cat("- Figure_2_AUC_Heatmap_ComplexHeatmap_5to1.tiff (600 DPI)\n")
cat("- Figure_2_AUC_Heatmap_ComplexHeatmap_5to1.eps\n")
cat("- Figure_2_AUC_Heatmap_ComplexHeatmap_5to1.png (300 DPI)\n")
cat("\nNormal proportion versions (with values):\n")
cat("- Figure_1_AUC_Heatmap_pheatmap_with_values.pdf\n")
cat("- Figure_2_AUC_Heatmap_ComplexHeatmap_with_values.pdf\n")
cat("\nColor schemes used:\n")
cat("- Sequential: Nature-style blue gradient / Viridis\n")
cat("- Diverging: Scientific red-white-blue (if applicable)\n")
cat("================================================================\n")
