# figure_AUC_SCI_publication_ggplot.R
################################################################################
# Publication-Quality AUC Heatmap Generation for Scientific Journals
# Author: Hydraulik (Modified for ggplot2)
# Date: 2025-01-30
# Description: Generate high-quality heatmaps for AUC values without cell 
#              annotations using ggplot2, suitable for SCI journal submission
################################################################################

# --- Package Management -------------------------------------------------------
# Function to install and load packages
load_packages <- function(packages) {
    for (pkg in packages) {
        if (!requireNamespace(pkg, quietly = TRUE)) {
            install.packages(pkg, repos = "https://cran.r-project.org")
        }
        library(pkg, character.only = TRUE)
    }
}

# Load required packages
required_packages <- c("ggplot2", "tidyr", "dplyr", "scales", "viridis", 
                       "RColorBrewer", "Cairo", "cowplot", "ggthemes")

load_packages(required_packages)

# --- Font Setup ---------------------------------------------------------------
# Define font family with fallback
font_family <- "sans"  # Use default sans-serif font for most formats
font_family_eps <- "Helvetica"  # Use Helvetica for EPS (postscript compatible)

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

# Convert to matrix and then to long format for ggplot2
mat <- as.matrix(data_ordered)

# Check for all NA matrix
if (all(is.na(mat))) {
    stop("Error: All values in the matrix are NA")
}

# Convert to long format for ggplot2
data_long <- data.frame(
    Group = rep(rownames(mat), ncol(mat)),
    Model = rep(colnames(mat), each = nrow(mat)),
    AUC = as.vector(mat)
)

# Factor levels to maintain order
data_long$Group <- factor(data_long$Group, levels = rev(rownames(mat)))
data_long$Model <- factor(data_long$Model, levels = colnames(mat))

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
# SECTION 1: Simple Publication-Quality Heatmap with ggplot2
################################################################################

# Calculate data range for consistent scaling
data_range <- range(mat, na.rm = TRUE)

# Create base ggplot heatmap
p1 <- ggplot(data_long, aes(x = Model, y = Group, fill = AUC)) +
    geom_tile(color = "white", size = 0.5) +  # White borders between cells
    scale_fill_gradientn(
        colors = nature_blue,
        limits = data_range,
        breaks = pretty(data_range, n = 5),
        labels = function(x) sprintf("%.2f", x),
        na.value = "#E5E5E5",  # Light gray for NA values
        name = "AUC"
    ) +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    theme_minimal(base_size = 12, base_family = font_family) +
    theme(
        # Panel settings
        panel.grid = element_blank(),
        panel.background = element_rect(fill = "white", color = NA),
        
        # Axis settings
        axis.text.x = element_blank(),  # Hide x-axis labels for cleaner look
        axis.text.y = element_text(size = 12, color = "black", hjust = 1),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        
        # Legend settings
        legend.position = "right",
        legend.title = element_text(size = 11, face = "bold"),
        legend.text = element_text(size = 10),
        legend.key.height = unit(0.8, "cm"),
        legend.key.width = unit(0.4, "cm"),
        
        # Plot margins
        plot.margin = margin(10, 10, 10, 10)
    ) +
    coord_fixed(ratio = ncol(mat) / nrow(mat) * 0.6)  # Adjust aspect ratio

# Calculate dimensions for optimal layout
n_cols <- ncol(mat)
n_rows <- nrow(mat)
fig_width <- max(8, n_cols * 0.5 + 2.5)
fig_height <- max(3, n_rows * 0.8 + 1.5)

# Save the first version
ggsave(
    filename = "Figure_1_AUC_Heatmap_ggplot_simple.pdf",
    plot = p1,
    device = cairo_pdf,
    width = fig_width,
    height = fig_height,
    units = "in",
    dpi = 300
)

################################################################################
# SECTION 2: Advanced Publication-Quality Heatmap with ggplot2
################################################################################

# Determine appropriate color scheme
has_negative <- any(mat < 0, na.rm = TRUE)
has_positive <- any(mat > 0, na.rm = TRUE)

if (has_negative && has_positive) {
    # Use diverging palette centered at zero
    max_abs <- max(abs(data_range))
    color_palette <- c("#053061", "#4393C3", "#F7F7F7", "#D6604D", "#67001F")
    breaks_seq <- seq(-max_abs, max_abs, length.out = 100)
    legend_breaks <- seq(-max_abs, max_abs, length.out = 5)
} else {
    # Use sequential palette
    color_palette <- viridis_palette
    breaks_seq <- seq(data_range[1], data_range[2], length.out = 100)
    legend_breaks <- pretty(data_range, n = 5)
}

# Create advanced heatmap with viridis palette
p2 <- ggplot(data_long, aes(x = Model, y = Group, fill = AUC)) +
    geom_tile(color = "white", size = 0.5) +
    scale_fill_gradientn(
        colors = color_palette,
        limits = range(breaks_seq),
        breaks = legend_breaks,
        labels = function(x) sprintf("%.2f", x),
        na.value = "#E5E5E5",
        name = "AUC"
    ) +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    theme_minimal(base_size = 12, base_family = font_family) +
    theme(
        # Panel settings
        panel.grid = element_blank(),
        panel.background = element_rect(fill = "white", color = NA),
        
        # Axis settings
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 12, color = "black", hjust = 1),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        
        # Legend settings
        legend.position = "right",
        legend.title = element_text(size = 11, face = "bold"),
        legend.text = element_text(size = 10),
        legend.key.height = unit(min(40, n_rows * 10), "mm"),
        legend.key.width = unit(4, "mm"),
        legend.background = element_rect(color = "black", size = 0.5),
        
        # Plot margins
        plot.margin = margin(5, 10, 5, 5, "mm")
    ) +
    coord_fixed(ratio = ncol(mat) / nrow(mat) * 0.6)

# Export to multiple formats

# 1. PDF (vector format, best for publications)
ggsave(
    filename = "Figure_2_AUC_Heatmap_ggplot_advanced.pdf",
    plot = p2,
    device = cairo_pdf,
    width = fig_width,
    height = fig_height,
    units = "in",
    dpi = 300
)

# 2. TIFF (600 DPI for print quality)
ggsave(
    filename = "Figure_2_AUC_Heatmap_ggplot_advanced.tiff",
    plot = p2,
    device = "tiff",
    width = fig_width,
    height = fig_height,
    units = "in",
    dpi = 600,
    compression = "lzw"
)

# 3. EPS (vector format for certain journals) - with special font handling
# Create version with EPS-compatible font
p2_eps <- p2 + theme(text = element_text(family = font_family_eps))

# Try cairo_ps first (better font support)
tryCatch({
    ggsave(
        filename = "Figure_2_AUC_Heatmap_ggplot_advanced.eps",
        plot = p2_eps,
        device = cairo_ps,
        width = fig_width,
        height = fig_height,
        units = "in"
    )
}, error = function(e) {
    # If cairo_ps fails, use regular postscript with Helvetica font
    message("cairo_ps failed, using standard postscript device")
    ggsave(
        filename = "Figure_2_AUC_Heatmap_ggplot_advanced.eps",
        plot = p2_eps,
        device = "eps",
        width = fig_width,
        height = fig_height,
        units = "in"
    )
})

# 4. PNG (for presentations/web, 300 DPI)
ggsave(
    filename = "Figure_2_AUC_Heatmap_ggplot_advanced.png",
    plot = p2,
    device = "png",
    width = fig_width,
    height = fig_height,
    units = "in",
    dpi = 300
)

################################################################################
# SECTION 3: Compact Version with 5:1 Aspect Ratio
################################################################################

# Create compact version with adjusted theme
p3 <- ggplot(data_long, aes(x = Model, y = Group, fill = AUC)) +
    geom_tile(color = "white", size = 0.5) +
    scale_fill_gradientn(
        colors = color_palette,
        limits = range(breaks_seq),
        breaks = legend_breaks,
        labels = function(x) sprintf("%.2f", x),
        na.value = "#E5E5E5",
        name = "AUC"
    ) +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    theme_minimal(base_size = 10, base_family = font_family) +
    theme(
        # Panel settings
        panel.grid = element_blank(),
        panel.background = element_rect(fill = "white", color = NA),
        # Axis settings
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = rel(1), color = "black", hjust = 1),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        # Legend settings
        legend.position = "right",
        legend.title = element_text(size = rel(0.9), face = "bold"),
        legend.text = element_text(size = rel(0.8)),
        # legend.background = element_rect(color = "black", size = 0.5),
        # Plot margins
        plot.margin = margin(0.02, 0.05, 0.02, 0.02, "npc")
    ) +
    coord_fixed(ratio = 5)

# Export compact version
ggsave(
    filename = "Figure_3_AUC_Heatmap_Compact_5to1.pdf",
    plot = p3,
    device = cairo_pdf,
    width = 10,
    height = 2,
    units = "in",
    dpi = 300
)

################################################################################
# SECTION 4: Alternative Style with Cell Borders
################################################################################

# Create version with pronounced cell borders
p4 <- ggplot(data_long, aes(x = Model, y = Group, fill = AUC)) +
    geom_tile(color = "black", size = 0.8) +  # Black borders for emphasis
    scale_fill_gradientn(
        colors = cell_palette,  # Use Cell/Science color palette
        limits = data_range,
        breaks = pretty(data_range, n = 5),
        labels = function(x) sprintf("%.2f", x),
        na.value = "#E5E5E5",
        name = "AUC"
    ) +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    theme_classic(base_size = 12, base_family = font_family) +
    theme(
        # Axis settings
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 12, color = "black"),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        
        # Legend settings
        legend.position = "right",
        legend.title = element_text(size = 11, face = "bold"),
        legend.text = element_text(size = 10),
        
        # Plot margins
        plot.margin = margin(10, 10, 10, 10)
    ) +
    coord_fixed(ratio = ncol(mat) / nrow(mat) * 0.6)

# Save alternative style
ggsave(
    filename = "Figure_4_AUC_Heatmap_ggplot_cell_style.pdf",
    plot = p4,
    device = cairo_pdf,
    width = fig_width,
    height = fig_height,
    units = "in",
    dpi = 300
)

################################################################################
# SECTION 5: Additional Export Options Without EPS
################################################################################

# If EPS is problematic, create SVG as an alternative vector format
ggsave(
    filename = "Figure_2_AUC_Heatmap_ggplot_advanced.svg",
    plot = p2,
    device = "svg",
    width = fig_width,
    height = fig_height,
    units = "in"
)

# --- Summary message ----------------------------------------------------------
cat("\n================== Heatmap Generation Complete ==================\n")
cat("Generated files (all without cell values, using ggplot2):\n")
cat("\nggplot2 versions:\n")
cat("- Figure_1_AUC_Heatmap_ggplot_simple.pdf (Nature blue palette)\n")
cat("- Figure_2_AUC_Heatmap_ggplot_advanced.pdf (Viridis palette)\n")
cat("- Figure_2_AUC_Heatmap_ggplot_advanced.tiff (600 DPI)\n")
if (file.exists("Figure_2_AUC_Heatmap_ggplot_advanced.eps")) {
    cat("- Figure_2_AUC_Heatmap_ggplot_advanced.eps\n")
} else {
    cat("- Figure_2_AUC_Heatmap_ggplot_advanced.svg (alternative to EPS)\n")
}
cat("- Figure_2_AUC_Heatmap_ggplot_advanced.png (300 DPI)\n")
cat("- Figure_3_AUC_Heatmap_Compact_5to1.pdf (5:1 aspect ratio)\n")
cat("- Figure_4_AUC_Heatmap_ggplot_cell_style.pdf (Cell/Science style)\n")
cat("\nFeatures:\n")
cat("- Row labels (Sex-combined, Male, Female) on the left side\n")
cat("- Optimized legend size and positioning\n")
cat("- No cell value annotations for cleaner appearance\n")
cat("- Professional color schemes: Nature blue / Viridis / Cell palette\n")
cat("- All plots created using ggplot2 for consistency\n")
cat("================================================================\n")