# =========================================================================
# Stage 1: Setup - Loading Required Libraries
# =========================================================================
# Core packages for data manipulation and visualization
library(readr)
library(dplyr)
library(ggplot2)
library(scales)
library(tidyr)  # Added explicitly for pivot_longer

# Optional: For better font handling (uncomment if needed)
# install.packages("showtext")
# library(showtext)
# showtext_auto()

# =========================================================================
# Stage 2: Data Loading
# =========================================================================
# Define file paths using forward slashes for cross-platform compatibility
path_all    <- "C:/Users/shaok/Desktop/ProteoNexus_Citation/GRABBING/merged_sig_summ_all2_with_geneid1_filtered_annotated.csv"
path_female <- "C:/Users/shaok/Desktop/ProteoNexus_Citation/GRABBING/merged_sig_summ_female2_with_geneid1_filtered_annotated.csv"
path_male   <- "C:/Users/shaok/Desktop/ProteoNexus_Citation/GRABBING/merged_sig_summ_male2_with_geneid1_filtered_annotated.csv"

# Load data with error handling
tryCatch({
    all_pqtl    <- read_csv(path_all, show_col_types = FALSE)
    female_pqtl <- read_csv(path_female, show_col_types = FALSE)
    male_pqtl   <- read_csv(path_male, show_col_types = FALSE)
}, error = function(e) {
    stop("Error loading files. Please check if the file paths are correct and the files exist.\nOriginal error: ", e$message)
})

# =========================================================================
# Stage 3: Data Processing and Transformation
# =========================================================================
# Calculate cis and trans pQTL counts for each dataset
summary_all <- all_pqtl %>% 
    summarise(cis = sum(cis, na.rm = TRUE), trans = sum(trans, na.rm = TRUE))

summary_female <- female_pqtl %>% 
    summarise(cis = sum(cis, na.rm = TRUE), trans = sum(trans, na.rm = TRUE))

summary_male <- male_pqtl %>% 
    summarise(cis = sum(cis, na.rm = TRUE), trans = sum(trans, na.rm = TRUE))

# Combine and reshape data for plotting
plot_data <- bind_rows(
    "Sex-combined" = summary_all,
    "Male"         = summary_male,
    "Female"       = summary_female,
    .id = "Dataset"
) %>%
    pivot_longer(
        cols = c(cis, trans),
        names_to = "pQTL_Type",
        values_to = "Count"
    )

# Convert to ordered factors
plot_data$Dataset <- factor(plot_data$Dataset, 
                            levels = c("Sex-combined", "Male", "Female"))
plot_data$pQTL_Type <- factor(plot_data$pQTL_Type, 
                              levels = c("cis", "trans"), 
                              labels = c("cis-SNPs", "trans-SNPs"))

# =========================================================================
# Stage 4: Create Publication-Quality Visualization
# =========================================================================
# Define publication-standard color palette
custom_colors <- c("cis-SNPs" = "#cbbe7b", "trans-SNPs" = "#4a6d9f")

# Create the bar chart
pQTL_bar_chart <- ggplot(plot_data, aes(x = Dataset, y = Count, fill = pQTL_Type)) +
    
    # Bar layer with optimized spacing
    geom_bar(stat = "identity", 
             position = position_dodge(width = 0.8), 
             width = 0.7,
             color = "black",       # Add black outline
             linewidth = 0.3) +     # Thin outline for definition
    
    # Add value labels above bars
    geom_text(aes(label = format(Count, big.mark = ",")), 
              position = position_dodge(width = 0.8), 
              vjust = -0.5, 
              size = 3.5,
              fontface = "plain") + # Remove Arial font specification
    
    # Apply custom colors
    scale_fill_manual(values = custom_colors,
                      name = "pQTL Type") +
    
    # Axis labels
    labs(
        x = "Dataset",
        y = "Number of Significant pQTLs",
        fill = "pQTL Type"
    ) +
    
    # Y-axis formatting
    scale_y_continuous(
        expand = expansion(mult = c(0, 0.1)),  # No expansion at bottom, 10% at top
        labels = function(x) format(x, big.mark = ",", scientific = FALSE),
        breaks = pretty_breaks(n = 5)
    ) +
    
    # Apply minimalist theme
    theme_classic(base_size = 12) +
    
    # Detailed theme customization for publication standards
    theme(
        # Text elements - using default fonts to avoid errors
        text = element_text(color = "black"),
        
        # Axis text
        axis.text.x = element_text(size = 11, color = "black", margin = margin(t = 5)),
        axis.text.y = element_text(size = 11, color = "black", margin = margin(r = 5)),
        
        # Axis titles
        axis.title.x = element_text(size = 12, face = "bold", margin = margin(t = 15)),
        axis.title.y = element_text(size = 12, face = "bold", margin = margin(r = 15)),
        
        # Axis appearance
        axis.line = element_line(colour = "black", linewidth = 0.5),
        axis.ticks = element_line(colour = "black", linewidth = 0.5),
        axis.ticks.length = unit(0.15, "cm"),
        
        # Legend at bottom with improved spacing
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.title = element_text(size = 11, face = "bold", margin = margin(b = 5)),
        legend.text = element_text(size = 10),
        legend.key.size = unit(0.6, "cm"),
        legend.key.width = unit(1.2, "cm"),
        legend.spacing.x = unit(0.3, "cm"),
        legend.margin = margin(t = 15),
        legend.box.margin = margin(0, 0, 0, 0),
        
        # Clean background
        panel.background = element_rect(fill = "white", colour = NA),
        panel.grid = element_blank(),
        
        # Plot margins
        plot.margin = margin(t = 10, r = 10, b = 10, l = 10, unit = "pt")
    )

# Display the plot
print(pQTL_bar_chart)

# =========================================================================
# Stage 5: Export for Scientific Publication
# =========================================================================
# Create output directory if it doesn't exist
output_dir <- "publication_figures"
if (!dir.exists(output_dir)) {
    dir.create(output_dir)
}

# Function to save plots with consistent settings
save_publication_plot <- function(plot, filename, width, height, dpi = 300) {
    ggsave(
        filename = file.path(output_dir, filename),
        plot = plot,
        width = width,
        height = height,
        units = "in",
        dpi = dpi
    )
    cat("Saved:", file.path(output_dir, filename), "\n")
}

# Save in multiple formats for different journal requirements

# 1. PDF format (vector graphics, preferred for most journals)
save_publication_plot(pQTL_bar_chart, 
                      "pQTL_counts_by_dataset.pdf", 
                      width = 7, 
                      height = 5)

# 2. High-resolution PNG (300 DPI)
save_publication_plot(pQTL_bar_chart, 
                      "pQTL_counts_by_dataset_300dpi.png", 
                      width = 7, 
                      height = 5, 
                      dpi = 300)

# 3. Extra high-resolution PNG (600 DPI for some journals)
save_publication_plot(pQTL_bar_chart, 
                      "pQTL_counts_by_dataset_600dpi.png", 
                      width = 7, 
                      height = 5, 
                      dpi = 600)

# 4. TIFF format with LZW compression
ggsave(
    filename = file.path(output_dir, "pQTL_counts_by_dataset.tiff"),
    plot = pQTL_bar_chart,
    width = 7,
    height = 5,
    units = "in",
    dpi = 300,
    compression = "lzw"
)

# 5. Single-column width (Nature/Science format: 89mm ≈ 3.5 inches)
save_publication_plot(pQTL_bar_chart + 
                          theme(legend.direction = "vertical",
                                legend.position = "right",
                                legend.box.spacing = unit(0, "pt")), 
                      "pQTL_counts_single_column.pdf", 
                      width = 3.5, 
                      height = 4)

# 6. Double-column width (standard: 183mm ≈ 7.2 inches)
save_publication_plot(pQTL_bar_chart, 
                      "pQTL_counts_double_column.pdf", 
                      width = 7.2, 
                      height = 5)

# 7. EPS format (some journals still require this)
ggsave(
    filename = file.path(output_dir, "pQTL_counts_by_dataset.eps"),
    plot = pQTL_bar_chart,
    width = 7,
    height = 5,
    units = "in",
    device = "eps"
)

# =========================================================================
# Stage 6: Generate Summary Statistics Table
# =========================================================================
# Create a summary table for manuscript
summary_table <- plot_data %>%
    pivot_wider(names_from = pQTL_Type, values_from = Count) %>%
    mutate(
        Total = `cis-SNPs` + `trans-SNPs`,
        `cis-SNPs (%)` = round(`cis-SNPs` / Total * 100, 1),
        `trans-SNPs (%)` = round(`trans-SNPs` / Total * 100, 1)
    ) %>%
    select(Dataset, `cis-SNPs`, `trans-SNPs`, Total, 
           `cis-SNPs (%)`, `trans-SNPs (%)`)

# Print summary table
cat("\n=== pQTL Summary Statistics ===\n")
print(summary_table, n = Inf)

# Save summary table as CSV
write_csv(summary_table, file.path(output_dir, "pQTL_summary_statistics.csv"))

cat("\n=== All files have been saved to:", output_dir, "===\n")