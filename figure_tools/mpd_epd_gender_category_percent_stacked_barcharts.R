# -------------------------------------------------------------------
# 1. Load required libraries
# -------------------------------------------------------------------
# install.packages(c("readr","dplyr","ggplot2","scales"))
library(readr)       # Fast TSV import
library(dplyr)       # Data manipulation
library(ggplot2)     # Publication-quality plotting
library(scales)      # Percent formatting

# -------------------------------------------------------------------
# 2. Define file paths and read data
# -------------------------------------------------------------------
mpd_path <- "C:/Users/shaok/Desktop/ProteoNexus_Citation/GRABBING/M-P-D_annotated.tsv"
epd_path <- "C:/Users/shaok/Desktop/ProteoNexus_Citation/GRABBING/E-P-D_annotated.tsv"

mpd <- read_tsv(mpd_path, col_types = cols()) 
epd <- read_tsv(epd_path, col_types = cols())

# -------------------------------------------------------------------
# 3. Pre-process: Set custom factor levels and define color schemes
# -------------------------------------------------------------------
# Dataset levels (new order: Sex-combined, Male, Female)
dataset_levels <- c("all", "male", "female")

# MPD categories in specified order
mpd_categories <- c("Anthropometric", 
                    "Blood and Unrine Test", 
                    "Psychosocial Factors", 
                    "Other")

# EPD categories in specified order
epd_categories <- c("Residential air pollution", 
                    "Greenspace and coastal proximity", 
                    "Sociodemographics", 
                    "Healthy lifestyle")

# Professional color schemes for publication
# Using colorblind-friendly palette with good contrast
mpd_colors <- c("#4A4A4A", "#EB8928", "#275892", "#B5B5B5")  # Okabe-Ito palette
epd_colors <- c("#4A4A4A", "#EB8928", "#275892", "#B5B5B5")  # Complementary professional colors

# Apply factor levels with new order
mpd <- mpd %>%
    mutate(
        gender   = factor(gender, levels = dataset_levels),
        Category = factor(Category, levels = mpd_categories)
    )

epd <- epd %>%
    mutate(
        gender   = factor(gender, levels = dataset_levels),
        Category = factor(Category, levels = epd_categories)
    )

# -------------------------------------------------------------------
# 4. Enhanced function for publication-quality percent-stacked bar charts
# -------------------------------------------------------------------
create_publication_plot <- function(df, category_colors, plot_title = NULL) {
    p <- ggplot(df, aes(x = gender, fill = Category)) +
        geom_bar(position = "fill", width = 0.7, color = "black", linewidth = 0.4) +
        scale_y_continuous(
            labels = percent_format(accuracy = 1),
            expand = c(0, 0),
            limits = c(0, 1),
            breaks = seq(0, 1, 0.2)
        ) +
        scale_fill_manual(values = category_colors) +
        scale_x_discrete(labels = c("Sex-combined", "Male", "Female")) +
        labs(
            x = "Dataset",
            y = "Proportion (%)",
            fill = NULL  # Remove legend title for cleaner appearance
        ) +
        theme_bw(base_size = 12) +
        theme(
            # Axis styling
            axis.title.x = element_text(face = "bold", size = 14, margin = margin(t = 12)),
            axis.title.y = element_text(face = "bold", size = 14, margin = margin(r = 12)),
            axis.text.x = element_text(color = "black", size = 12, face = "plain"),
            axis.text.y = element_text(color = "black", size = 11),
            axis.line = element_line(color = "black", linewidth = 0.5),
            axis.ticks = element_line(color = "black", linewidth = 0.5),
            axis.ticks.length = unit(0.15, "cm"),
            
            # Legend styling - positioned at bottom
            legend.position = "bottom",
            legend.direction = "horizontal",
            legend.text = element_text(size = 11),
            legend.key.size = unit(0.9, "lines"),
            legend.key = element_rect(color = "black", linewidth = 0.3),
            legend.margin = margin(t = 15),
            legend.box.margin = margin(t = 5),
            legend.background = element_rect(fill = "white", color = NA),
            
            # Panel styling
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_rect(fill = "white"),
            panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
            plot.background = element_rect(fill = "white", color = NA),
            
            # Plot margins
            plot.margin = margin(t = 15, r = 15, b = 15, l = 15)
        ) +
        guides(fill = guide_legend(
            nrow = 2,  # Arrange legend items in 2 rows
            byrow = TRUE,
            keywidth = unit(1.5, "lines"),
            keyheight = unit(0.9, "lines"),
            label.position = "right"
        ))
    
    # Add title if provided
    if (!is.null(plot_title)) {
        p <- p + ggtitle(plot_title) +
            theme(plot.title = element_text(face = "bold", size = 16, hjust = 0.5, margin = margin(b = 15)))
    }
    
    return(p)
}

# -------------------------------------------------------------------
# 5. Generate publication-quality plots
# -------------------------------------------------------------------
plot_mpd <- create_publication_plot(mpd, mpd_colors)
plot_epd <- create_publication_plot(epd, epd_colors)

# Display plots
print(plot_mpd)
print(plot_epd)

# -------------------------------------------------------------------
# 6. Save high-resolution figures for journal submission
# -------------------------------------------------------------------
# Create output directory if it doesn't exist
dir.create("figures", showWarnings = FALSE)

# Save MPD plot
ggsave(
    filename = "figures/mpd_percent_stacked_bar_publication.png",
    plot = plot_mpd,
    width = 8,      # Optimal width for journal columns
    height = 5.5,   # Adjusted for legend at bottom
    dpi = 600,      # High resolution for publication
    bg = "white"    # White background
)

ggsave(
    filename = "figures/mpd_percent_stacked_bar_publication.pdf",
    plot = plot_mpd,
    width = 8,
    height = 5.5,
    device = cairo_pdf,  # Vector format for publication
    bg = "white"
)

ggsave(
    filename = "figures/mpd_percent_stacked_bar_publication.tiff",
    plot = plot_mpd,
    width = 8,
    height = 5.5,
    dpi = 600,
    compression = "lzw",  # TIFF compression for smaller file size
    bg = "white"
)

# Save EPD plot
ggsave(
    filename = "figures/epd_percent_stacked_bar_publication.png",
    plot = plot_epd,
    width = 8,
    height = 5.5,
    dpi = 600,
    bg = "white"
)

ggsave(
    filename = "figures/epd_percent_stacked_bar_publication.pdf",
    plot = plot_epd,
    width = 8,
    height = 5.5,
    device = cairo_pdf,
    bg = "white"
)

ggsave(
    filename = "figures/epd_percent_stacked_bar_publication.tiff",
    plot = plot_epd,
    width = 8,
    height = 5.5,
    dpi = 600,
    compression = "lzw",
    bg = "white"
)

# -------------------------------------------------------------------
# 7. Optional: Create combined figure for main manuscript or supplementary materials
# -------------------------------------------------------------------
library(patchwork)  # For combining plots

# Add panel labels for combined figure
plot_mpd_labeled <- plot_mpd + 
    labs(subtitle = "A. Medical/Physiological/Demographic (MPD) Variables") +
    theme(plot.subtitle = element_text(face = "bold", size = 13, hjust = 0))

plot_epd_labeled <- plot_epd + 
    labs(subtitle = "B. Environmental/Psychological/Demographic (EPD) Variables") +
    theme(plot.subtitle = element_text(face = "bold", size = 13, hjust = 0))

# Combine plots vertically
combined_plot <- plot_mpd_labeled / plot_epd_labeled +
    plot_layout(heights = c(1, 1))

# Save combined figure
ggsave(
    filename = "figures/combined_stacked_bars_publication.png",
    plot = combined_plot,
    width = 8,
    height = 11,
    dpi = 600,
    bg = "white"
)

ggsave(
    filename = "figures/combined_stacked_bars_publication.pdf",
    plot = combined_plot,
    width = 8,
    height = 11,
    device = cairo_pdf,
    bg = "white"
)

ggsave(
    filename = "figures/combined_stacked_bars_publication.tiff",
    plot = combined_plot,
    width = 8,
    height = 11,
    dpi = 600,
    compression = "lzw",
    bg = "white"
)

# -------------------------------------------------------------------
# 8. Generate summary statistics table for publication
# -------------------------------------------------------------------
# Function to calculate category proportions
calculate_proportions <- function(df, dataset_name) {
    df %>%
        group_by(gender, Category) %>%
        summarise(n = n(), .groups = "drop_last") %>%
        mutate(proportion = n / sum(n) * 100) %>%
        ungroup() %>%
        mutate(
            dataset = dataset_name,
            gender = factor(gender, 
                            levels = c("all", "male", "female"),
                            labels = c("Sex-combined", "Male", "Female"))
        ) %>%
        select(dataset, gender, Category, n, proportion) %>%
        arrange(gender, Category)
}

# Calculate proportions for both datasets
mpd_props <- calculate_proportions(mpd, "MPD")
epd_props <- calculate_proportions(epd, "EPD")

# Combine and save as CSV for supplementary table
all_props <- bind_rows(mpd_props, epd_props)
write_csv(all_props, "figures/category_proportions_table.csv")

# Print summary for verification
print("MPD Category Proportions:")
print(mpd_props)
print("\nEPD Category Proportions:")
print(epd_props)