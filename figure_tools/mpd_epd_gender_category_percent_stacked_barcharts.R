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
# Gender levels
gender_levels <- c("male", "female", "all")

# MPD categories in specified order
mpd_categories <- c("Anthropometric", 
                    "Blood and Unrine Test", 
                    "Psychosocial Factors", 
                    "Other")
mpd_colors <- c("#ffd7a5", "#e1c0a3", "#ffa07a", "#a9d7c0")

# EPD categories in specified order
epd_categories <- c("Residential air pollution", 
                    "Greenspace and coastal proximity", 
                    "Sociodemographics", 
                    "Healthy lifestyle")
epd_colors <- c("#ffd7a5", "#e1c0a3", "#ffa07a", "#a9d7c0")

# Apply factor levels
mpd <- mpd %>%
    mutate(
        gender   = factor(gender, levels = gender_levels),
        Category = factor(Category, levels = mpd_categories)
    )

epd <- epd %>%
    mutate(
        gender   = factor(gender, levels = gender_levels),
        Category = factor(Category, levels = epd_categories)
    )

# -------------------------------------------------------------------
# 4. Enhanced function for publication-quality percent-stacked bar charts
# -------------------------------------------------------------------
create_publication_plot <- function(df, category_colors) {
    ggplot(df, aes(x = gender, fill = Category)) +
        geom_bar(position = "fill", width = 0.6, color = "black", size = 0.3) +
        scale_y_continuous(
            labels = percent_format(accuracy = 1),
            expand = c(0, 0),
            limits = c(0, 1)
        ) +
        scale_fill_manual(values = category_colors) +
        scale_x_discrete(labels = c("Male", "Female", "All")) +
        labs(
            x = "Gender",
            y = "Proportion (%)",
            fill = NULL  # Remove legend title for cleaner appearance
        ) +
        theme_classic(base_size = 12) +
        theme(
            # Axis styling
            axis.title.x = element_text(face = "bold", size = 13, margin = margin(t = 10)),
            axis.title.y = element_text(face = "bold", size = 13, margin = margin(r = 10)),
            axis.text = element_text(color = "black", size = 11),
            axis.line = element_line(color = "black", size = 0.5),
            axis.ticks = element_line(color = "black", size = 0.5),
            
            # Legend styling - positioned at bottom
            legend.position = "bottom",
            legend.direction = "horizontal",
            legend.text = element_text(size = 10),
            legend.key.size = unit(0.8, "lines"),
            legend.margin = margin(t = 10),
            legend.box.margin = margin(t = 5),
            
            # Panel styling
            panel.grid = element_blank(),
            panel.background = element_rect(fill = "white"),
            plot.background = element_rect(fill = "white"),
            
            # Remove extra space
            plot.margin = margin(t = 10, r = 10, b = 10, l = 10)
        ) +
        guides(fill = guide_legend(
            nrow = 2,  # Arrange legend items in 2 rows
            byrow = TRUE,
            keywidth = unit(1.2, "lines"),
            keyheight = unit(0.8, "lines")
        ))
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
# Save MPD plot
ggsave(
    filename = "mpd_percent_stacked_bar_publication.png",
    plot = plot_mpd,
    width = 7,     # Slightly wider to accommodate legend
    height = 5,    # Adjusted for legend at bottom
    dpi = 600,     # Higher DPI for publication
    bg = "white"   # White background
)

ggsave(
    filename = "mpd_percent_stacked_bar_publication.pdf",
    plot = plot_mpd,
    width = 7,
    height = 5,
    device = cairo_pdf,  # Vector format for publication
    bg = "white"
)

# Save EPD plot
ggsave(
    filename = "epd_percent_stacked_bar_publication.png",
    plot = plot_epd,
    width = 7,
    height = 5,
    dpi = 600,
    bg = "white"
)

ggsave(
    filename = "epd_percent_stacked_bar_publication.pdf",
    plot = plot_epd,
    width = 7,
    height = 5,
    device = cairo_pdf,
    bg = "white"
)

# -------------------------------------------------------------------
# 7. Optional: Create combined figure for supplementary materials
# -------------------------------------------------------------------
library(patchwork)  # For combining plots

combined_plot <- plot_mpd / plot_epd +
    plot_layout(heights = c(1, 1)) +
    plot_annotation(tag_levels = 'A') &
    theme(plot.tag = element_text(size = 14, face = "bold"))

ggsave(
    filename = "combined_stacked_bars_publication.png",
    plot = combined_plot,
    width = 7,
    height = 10,
    dpi = 600,
    bg = "white"
)

ggsave(
    filename = "combined_stacked_bars_publication.pdf",
    plot = combined_plot,
    width = 7,
    height = 10,
    device = cairo_pdf,
    bg = "white"
)