# -------------------------------------------------------------------
# 1. Load required libraries
# -------------------------------------------------------------------
# install.packages(c("readr","dplyr","ggplot2","scales","patchwork"))
library(readr)       # Fast TSV import
library(dplyr)       # Data manipulation
library(ggplot2)     # Publication-quality plotting
library(scales)      # Number formatting
library(patchwork)   # Combining plots

# -------------------------------------------------------------------
# 2. Define file paths and read data
# -------------------------------------------------------------------
mpd_path <- "C:/Users/shaok/Desktop/ProteoNexus_Citation/GRABBING/M-P-D_annotated.tsv"
epd_path <- "C:/Users/shaok/Desktop/ProteoNexus_Citation/GRABBING/E-P-D_annotated.tsv"

mpd <- read_tsv(mpd_path, col_types = cols())
epd <- read_tsv(epd_path, col_types = cols())

# -------------------------------------------------------------------
# 3. Pre-process: factor levels and color schemes
# -------------------------------------------------------------------
# Order of gender levels
dataset_levels <- c("all", "male", "female")

# Category orders
mpd_categories <- c("Anthropometric",
                    "Blood and Unrine Test",
                    "Psychosocial Factors",
                    "Other")
epd_categories <- c("Residential air pollution",
                    "Greenspace and coastal proximity",
                    "Sociodemographics",
                    "Healthy lifestyle")

# Color-blind‐friendly palette (Okabe-Ito palette)
mpd_colors <- c("#B5B5B5", "#4A4A4A", "#EB8928", "#275892")
epd_colors <- mpd_colors

# Apply ordered factors
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
# 4. Publication-quality percent-stacked bar function
# -------------------------------------------------------------------
create_publication_plot <- function(df, category_colors) {
    ggplot(df, aes(x = gender, fill = Category)) +
        geom_bar(
            position    = "fill",
            width       = 0.7,
            color       = "black",
            linewidth   = 0.4
        ) +
        # Y axis: fractions → 0-100, no "%" suffix
        scale_y_continuous(
            labels = label_number(scale = 100, accuracy = 1),
            breaks = seq(0, 1, 0.2),
            expand = c(0, 0),
            limits = c(0, 1)
        ) +
        scale_x_discrete(labels = c("Sex-combined", "Male", "Female")) +
        scale_fill_manual(values = category_colors) +
        labs(
            x    = "Dataset",
            y    = "Proportion of pathway (%)",
            fill = "Category"  # Added legend title
        ) +
        theme_bw(base_size = 12) +
        theme(
            # Axis formatting
            axis.title.x      = element_text(face = "bold", size = 14, margin = margin(t = 12)),
            axis.title.y      = element_text(face = "bold", size = 14, margin = margin(r = 12)),
            axis.text.x       = element_text(color = "black", size = 12),
            axis.text.y       = element_text(color = "black", size = 11),
            axis.line         = element_line(color = "black", linewidth = 0.5),
            axis.ticks        = element_line(color = "black", linewidth = 0.5),
            axis.ticks.length = unit(0.15, "cm"),
            
            # Legend formatting
            legend.position   = "bottom",
            legend.direction  = "horizontal",
            legend.title      = element_text(face = "bold", size = 12, margin = margin(r = 10)),
            legend.text       = element_text(size = 11),
            legend.key.size   = unit(1, "lines"),
            legend.key        = element_rect(color = "black", linewidth = 0.3),
            legend.margin     = margin(t = 10),
            legend.box        = "horizontal",
            legend.box.margin = margin(0, 0, 0, 0),
            legend.spacing.x  = unit(0.5, "cm"),
            
            # Panel formatting
            panel.grid.major  = element_blank(),
            panel.grid.minor  = element_blank(),
            panel.background  = element_rect(fill = "white"),
            panel.border      = element_rect(color = "black", fill = NA, linewidth = 0.8),
            
            # Plot formatting
            plot.background   = element_rect(fill = "white", color = NA),
            plot.margin       = margin(t = 15, r = 15, b = 15, l = 15)
        ) +
        guides(
            fill = guide_legend(
                title.position = "left",     # Position title on the left
                title.hjust    = 0,          # Align title to the left
                nrow           = 2,
                byrow          = TRUE,
                keywidth       = unit(1.5, "lines"),
                keyheight      = unit(0.9, "lines")
            )
        )
}

# -------------------------------------------------------------------
# 5. Generate individual plots
# -------------------------------------------------------------------
plot_mpd <- create_publication_plot(mpd, mpd_colors)
plot_epd <- create_publication_plot(epd, epd_colors)

# Display on screen
print(plot_mpd)
print(plot_epd)

# -------------------------------------------------------------------
# 6. Save high-resolution figures
# -------------------------------------------------------------------
dir.create("figures", showWarnings = FALSE)

save_pub_plot <- function(plot_obj, basename, width = 8, height = 5.5) {
    # PNG format for digital viewing
    ggsave(
        filename = paste0("figures/", basename, ".png"),
        plot     = plot_obj,
        width    = width,
        height   = height,
        dpi      = 600,
        bg       = "white"
    )
    
    # PDF format for vector graphics (preferred for publication)
    ggsave(
        filename = paste0("figures/", basename, ".pdf"),
        plot     = plot_obj,
        width    = width,
        height   = height,
        device   = cairo_pdf,
        bg       = "white"
    )
    
    # TIFF format for some journals
    ggsave(
        filename     = paste0("figures/", basename, ".tiff"),
        plot         = plot_obj,
        width        = width,
        height       = height,
        dpi          = 600,
        compression  = "lzw",
        bg           = "white"
    )
}

save_pub_plot(plot_mpd, "mpd_percent_stacked_bar_publication")
save_pub_plot(plot_epd, "epd_percent_stacked_bar_publication")

# -------------------------------------------------------------------
# 7. Combine plots with proper legend alignment
# -------------------------------------------------------------------
# Method 1: Using patchwork with collected guides
tryCatch({
    combined_plot <- (plot_mpd + theme(legend.position = "none")) / 
        (plot_epd + theme(legend.position = "none")) +
        plot_layout(heights = c(1, 1), guides = "collect") &
        theme(legend.position = "bottom")
    
    save_pub_plot(combined_plot, "combined_stacked_bars_publication", width = 8, height = 11)
}, error = function(e) {
    message("Combined plot with patchwork failed. Trying alternative method...")
    
    # Method 2: Simple combination
    tryCatch({
        combined_plot <- plot_mpd / plot_epd
        save_pub_plot(combined_plot, "combined_stacked_bars_publication", width = 8, height = 11)
    }, error = function(e2) {
        message("Alternative method also failed. Saving plots separately.")
        
        # Method 3: Using cowplot as fallback
        if(requireNamespace("cowplot", quietly = TRUE)) {
            library(cowplot)
            combined_plot_cowplot <- plot_grid(
                plot_mpd + theme(legend.position = "none"), 
                plot_epd + theme(legend.position = "none"),
                ncol = 1,
                align = "v",
                rel_heights = c(1, 1),
                labels = c("A", "B"),
                label_size = 16,
                label_fontface = "bold"
            )
            
            # Extract legend from one plot
            legend <- get_legend(plot_mpd)
            
            # Combine plots with legend
            final_plot <- plot_grid(
                combined_plot_cowplot,
                legend,
                ncol = 1,
                rel_heights = c(10, 1)
            )
            
            save_pub_plot(final_plot, "combined_stacked_bars_publication_cowplot", 
                          width = 8, height = 11)
        }
    })
})

# -------------------------------------------------------------------
# 8. Export category proportions table
# -------------------------------------------------------------------
calculate_proportions <- function(df, dataset_name) {
    df %>%
        group_by(gender, Category) %>%
        tally(name = "n") %>%
        mutate(
            proportion = n / sum(n) * 100,
            dataset    = dataset_name,
            gender     = factor(gender,
                                levels = c("all","male","female"),
                                labels = c("Sex-combined","Male","Female"))
        ) %>%
        select(dataset, gender, Category, n, proportion) %>%
        arrange(gender, Category)
}

mpd_props <- calculate_proportions(mpd, "MPD")
epd_props <- calculate_proportions(epd, "EPD")
all_props <- bind_rows(mpd_props, epd_props)

write_csv(all_props, "figures/category_proportions_table.csv")

# Print summary statistics
message("MPD proportions:")
print(mpd_props)
message("\nEPD proportions:")
print(epd_props)

# -------------------------------------------------------------------
# 9. Print plot dimensions for publication
# -------------------------------------------------------------------
message("\n=== Figure Information for Publication ===")
message("Individual plot dimensions: 8 × 5.5 inches")
message("Combined plot dimensions: 8 × 11 inches")
message("Resolution: 600 DPI")
message("Formats saved: PNG, PDF, TIFF")
message("Color scheme: Okabe-Ito color-blind friendly palette")