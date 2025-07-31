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
gpd_path <- "C:/Users/shaok/Desktop/ProteoNexus_Citation/GRABBING/GPD_enriched_annotated.tsv"

# Read GPD data
gpd <- read_tsv(gpd_path, col_types = cols())

# -------------------------------------------------------------------
# 3. Create SNP Type classification based on cis column
# -------------------------------------------------------------------
# Create SNP Type based on cis column values
gpd <- gpd %>%
    mutate(
        `SNP Type` = case_when(
            cis == 1 ~ "cis-SNPs",
            cis == 0 ~ "trans-SNPs",
            TRUE ~ NA_character_  # Handle any unexpected values
        )
    ) %>%
    filter(!is.na(`SNP Type`))  # Remove any rows with NA SNP Type

# -------------------------------------------------------------------
# 4. Deduplicate proteins within each SNP Type and gender
# -------------------------------------------------------------------
# Remove duplicate proteins within each gender and SNP Type combination
gpd_unique <- gpd %>%
    group_by(gender, `SNP Type`, protein) %>%
    slice(1) %>%  # Keep only the first occurrence of each protein
    ungroup()

message(paste("Original number of rows:", nrow(gpd)))
message(paste("Number of unique protein-SNP Type combinations:", nrow(gpd_unique)))

# -------------------------------------------------------------------
# 5. Pre-process: factor levels and color schemes
# -------------------------------------------------------------------
# Order of gender levels
dataset_levels <- c("all", "male", "female")

# SNP Type orders
snp_type_levels <- c("cis-SNPs", "trans-SNPs")

# Color-blind-friendly palette (Okabe-Ito palette) - using 2 colors
gpd_colors <- c("#cbbe7b", "#4a6d9f")  # Orange for cis-SNP, Sky blue for trans-SNP

# Apply ordered factors
gpd_unique <- gpd_unique %>%
    mutate(
        gender = factor(gender, levels = dataset_levels),
        `SNP Type` = factor(`SNP Type`, levels = snp_type_levels)
    )

# -------------------------------------------------------------------
# 6. Publication-quality percent-stacked bar function
# -------------------------------------------------------------------
create_publication_plot <- function(df, snp_colors) {
    ggplot(df, aes(x = gender, fill = `SNP Type`)) +
        geom_bar(
            position = "fill",
            width = 0.7,
            color = "black",
            linewidth = 0.4
        ) +
        # Y axis: fractions → 0-100, no "%" suffix
        scale_y_continuous(
            labels = label_number(scale = 100, accuracy = 1),
            breaks = seq(0, 1, 0.2),
            expand = c(0, 0),
            limits = c(0, 1)
        ) +
        scale_x_discrete(labels = c("Sex-combined", "Male", "Female")) +
        scale_fill_manual(values = snp_colors) +
        labs(
            x = "Dataset",
            y = "Proportion of protein (%)",
            fill = "SNP Type"  # Legend title
        ) +
        theme_bw(base_size = 12) +
        theme(
            # Axis formatting
            axis.title.x = element_text(face = "bold", size = 14, margin = margin(t = 12)),
            axis.title.y = element_text(face = "bold", size = 14, margin = margin(r = 12)),
            axis.text.x = element_text(color = "black", size = 12),
            axis.text.y = element_text(color = "black", size = 11),
            axis.line = element_line(color = "black", linewidth = 0.5),
            axis.ticks = element_line(color = "black", linewidth = 0.5),
            axis.ticks.length = unit(0.15, "cm"),
            
            # Legend formatting
            legend.position = "bottom",
            legend.direction = "horizontal",
            legend.title = element_text(face = "bold", size = 12, margin = margin(r = 10)),
            legend.text = element_text(size = 11),
            legend.key.size = unit(1, "lines"),
            legend.key = element_rect(color = "black", linewidth = 0.3),
            legend.margin = margin(t = 10),
            legend.box = "horizontal",
            legend.box.margin = margin(0, 0, 0, 0),
            legend.spacing.x = unit(0.5, "cm"),
            
            # Panel formatting
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_rect(fill = "white"),
            panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
            
            # Plot formatting
            plot.background = element_rect(fill = "white", color = NA),
            plot.margin = margin(t = 15, r = 15, b = 15, l = 15)
        ) +
        guides(
            fill = guide_legend(
                title.position = "left",     # Position title on the left
                title.hjust = 0,            # Align title to the left
                nrow = 1,                   # Single row for 2 categories
                byrow = TRUE,
                keywidth = unit(1.5, "lines"),
                keyheight = unit(0.9, "lines")
            )
        )
}

# -------------------------------------------------------------------
# 7. Generate plot
# -------------------------------------------------------------------
plot_gpd <- create_publication_plot(gpd_unique, gpd_colors)

# Display on screen
print(plot_gpd)

# -------------------------------------------------------------------
# 8. Save high-resolution figures
# -------------------------------------------------------------------
dir.create("figures", showWarnings = FALSE)

save_pub_plot <- function(plot_obj, basename, width = 8, height = 5.5) {
    # PNG format for digital viewing
    ggsave(
        filename = paste0("figures/", basename, ".png"),
        plot = plot_obj,
        width = width,
        height = height,
        dpi = 600,
        bg = "white"
    )
    
    # PDF format for vector graphics (preferred for publication)
    ggsave(
        filename = paste0("figures/", basename, ".pdf"),
        plot = plot_obj,
        width = width,
        height = height,
        device = cairo_pdf,
        bg = "white"
    )
    
    # TIFF format for some journals
    ggsave(
        filename = paste0("figures/", basename, ".tiff"),
        plot = plot_obj,
        width = width,
        height = height,
        dpi = 600,
        compression = "lzw",
        bg = "white"
    )
}

save_pub_plot(plot_gpd, "gpd_protein_regulation_stacked_bar_publication")

# -------------------------------------------------------------------
# 9. Export protein proportions table (based on unique proteins)
# -------------------------------------------------------------------
calculate_protein_proportions <- function(df, dataset_name) {
    df %>%
        group_by(gender, `SNP Type`) %>%
        tally(name = "n_proteins") %>%
        mutate(
            proportion = n_proteins / sum(n_proteins) * 100,
            dataset = dataset_name,
            gender = factor(gender,
                            levels = c("all", "male", "female"),
                            labels = c("Sex-combined", "Male", "Female"))
        ) %>%
        select(dataset, gender, `SNP Type`, n_proteins, proportion) %>%
        arrange(gender, `SNP Type`)
}

gpd_protein_props <- calculate_protein_proportions(gpd_unique, "GPD")

write_csv(gpd_protein_props, "figures/protein_regulation_proportions_table.csv")

# Print summary statistics
message("GPD protein regulation proportions (unique proteins):")
print(gpd_protein_props)

# -------------------------------------------------------------------
# 10. Additional analysis: Overall protein regulation counts
# -------------------------------------------------------------------
overall_protein_counts <- gpd_unique %>%
    group_by(`SNP Type`) %>%
    summarise(
        unique_proteins = n(),
        percentage = n() / nrow(gpd_unique) * 100
    )

message("\nOverall protein regulation distribution:")
print(overall_protein_counts)

# Gender-specific protein counts
gender_protein_counts <- gpd_unique %>%
    group_by(gender, `SNP Type`) %>%
    summarise(
        protein_count = n(),
        .groups = "drop"
    ) %>%
    group_by(gender) %>%
    mutate(
        percentage = protein_count / sum(protein_count) * 100
    )

message("\nGender-specific protein regulation distribution:")
print(gender_protein_counts)

# Protein regulation diversity analysis
protein_regulation_stats <- gpd %>%
    group_by(gender, `SNP Type`) %>%
    summarise(
        unique_proteins = n_distinct(protein),
        total_snp_protein_pairs = n(),
        avg_snps_per_protein = total_snp_protein_pairs / unique_proteins,
        .groups = "drop"
    )

message("\nProtein regulation statistics by gender and SNP Type:")
print(protein_regulation_stats)

# -------------------------------------------------------------------
# 11. Print plot information for publication
# -------------------------------------------------------------------
message("\n=== Figure Information for Publication ===")
message("Plot dimensions: 8 × 5.5 inches")
message("Resolution: 600 DPI")
message("Formats saved: PNG, PDF, TIFF")
message("Color scheme: Okabe-Ito color-blind friendly palette")
message("SNP Types: cis-SNP (orange), trans-SNP (sky blue)")
message(paste("Total number of SNP-protein associations:", nrow(gpd)))
message(paste("Number of unique proteins analyzed:", nrow(gpd_unique)))
message(paste("Number of unique cis-regulated proteins:", sum(gpd_unique$`SNP Type` == "cis-SNP")))
message(paste("Number of unique trans-regulated proteins:", sum(gpd_unique$`SNP Type` == "trans-SNP")))