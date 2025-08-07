# =============================================================================
# Title:    Bar plot of protein cis/trans pQTL count distributions
# Author:   Hydraulik
# Date:     2025-07-30
# Modified: 2025-08-07
# =============================================================================

# 1. LOAD REQUIRED LIBRARIES
library(tidyverse)    # For data import, manipulation, and plotting
library(tidyr)        # For complete() function
library(scales)       # For axis formatting
library(RColorBrewer) # For color palettes

# 2. DEFINE FILE PATHS
base_dir <- "C:/Users/shaok/Desktop/ProteoNexus_Citation/GRABBING"
files <- list(
    all    = file.path(base_dir, "pip075_merged_sig_summ_all2_with_geneid1_filtered_annotated.csv"),
    female = file.path(base_dir, "pip075_merged_sig_summ_female2_with_geneid1_filtered_annotated.csv"),
    male   = file.path(base_dir, "pip075_merged_sig_summ_male2_with_geneid1_filtered_annotated.csv")
)

# 3. READ DATA
all_pqtl    <- read_csv(files$all)
female_pqtl <- read_csv(files$female)
male_pqtl   <- read_csv(files$male)

# 4. FUNCTION TO COMPUTE COUNT DISTRIBUTION OF CIS/TRANS PQTLS
compute_count_distribution <- function(df, label) {
    # Calculate counts per protein
    protein_counts <- df %>%
        group_by(protein_name) %>%
        summarise(
            cis_count   = sum(cis,   na.rm = TRUE),
            trans_count = sum(trans, na.rm = TRUE),
            .groups = "drop"
        )
    
    # Create count distribution for each type
    create_count_distribution <- function(counts, type_label) {
        data.frame(count = counts) %>%
            filter(count > 0) %>%
            count(count, name = "n_proteins") %>%
            mutate(
                type = type_label,
                dataset = label
            )
    }
    
    # Apply to both cis and trans counts
    cis_dist   <- create_count_distribution(protein_counts$cis_count,   "cis")
    trans_dist <- create_count_distribution(protein_counts$trans_count, "trans")
    
    # Combine results
    bind_rows(cis_dist, trans_dist)
}

# 5. COMPUTE COUNT DISTRIBUTIONS
dist_all    <- compute_count_distribution(all_pqtl,    "Sex-combined")
dist_female <- compute_count_distribution(female_pqtl, "Female")
dist_male   <- compute_count_distribution(male_pqtl,   "Male")

# 6. COMBINE ALL DISTRIBUTIONS
dist_combined <- bind_rows(dist_all, dist_female, dist_male) %>%
    # Set factor levels for desired order
    mutate(dataset = factor(dataset, levels = c("Sex-combined", "Male", "Female")))

# Determine the range of counts for consistent x-axis
max_count <- max(dist_combined$count)

# 7. CREATE PUBLICATION-READY BAR PLOT
# Define professional color palette
color_palette <- c(
    "Sex-combined" = "#F29B76",  # Coral
    "Male"         = "#5F8FCA",  # Blue
    "Female"       = "#9EB4DC"   # Light blue
)

# Create the plot with continuous x-axis
p <- ggplot(dist_combined,
            aes(x = count, y = n_proteins, fill = dataset)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.9), 
             width = 0.8, alpha = 0.9) +
    facet_wrap(~ type, scales = "free_y", ncol = 1,
               labeller = labeller(type = c(cis = "cis-pQTLs", 
                                            trans = "trans-pQTLs"))) +
    scale_fill_manual(values = color_palette, name = "Datasets") +
    scale_x_continuous(
        breaks = pretty_breaks(n = 10),
        expand = expansion(mult = c(0.02, 0.02))
    ) +
    scale_y_continuous(
        expand = expansion(mult = c(0, 0.1)),
        labels = comma
    ) +
    labs(
        x = "Number of pQTLs per protein",
        y = "Number of proteins"
    ) +
    theme_bw(base_size = 14) +
    theme(
        # Axis formatting
        axis.title = element_text(size = 18, face = "bold"),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.ticks = element_line(size = 0.5),
        
        # Panel and grid formatting
        panel.grid.major.x = element_line(color = "gray90", size = 0.5),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.border = element_rect(size = 1),
        
        # Facet formatting
        strip.background = element_rect(fill = "gray95", color = "gray20", size = 1),
        strip.text = element_text(size = 11, face = "bold"),
        
        # Legend formatting
        legend.position = "bottom",
        legend.title = element_text(size = 11, face = "bold"),
        legend.text = element_text(size = 10),
        legend.key.size = unit(0.8, "cm"),
        
        # Overall appearance
        plot.margin = margin(10, 10, 10, 10)
    )

# 8. SAVE PLOT IN HIGH RESOLUTION
ggsave(
    filename = "protein_pqtl_distribution_continuous.pdf",
    plot = p,
    width = 10,
    height = 10,
    dpi = 300,
    device = "pdf"
)

# Also save as PNG for quick viewing
ggsave(
    filename = "protein_pqtl_distribution_continuous.png",
    plot = p,
    width = 10,
    height = 10,
    dpi = 300,
    device = "png"
)

# 9. DISPLAY PLOT
print(p)

# 10. GENERATE SUMMARY STATISTICS TABLE
summary_stats <- dist_combined %>%
    group_by(dataset, type) %>%
    summarise(
        total_proteins = sum(n_proteins),
        mean_pqtl_count = weighted.mean(count, n_proteins),
        median_pqtl_count = median(rep(count, n_proteins)),
        max_pqtl_count = max(count),
        proteins_with_1_pqtl = sum(n_proteins[count == 1]),
        proteins_with_over_10_pqtls = sum(n_proteins[count > 10]),
        .groups = "drop"
    ) %>%
    arrange(dataset, type)

print("Summary Statistics:")
print(summary_stats)

# 11. CREATE ALTERNATIVE VISUALIZATION WITH LOG-TRANSFORMED X-AXIS (OPTIONAL)
# This can be useful if the distribution is highly skewed
p_log <- ggplot(dist_combined,
                aes(x = count, y = n_proteins, fill = dataset)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.9), 
             width = 0.8, alpha = 0.9) +
    facet_wrap(~ type, scales = "free_y", ncol = 1,
               labeller = labeller(type = c(cis = "cis-pQTLs", 
                                            trans = "trans-pQTLs"))) +
    scale_fill_manual(values = color_palette, name = "Datasets") +
    scale_x_log10(
        breaks = c(1, 2, 5, 10, 20, 50, 100, 200, 500),
        labels = c(1, 2, 5, 10, 20, 50, 100, 200, 500)
    ) +
    scale_y_continuous(
        expand = expansion(mult = c(0, 0.1)),
        labels = comma
    ) +
    labs(
        x = "Number of pQTLs per protein (log scale)",
        y = "Number of proteins"
    ) +
    theme_bw(base_size = 14) +
    theme(
        # Axis formatting
        axis.title = element_text(size = 18, face = "bold"),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.ticks = element_line(size = 0.5),
        
        # Panel and grid formatting
        panel.grid.major.x = element_line(color = "gray90", size = 0.5),
        panel.grid.minor.x = element_line(color = "gray95", size = 0.3),
        panel.grid.minor.y = element_blank(),
        panel.border = element_rect(size = 1),
        
        # Facet formatting
        strip.background = element_rect(fill = "gray95", color = "gray20", size = 1),
        strip.text = element_text(size = 11, face = "bold"),
        
        # Legend formatting
        legend.position = "bottom",
        legend.title = element_text(size = 11, face = "bold"),
        legend.text = element_text(size = 10),
        legend.key.size = unit(0.8, "cm"),
        
        # Overall appearance
        plot.margin = margin(10, 10, 10, 10)
    )

# Save log-scale version
ggsave(
    filename = "protein_pqtl_distribution_log_scale.pdf",
    plot = p_log,
    width = 10,
    height = 10,
    dpi = 300,
    device = "pdf"
)