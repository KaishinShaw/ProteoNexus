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
    
    # Create count distribution for each type with grouping >20
    create_count_distribution <- function(counts, type_label) {
        # Create initial distribution
        count_df <- data.frame(count = counts) %>%
            filter(count > 0)
        
        # Group counts into categories: 1-20 as individual values, >20 as one group
        count_summary <- count_df %>%
            mutate(
                count_category = ifelse(count > 20, ">20", as.character(count))
            ) %>%
            count(count_category, name = "n_proteins") %>%
            mutate(
                type = type_label,
                dataset = label
            )
        
        return(count_summary)
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

# 6. COMBINE ALL DISTRIBUTIONS AND PREPARE FOR PLOTTING
# Create ordered factor levels for x-axis
x_levels <- c(as.character(1:20), ">20")

dist_combined <- bind_rows(dist_all, dist_female, dist_male) %>%
    # Set factor levels for desired order
    mutate(
        dataset = factor(dataset, levels = c("Sex-combined", "Male", "Female")),
        count_category = factor(count_category, levels = x_levels)
    ) %>%
    # Ensure all combinations exist (fill missing with 0)
    complete(
        dataset = unique(dataset),
        type = unique(type),
        count_category = factor(x_levels, levels = x_levels),
        fill = list(n_proteins = 0)
    )

# 7. CREATE PUBLICATION-READY BAR PLOT
# Define professional color palette
color_palette <- c(
    "Sex-combined" = "#F29B76",  # Coral
    "Male"         = "#5F8FCA",  # Blue
    "Female"       = "#9EB4DC"   # Light blue
)

# Create the plot
p <- ggplot(dist_combined,
            aes(x = count_category, y = n_proteins, fill = dataset)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.9), 
             width = 0.8, alpha = 0.9) +
    facet_wrap(~ type, scales = "free_y", ncol = 1,
               labeller = labeller(type = c(cis = "cis-pQTLs", 
                                            trans = "trans-pQTLs"))) +
    scale_fill_manual(values = color_palette, name = "Datasets") +
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
        axis.text.x = element_text(size = 10, angle = 0, hjust = 0.5),
        axis.text.y = element_text(size = 10),
        axis.ticks = element_line(size = 0.5),
        
        # Panel and grid formatting
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
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
    filename = "protein_pqtl_distribution_grouped.pdf",
    plot = p,
    width = 10,
    height = 10,
    dpi = 300,
    device = "pdf"
)

# Also save as PNG for quick viewing
ggsave(
    filename = "protein_pqtl_distribution_grouped.png",
    plot = p,
    width = 10,
    height = 10,
    dpi = 300,
    device = "png"
)

# 9. DISPLAY PLOT
print(p)

# 10. GENERATE DETAILED SUMMARY STATISTICS TABLE
# First, get the original count data for accurate statistics
compute_detailed_stats <- function(df, label) {
    protein_counts <- df %>%
        group_by(protein_name) %>%
        summarise(
            cis_count   = sum(cis,   na.rm = TRUE),
            trans_count = sum(trans, na.rm = TRUE),
            .groups = "drop"
        )
    
    stats <- bind_rows(
        data.frame(
            dataset = label,
            type = "cis",
            total_proteins = sum(protein_counts$cis_count > 0),
            mean_pqtl_count = mean(protein_counts$cis_count[protein_counts$cis_count > 0]),
            median_pqtl_count = median(protein_counts$cis_count[protein_counts$cis_count > 0]),
            max_pqtl_count = max(protein_counts$cis_count),
            proteins_with_1_pqtl = sum(protein_counts$cis_count == 1),
            proteins_with_2_5_pqtls = sum(protein_counts$cis_count >= 2 & protein_counts$cis_count <= 5),
            proteins_with_6_10_pqtls = sum(protein_counts$cis_count >= 6 & protein_counts$cis_count <= 10),
            proteins_with_11_20_pqtls = sum(protein_counts$cis_count >= 11 & protein_counts$cis_count <= 20),
            proteins_with_over_20_pqtls = sum(protein_counts$cis_count > 20)
        ),
        data.frame(
            dataset = label,
            type = "trans",
            total_proteins = sum(protein_counts$trans_count > 0),
            mean_pqtl_count = mean(protein_counts$trans_count[protein_counts$trans_count > 0]),
            median_pqtl_count = median(protein_counts$trans_count[protein_counts$trans_count > 0]),
            max_pqtl_count = max(protein_counts$trans_count),
            proteins_with_1_pqtl = sum(protein_counts$trans_count == 1),
            proteins_with_2_5_pqtls = sum(protein_counts$trans_count >= 2 & protein_counts$trans_count <= 5),
            proteins_with_6_10_pqtls = sum(protein_counts$trans_count >= 6 & protein_counts$trans_count <= 10),
            proteins_with_11_20_pqtls = sum(protein_counts$trans_count >= 11 & protein_counts$trans_count <= 20),
            proteins_with_over_20_pqtls = sum(protein_counts$trans_count > 20)
        )
    )
    return(stats)
}

# Compute statistics for all datasets
stats_all <- compute_detailed_stats(all_pqtl, "Sex-combined")
stats_female <- compute_detailed_stats(female_pqtl, "Female")
stats_male <- compute_detailed_stats(male_pqtl, "Male")

summary_stats <- bind_rows(stats_all, stats_female, stats_male) %>%
    mutate(dataset = factor(dataset, levels = c("Sex-combined", "Male", "Female"))) %>%
    arrange(dataset, type)

print("Detailed Summary Statistics:")
print(summary_stats)

# 11. EXPORT SUMMARY STATISTICS TO CSV
write_csv(summary_stats, "protein_pqtl_summary_statistics.csv")

# 12. CREATE A SUPPLEMENTARY PLOT SHOWING THE TAIL DISTRIBUTION (OPTIONAL)
# This focuses on proteins with high pQTL counts
tail_data <- bind_rows(
    all_pqtl %>% 
        group_by(protein_name) %>%
        summarise(
            cis_count = sum(cis, na.rm = TRUE),
            trans_count = sum(trans, na.rm = TRUE),
            .groups = "drop"
        ) %>%
        mutate(dataset = "Sex-combined"),
    female_pqtl %>% 
        group_by(protein_name) %>%
        summarise(
            cis_count = sum(cis, na.rm = TRUE),
            trans_count = sum(trans, na.rm = TRUE),
            .groups = "drop"
        ) %>%
        mutate(dataset = "Female"),
    male_pqtl %>% 
        group_by(protein_name) %>%
        summarise(
            cis_count = sum(cis, na.rm = TRUE),
            trans_count = sum(trans, na.rm = TRUE),
            .groups = "drop"
        ) %>%
        mutate(dataset = "Male")
) %>%
    pivot_longer(cols = c(cis_count, trans_count), 
                 names_to = "type", 
                 values_to = "count") %>%
    filter(count > 20) %>%
    mutate(
        type = str_replace(type, "_count", ""),
        dataset = factor(dataset, levels = c("Sex-combined", "Male", "Female"))
    )

if (nrow(tail_data) > 0) {
    p_tail <- ggplot(tail_data, aes(x = count, fill = dataset)) +
        geom_histogram(binwidth = 5, position = "dodge", alpha = 0.9) +
        facet_wrap(~ type, scales = "free", ncol = 1,
                   labeller = labeller(type = c(cis = "cis-pQTLs (>20)", 
                                                trans = "trans-pQTLs (>20)"))) +
        scale_fill_manual(values = color_palette, name = "Datasets") +
        scale_x_continuous(breaks = seq(25, max(tail_data$count) + 5, by = 10)) +
        labs(
            x = "Number of pQTLs per protein",
            y = "Number of proteins",
            title = "Distribution of proteins with >20 pQTLs"
        ) +
        theme_bw(base_size = 14) +
        theme(
            plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
            axis.title = element_text(size = 14, face = "bold"),
            strip.text = element_text(size = 11, face = "bold"),
            legend.position = "bottom"
        )
    
    ggsave(
        filename = "protein_pqtl_distribution_tail.pdf",
        plot = p_tail,
        width = 8,
        height = 8,
        dpi = 300,
        device = "pdf"
    )
}