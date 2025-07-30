# =============================================================================
# Title:    Bar plot of protein cis/trans pQTL count distributions by intervals
# Author:   Hydraulik
# Date:     2025-07-30
# =============================================================================

# 1. LOAD REQUIRED LIBRARIES
library(tidyverse)    # for data import, manipulation, and plotting
library(tidyr)        # for complete()
library(scales)       # for formatting
library(RColorBrewer) # for color palettes

# 2. DEFINE FILE PATHS
base_dir <- "C:/Users/shaok/Desktop/ProteoNexus_Citation/GRABBING"
files <- list(
    all    = file.path(base_dir, "merged_sig_summ_all2_with_geneid1_filtered_annotated.csv"),
    female = file.path(base_dir, "merged_sig_summ_female2_with_geneid1_filtered_annotated.csv"),
    male   = file.path(base_dir, "merged_sig_summ_male2_with_geneid1_filtered_annotated.csv")
)

# 3. READ DATA
all_pqtl    <- read_csv(files$all)
female_pqtl <- read_csv(files$female)
male_pqtl   <- read_csv(files$male)

# 4. FUNCTION TO COMPUTE BINNED DISTRIBUTION OF CIS/TRANS COUNTS
compute_binned_distribution <- function(df, label) {
    # Calculate counts per protein
    protein_counts <- df %>%
        group_by(protein_name) %>%
        summarise(
            cis_count   = sum(cis,   na.rm = TRUE),
            trans_count = sum(trans, na.rm = TRUE),
            .groups = "drop"
        )
    
    # Create bins for counts
    create_bins <- function(counts, type_label) {
        data.frame(count = counts) %>%
            filter(count > 0) %>%
            mutate(
                bin = case_when(
                    count <= 10  ~ "1-10",
                    count <= 20  ~ "11-20",
                    count <= 30  ~ "21-30",
                    count <= 40  ~ "31-40",
                    count <= 50  ~ "41-50",
                    count <= 60  ~ "51-60",
                    count <= 70  ~ "61-70",
                    count <= 80  ~ "71-80",
                    count <= 90  ~ "81-90",
                    count <= 100 ~ "91-100",
                    TRUE         ~ ">100"
                ),
                bin = factor(bin, levels = c("1-10", "11-20", "21-30", "31-40", "41-50",
                                             "51-60", "61-70", "71-80", "81-90", "91-100", ">100"))
            ) %>%
            count(bin, name = "n_proteins") %>%
            mutate(
                type = type_label,
                dataset = label
            )
    }
    
    # Apply binning to both cis and trans counts
    cis_bins   <- create_bins(protein_counts$cis_count,   "cis")
    trans_bins <- create_bins(protein_counts$trans_count, "trans")
    
    # Combine results
    bind_rows(cis_bins, trans_bins)
}

# 5. COMPUTE BINNED DISTRIBUTIONS
dist_all    <- compute_binned_distribution(all_pqtl,    "Sex-combined")
dist_female <- compute_binned_distribution(female_pqtl, "Female")
dist_male   <- compute_binned_distribution(male_pqtl,   "Male")

# 6. COMBINE ALL DISTRIBUTIONS
dist_combined <- bind_rows(dist_all, dist_female, dist_male) %>%
    # Set factor levels for desired order
    mutate(dataset = factor(dataset, levels = c("Sex-combined", "Male", "Female"))) %>%
    # Ensure all combinations exist (fill missing with 0)
    complete(
        dataset = unique(dataset),
        type    = unique(type),
        bin     = unique(bin),
        fill    = list(n_proteins = 0)
    )

# 7. CREATE PUBLICATION-READY BAR PLOT
# Define professional color palette
color_palette <- c(
    "Sex-combined" = "#F29B76",  # Purple
    "Male"         = "#5F8FCA",  # Teal
    "Female"       = "#9EB4DC"   # Yellow
)

# Create the plot
p <- ggplot(dist_combined,
            aes(x = bin, y = n_proteins, fill = dataset)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.9), 
             width = 0.8, alpha = 0.9) +
    facet_wrap(~ type, scales = "free_y", ncol = 1,
               labeller = labeller(type = c(cis = "cis-SNPs", 
                                            trans = "trans-SNPs"))) +
    scale_fill_manual(values = color_palette, name = "Dataset") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.1)),
                       labels = comma) +
    labs(
        x = "Number of pQTLs",
        y = "Number of proteins"
    ) +
    theme_bw(base_size = 12) +
    theme(
        # Axis formatting
        axis.title = element_text(size = 11, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 10),
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
        # legend.box.background = element_rect(color = "gray80", size = 0.5),
        
        # Overall appearance
        plot.margin = margin(10, 10, 10, 10)
    )

# 8. SAVE PLOT IN HIGH RESOLUTION
ggsave(
    filename = "protein_pqtl_distribution_barplot.pdf",
    plot = p,
    width = 8,
    height = 10,
    dpi = 300,
    device = "pdf"
)

# Also save as PNG for quick viewing
ggsave(
    filename = "protein_pqtl_distribution_barplot.png",
    plot = p,
    width = 8,
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
        median_bin = bin[which.max(n_proteins)],
        .groups = "drop"
    ) %>%
    pivot_wider(
        names_from = type,
        values_from = c(total_proteins, median_bin),
        names_glue = "{type}_{.value}"
    )

print("Summary Statistics:")
print(summary_stats)