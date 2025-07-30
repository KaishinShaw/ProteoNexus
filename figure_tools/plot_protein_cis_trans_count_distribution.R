# =============================================================================
# Title:    Plot protein cis/trans count distributions for all, female, and male
# Author:   Hydraulik
# Date:     2025-07-30
# =============================================================================

# 1. LOAD REQUIRED LIBRARIES
library(tidyverse)    # for data import, manipulation, and plotting
library(tidyr)        # for completing missing count bins

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

# 4. FUNCTION TO COMPUTE DISTRIBUTION OF CIS/TRANS COUNTS
compute_count_distribution <- function(df, dataset_label) {
    df %>%
        # Sum cis/trans flags per protein
        group_by(protein_name) %>%
        summarise(
            cis_count   = sum(cis,   na.rm = TRUE),
            trans_count = sum(trans, na.rm = TRUE),
            .groups = "drop"
        ) %>%
        # Pivot longer so we have one row per (protein, type)
        pivot_longer(
            cols      = c(cis_count, trans_count),
            names_to  = "type",
            values_to = "count"
        ) %>%
        # We only care about proteins with at least one event
        filter(count > 0) %>%
        # Count how many proteins have each count value
        count(type, count, name = "n_proteins") %>%
        # Tag with dataset label
        mutate(dataset = dataset_label)
}

# 5. COMPUTE DISTRIBUTIONS FOR EACH DATASET
dist_all    <- compute_count_distribution(all_pqtl,    "all")
dist_female <- compute_count_distribution(female_pqtl, "female")
dist_male   <- compute_count_distribution(male_pqtl,   "male")

# 6. COMBINE AND FILL MISSING BINS (so x-axis is uniform)
dist_combined <- bind_rows(dist_all, dist_female, dist_male) %>%
    # Ensure that for each (dataset, type) we have counts 1,2,... up to the max observed
    group_by(dataset, type) %>%
    complete(
        count = full_seq(count, 1),     # generates sequence 1,2,3,...,max
        fill  = list(n_proteins = 0)
    ) %>%
    ungroup()

# 7. PLOT SIX LINES USING ggplot2
p <- ggplot(dist_combined,
            aes(x = count,
                y = n_proteins,
                color = dataset,
                linetype = type)
) +
    geom_line(size = 1) +
    scale_color_brewer(palette = "Dark2", name = "Dataset") +
    scale_linetype_manual(
        values = c(cis_count   = "solid",
                   trans_count = "dashed"),
        labels = c(cis_count   = "cis",
                   trans_count = "trans"),
        name   = "Association Type"
    ) +
    scale_x_continuous(breaks = scales::pretty_breaks()) +
    labs(
        title    = "Distribution of Protein-level cis/trans pQTL Counts",
        subtitle = "Number of proteins with exactly 1, 2, 3, … cis or trans associations",
        x        = "Number of cis/trans associations per protein",
        y        = "Number of proteins",
        caption  = "Data: merged_sig_summ_*_with_geneid1_filtered_annotated.csv"
    ) +
    theme_minimal(base_size = 14) +
    theme(
        legend.position = "bottom",
        legend.box      = "vertical",
        plot.title      = element_text(face = "bold", size = 16),
        plot.subtitle   = element_text(size = 12)
    )

# 8. RENDER PLOT
print(p)