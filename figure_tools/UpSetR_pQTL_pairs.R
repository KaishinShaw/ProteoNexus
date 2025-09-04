# 1. Load Required Libraries ------------------------------------------------
library(readr)      # Fast CSV import
library(dplyr)      # Data manipulation
library(tidyr)      # Data reshaping
library(ggplot2)    # Plotting
library(scales)     # Axis formatting

# 2. Define File Paths & Import Data ----------------------------------------
base_dir    <- "C:/Users/shaok/Desktop/ProteoNexus_Citation/GRABBING"
file_all    <- file.path(base_dir, "merged_sig_summ_all2_with_geneid1_filtered_annotated.csv")
file_female <- file.path(base_dir, "merged_sig_summ_female2_with_geneid1_filtered_annotated.csv")
file_male   <- file.path(base_dir, "merged_sig_summ_male2_with_geneid1_filtered_annotated.csv")

# Read each dataset with error handling
tryCatch({
    df_all    <- read_csv(file_all,    show_col_types = FALSE)
    df_female <- read_csv(file_female, show_col_types = FALSE)
    df_male   <- read_csv(file_male,   show_col_types = FALSE)
}, error = function(e) {
    stop("Failed to load data. Please verify file paths.\nDetails: ", e$message)
})

# 3. Install and Load Additional Required Libraries -------------------------
if (!require("UpSetR")) {
    install.packages("UpSetR")
    library(UpSetR)
}

# 4. Data Validation and Preprocessing --------------------------------------
# Verify required columns exist in all datasets
required_cols <- c("protein_name", "rs", "cis", "trans")
datasets <- list("Sex-combined" = df_all, "Female" = df_female, "Male" = df_male)

for (name in names(datasets)) {
    missing_cols <- setdiff(required_cols, names(datasets[[name]]))
    if (length(missing_cols) > 0) {
        stop(sprintf("Dataset '%s' is missing required columns: %s", 
                     name, paste(missing_cols, collapse = ", ")))
    }
}

# 5. Create Protein-RS Pair Identifiers -------------------------------------
# Function to create unique protein_name+rs combinations
create_protein_rs_pairs <- function(df, cis_value) {
    df %>%
        filter(cis == cis_value) %>%
        mutate(protein_rs_pair = paste(protein_name, rs, sep = "_")) %>%
        select(protein_rs_pair) %>%
        distinct() %>%
        pull(protein_rs_pair)
}

# 6. Generate Pair Sets for trans-SNPs (cis = 0) ----------------------------
cat("Processing trans-SNPs subsets...\n")
pairs_sex_combined_trans <- create_protein_rs_pairs(df_all, 0)
pairs_female_trans       <- create_protein_rs_pairs(df_female, 0)
pairs_male_trans         <- create_protein_rs_pairs(df_male, 0)

# Report subset sizes
cat(sprintf("  Sex-combined dataset (trans-SNPs): %d unique protein-rs pairs\n", length(pairs_sex_combined_trans)))
cat(sprintf("  Female dataset (trans-SNPs): %d unique protein-rs pairs\n", length(pairs_female_trans)))
cat(sprintf("  Male dataset (trans-SNPs): %d unique protein-rs pairs\n", length(pairs_male_trans)))

# 7. Generate Pair Sets for cis-SNPs (cis = 1) ------------------------------
cat("\nProcessing cis-SNPs subsets...\n")
pairs_sex_combined_cis <- create_protein_rs_pairs(df_all, 1)
pairs_female_cis       <- create_protein_rs_pairs(df_female, 1)
pairs_male_cis         <- create_protein_rs_pairs(df_male, 1)

# Report subset sizes
cat(sprintf("  Sex-combined dataset (cis-SNPs): %d unique protein-rs pairs\n", length(pairs_sex_combined_cis)))
cat(sprintf("  Female dataset (cis-SNPs): %d unique protein-rs pairs\n", length(pairs_female_cis)))
cat(sprintf("  Male dataset (cis-SNPs): %d unique protein-rs pairs\n", length(pairs_male_cis)))

# 8. Prepare Data for UpSet Visualization -----------------------------------
# Function to create binary matrix for UpSet plot
prepare_upset_data <- function(pairs_sex_combined, pairs_female, pairs_male, label) {
    # Combine all unique pairs
    all_unique_pairs <- unique(c(pairs_sex_combined, pairs_female, pairs_male))
    
    # Create binary membership matrix
    upset_matrix <- data.frame(
        protein_rs_pair = all_unique_pairs,
        `Sex-combined` = as.integer(all_unique_pairs %in% pairs_sex_combined),
        Female = as.integer(all_unique_pairs %in% pairs_female),
        Male   = as.integer(all_unique_pairs %in% pairs_male),
        stringsAsFactors = FALSE,
        check.names = FALSE
    )
    
    # Remove pairs not present in any dataset (should not happen, but safety check)
    upset_matrix <- upset_matrix[rowSums(upset_matrix[, -1]) > 0, ]
    
    return(upset_matrix)
}

# 9. Create UpSet Data Matrices ----------------------------------------------
upset_data_trans <- prepare_upset_data(pairs_sex_combined_trans, pairs_female_trans, pairs_male_trans, "trans-SNPs")
upset_data_cis <- prepare_upset_data(pairs_sex_combined_cis, pairs_female_cis, pairs_male_cis, "cis-SNPs")

# 10. Generate UpSet Plot for trans-SNPs ------------------------------------
cat("\nGenerating UpSet plot for trans-SNPs...\n")

# Determine appropriate number of intersections to show
n_intersections_trans <- min(40, length(unique(apply(upset_data_trans[, -1], 1, paste, collapse=""))))

pdf(file.path(base_dir, "upset_plot_trans_SNPs_protein_rs_pairs.pdf"), width = 12, height = 8)
upset(
    upset_data_trans[, -1],  # Exclude the pair identifier column
    sets = c("Sex-combined", "Female", "Male"),
    sets.bar.color = "#56B4E9",
    order.by = "freq",
    decreasing = TRUE,
    empty.intersections = NULL,
    nintersects = n_intersections_trans,
    mainbar.y.label = "Intersection of putative causal pQTLs",
    sets.x.label = "Number of putative causal pQTLs",
    text.scale = c(1.5, 1.3, 1.3, 1.3, 1.5, 1.3),
    mb.ratio = c(0.65, 0.35),
    point.size = 3.5,
    line.size = 1.2,
    shade.color = "gray88",
    matrix.color = "gray30",
    main.bar.color = "gray30",
    set_size.angles = 0,
    set_size.show = TRUE,
    set_size.scale_max = ceiling(max(colSums(upset_data_trans[, -1])) * 1.1)
)
dev.off()
cat("  Saved: upset_plot_trans_SNPs_protein_rs_pairs.pdf\n")

# 11. Generate UpSet Plot for cis-SNPs --------------------------------------
cat("\nGenerating UpSet plot for cis-SNPs...\n")

# Determine appropriate number of intersections to show
n_intersections_cis <- min(40, length(unique(apply(upset_data_cis[, -1], 1, paste, collapse=""))))

pdf(file.path(base_dir, "upset_plot_cis_SNPs_protein_rs_pairs.pdf"), width = 12, height = 8)
upset(
    upset_data_cis[, -1],  # Exclude the pair identifier column
    sets = c("Sex-combined", "Female", "Male"),
    sets.bar.color = "#E69F00",
    order.by = "freq",
    decreasing = TRUE,
    empty.intersections = NULL,
    nintersects = n_intersections_cis,
    mainbar.y.label = "Intersection of putative causal pQTLs",
    sets.x.label = "Number of putative causal pQTLs",
    text.scale = c(1.5, 1.3, 1.3, 1.3, 1.5, 1.3),
    mb.ratio = c(0.65, 0.35),
    point.size = 3.5,
    line.size = 1.2,
    shade.color = "gray88",
    matrix.color = "gray30",
    main.bar.color = "gray30",
    set_size.angles = 0,
    set_size.show = TRUE,
    set_size.scale_max = ceiling(max(colSums(upset_data_cis[, -1])) * 1.1)
)
dev.off()
cat("  Saved: upset_plot_cis_SNPs_protein_rs_pairs.pdf\n")

# 12. Generate Summary Statistics --------------------------------------------
cat("\n" , strrep("=", 60), "\n")
cat("SUMMARY STATISTICS\n")
cat(strrep("=", 60), "\n")

# Function to calculate intersection statistics
calculate_intersection_stats <- function(upset_data, snp_type_label) {
    cat(sprintf("\n%s:\n", snp_type_label))
    cat(strrep("-", 40), "\n")
    
    # Calculate intersection counts
    only_sex_combined <- sum(upset_data$`Sex-combined` == 1 & upset_data$Female == 0 & upset_data$Male == 0)
    only_female       <- sum(upset_data$`Sex-combined` == 0 & upset_data$Female == 1 & upset_data$Male == 0)
    only_male         <- sum(upset_data$`Sex-combined` == 0 & upset_data$Female == 0 & upset_data$Male == 1)
    sex_combined_female <- sum(upset_data$`Sex-combined` == 1 & upset_data$Female == 1 & upset_data$Male == 0)
    sex_combined_male   <- sum(upset_data$`Sex-combined` == 1 & upset_data$Female == 0 & upset_data$Male == 1)
    female_male       <- sum(upset_data$`Sex-combined` == 0 & upset_data$Female == 1 & upset_data$Male == 1)
    all_three         <- sum(upset_data$`Sex-combined` == 1 & upset_data$Female == 1 & upset_data$Male == 1)
    
    cat(sprintf("  Unique to Sex-combined dataset: %d\n", only_sex_combined))
    cat(sprintf("  Unique to Female dataset: %d\n", only_female))
    cat(sprintf("  Unique to Male dataset: %d\n", only_male))
    cat(sprintf("  Sex-combined ∩ Female only: %d\n", sex_combined_female))
    cat(sprintf("  Sex-combined ∩ Male only: %d\n", sex_combined_male))
    cat(sprintf("  Female ∩ Male only: %d\n", female_male))
    cat(sprintf("  All three datasets: %d\n", all_three))
    cat(sprintf("  Total unique pairs: %d\n", nrow(upset_data)))
    
    # Calculate percentages
    total <- nrow(upset_data)
    cat(sprintf("\n  Percentage breakdown:\n"))
    cat(sprintf("    Unique to Sex-combined: %.1f%%\n", 100 * only_sex_combined / total))
    cat(sprintf("    Unique to Female: %.1f%%\n", 100 * only_female / total))
    cat(sprintf("    Unique to Male: %.1f%%\n", 100 * only_male / total))
    cat(sprintf("    Sex-combined ∩ Female: %.1f%%\n", 100 * sex_combined_female / total))
    cat(sprintf("    Sex-combined ∩ Male: %.1f%%\n", 100 * sex_combined_male / total))
    cat(sprintf("    Female ∩ Male: %.1f%%\n", 100 * female_male / total))
    cat(sprintf("    All three: %.1f%%\n", 100 * all_three / total))
    
    # Return statistics as a list for potential further use
    return(list(
        only_sex_combined = only_sex_combined,
        only_female = only_female,
        only_male = only_male,
        sex_combined_female = sex_combined_female,
        sex_combined_male = sex_combined_male,
        female_male = female_male,
        all_three = all_three,
        total = total
    ))
}

stats_trans <- calculate_intersection_stats(upset_data_trans, "TRANS-SNPs INTERSECTIONS")
stats_cis <- calculate_intersection_stats(upset_data_cis, "CIS-SNPs INTERSECTIONS")

# 13. Export Intersection Data for Further Analysis -------------------------
cat("\n" , strrep("=", 60), "\n")
cat("EXPORTING INTERSECTION DATA\n")
cat(strrep("=", 60), "\n")

# Save the upset data matrices for potential further analysis
write_csv(upset_data_trans, file.path(base_dir, "upset_data_trans_SNPs_protein_rs_pairs.csv"))
write_csv(upset_data_cis, file.path(base_dir, "upset_data_cis_SNPs_protein_rs_pairs.csv"))
cat("  Saved: upset_data_trans_SNPs_protein_rs_pairs.csv\n")
cat("  Saved: upset_data_cis_SNPs_protein_rs_pairs.csv\n")

# 14. Create Detailed Intersection Report -----------------------------------
# Create a summary table of all statistics
summary_table <- data.frame(
    Intersection = c("Unique to Sex-combined", "Unique to Female", "Unique to Male",
                     "Sex-combined ∩ Female", "Sex-combined ∩ Male", "Female ∩ Male", 
                     "All Three", "Total"),
    `Trans_SNPs_Count` = c(stats_trans$only_sex_combined, stats_trans$only_female, stats_trans$only_male,
                           stats_trans$sex_combined_female, stats_trans$sex_combined_male, stats_trans$female_male,
                           stats_trans$all_three, stats_trans$total),
    `Trans_SNPs_Percent` = sprintf("%.1f%%", 100 * c(stats_trans$only_sex_combined, stats_trans$only_female, 
                                                     stats_trans$only_male, stats_trans$sex_combined_female, 
                                                     stats_trans$sex_combined_male, stats_trans$female_male,
                                                     stats_trans$all_three, stats_trans$total) / stats_trans$total),
    `Cis_SNPs_Count` = c(stats_cis$only_sex_combined, stats_cis$only_female, stats_cis$only_male,
                         stats_cis$sex_combined_female, stats_cis$sex_combined_male, stats_cis$female_male,
                         stats_cis$all_three, stats_cis$total),
    `Cis_SNPs_Percent` = sprintf("%.1f%%", 100 * c(stats_cis$only_sex_combined, stats_cis$only_female, 
                                                   stats_cis$only_male, stats_cis$sex_combined_female, 
                                                   stats_cis$sex_combined_male, stats_cis$female_male,
                                                   stats_cis$all_three, stats_cis$total) / stats_cis$total),
    check.names = FALSE
)

write_csv(summary_table, file.path(base_dir, "intersection_summary_statistics.csv"))
cat("  Saved: intersection_summary_statistics.csv\n")

# 15. Generate Venn Diagram as Alternative Visualization --------------------
if (require("VennDiagram")) {
    cat("\nGenerating Venn diagrams as alternative visualization...\n")
    
    # Function to create Venn diagram
    create_venn <- function(pairs_sex_combined, pairs_female, pairs_male, snp_type, color_scheme) {
        venn.plot <- venn.diagram(
            x = list(
                `Sex-combined` = pairs_sex_combined,
                Female = pairs_female,
                Male = pairs_male
            ),
            filename = NULL,
            col = "black",
            fill = color_scheme,
            alpha = 0.5,
            cex = 1.5,
            cat.cex = 1.5,
            cat.pos = c(-20, 20, 180),
            cat.dist = c(0.05, 0.05, 0.025),
            main = sprintf("Protein-RS Pairs (%s)", snp_type),
            main.cex = 2
        )
        
        # Save to file
        pdf(file.path(base_dir, sprintf("venn_diagram_%s_protein_rs_pairs.pdf", gsub("-", "_", snp_type))), 
            width = 8, height = 8)
        grid::grid.draw(venn.plot)
        dev.off()
    }
    
    # Create Venn diagrams
    create_venn(pairs_sex_combined_trans, pairs_female_trans, pairs_male_trans, "trans-SNPs", 
                c("#56B4E9", "#009E73", "#F0E442"))
    create_venn(pairs_sex_combined_cis, pairs_female_cis, pairs_male_cis, "cis-SNPs", 
                c("#E69F00", "#CC79A7", "#0072B2"))
    
    cat("  Saved: venn_diagram_trans_SNPs_protein_rs_pairs.pdf\n")
    cat("  Saved: venn_diagram_cis_SNPs_protein_rs_pairs.pdf\n")
} else {
    cat("\nNote: Install 'VennDiagram' package for additional Venn diagram visualizations.\n")
}

# 16. Export Lists of Unique and Shared Pairs -------------------------------
cat("\nExporting detailed pair lists...\n")

# Function to export specific intersection pairs
export_intersection_pairs <- function(upset_data, snp_type) {
    base_name <- sprintf("protein_rs_pairs_%s", gsub("-", "_", snp_type))
    
    # Extract pairs for each intersection type
    pairs_only_sex_combined <- upset_data[upset_data$`Sex-combined` == 1 & upset_data$Female == 0 & upset_data$Male == 0, "protein_rs_pair"]
    pairs_only_female <- upset_data[upset_data$`Sex-combined` == 0 & upset_data$Female == 1 & upset_data$Male == 0, "protein_rs_pair"]
    pairs_only_male <- upset_data[upset_data$`Sex-combined` == 0 & upset_data$Female == 0 & upset_data$Male == 1, "protein_rs_pair"]
    pairs_all_three <- upset_data[upset_data$`Sex-combined` == 1 & upset_data$Female == 1 & upset_data$Male == 1, "protein_rs_pair"]
    
    # Write to files if non-empty
    if (length(pairs_only_sex_combined) > 0) {
        write_lines(pairs_only_sex_combined, file.path(base_dir, paste0(base_name, "_unique_to_sex_combined.txt")))
    }
    if (length(pairs_only_female) > 0) {
        write_lines(pairs_only_female, file.path(base_dir, paste0(base_name, "_unique_to_female.txt")))
    }
    if (length(pairs_only_male) > 0) {
        write_lines(pairs_only_male, file.path(base_dir, paste0(base_name, "_unique_to_male.txt")))
    }
    if (length(pairs_all_three) > 0) {
        write_lines(pairs_all_three, file.path(base_dir, paste0(base_name, "_shared_all_three.txt")))
    }
}

export_intersection_pairs(upset_data_trans, "trans_SNPs")
export_intersection_pairs(upset_data_cis, "cis_SNPs")
cat("  Detailed pair lists exported.\n")

cat("\n" , strrep("=", 60), "\n")
cat("ANALYSIS COMPLETE\n")
cat(strrep("=", 60), "\n")
cat("\nAll output files have been saved to:\n")
cat(sprintf("  %s\n", base_dir))