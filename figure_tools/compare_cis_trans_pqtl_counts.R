# ============================================================
# R Script for Publication‐Quality pQTL Counts Visualization
# Author: Hydraulik (modified)
# Date: 2025‐07‐30
# Description: 
#   1. Load significant pQTL datasets (combined, male, female)
#   2. Summarise cis vs. trans counts
#   3. Create an SCI‐standard bar chart with custom colours
#   4. Export figures in multiple formats
#   5. Generate and save a summary statistics table
# Modifications:
#   • Rename x-axis label from "Cohort" to "Dataset"
#   • Rename "Combined" group to "Sex-combined"
#   • Retain SCI-approved colours: cis = "#cbbe7b", trans = "#4a6d9f"
#   • Add complete border frame around the plot
# ============================================================

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

# 3. Summarise cis/trans Counts ---------------------------------------------
summ_all    <- df_all    %>% summarise(cis = sum(cis, na.rm = TRUE),
                                       trans = sum(trans, na.rm = TRUE))
summ_male   <- df_male   %>% summarise(cis = sum(cis, na.rm = TRUE),
                                       trans = sum(trans, na.rm = TRUE))
summ_female <- df_female %>% summarise(cis = sum(cis, na.rm = TRUE),
                                       trans = sum(trans, na.rm = TRUE))

# 4. Combine and Reshape for Plotting ---------------------------------------
plot_df <- bind_rows(
    "Sex-combined" = summ_all,
    "Male"         = summ_male,
    "Female"       = summ_female,
    .id = "Dataset"
) %>%
    pivot_longer(cols      = c(cis, trans),
                 names_to  = "Type",
                 values_to = "Count") %>%
    mutate(
        Dataset = factor(Dataset, levels = c("Sex-combined", "Male", "Female")),
        Type    = factor(Type, levels = c("cis", "trans"),
                         labels = c("cis-SNPs", "trans-SNPs"))
    )

# 5. Create Publication‐Quality Bar Chart -----------------------------------
# Define SCI-approved colour palette
custom_cols <- c("cis-SNPs"   = "#cbbe7b",
                 "trans-SNPs" = "#4a6d9f")

p <- ggplot(plot_df, aes(x = Dataset, y = Count, fill = Type)) +
    geom_bar(stat     = "identity",
             position = position_dodge(width = 0.7),
             width    = 0.6,
             colour   = "black",
             size     = 0.3) +
    geom_text(aes(label = format(Count, big.mark = ",", scientific = FALSE)),
              position = position_dodge(width = 0.7),
              vjust    = -0.5,
              size     = 3.5,
              colour   = "black") +
    scale_fill_manual(values = custom_cols, name = "pQTL Type") +
    scale_y_continuous(
        labels = function(x) format(x, big.mark = ",", scientific = FALSE),
        breaks = pretty_breaks(n = 5),
        expand = expansion(mult = c(0, 0.05))
    ) +
    labs(
        x = "Dataset",
        y = "Number of significant pQTLs"
    ) +
    theme_classic(base_size = 16) +
    theme(
        axis.title       = element_text(face = "bold"),
        axis.text        = element_text(size = 14),
        axis.ticks       = element_line(colour = "black"),
        axis.line        = element_line(colour = "black"),
        legend.position  = "bottom",
        legend.title     = element_text(face = "bold"),
        legend.text      = element_text(size = 10),
        legend.key.size  = unit(0.6, "cm"),
        panel.background = element_rect(fill = "white", colour = NA),
        panel.border     = element_rect(colour = "black", fill = NA, size = 0.8),  # 添加完整外框
        panel.grid       = element_blank(),
        plot.margin      = margin(10, 10, 10, 10)
    )

# Display the plot
print(p)

# 6. Export Figures for Journal Submission ----------------------------------
output_dir <- "publication_figures"
if (!dir.exists(output_dir)) dir.create(output_dir)

save_plot <- function(plot, filename, width, height, dpi = 300, device = NULL, ...) {
    ggsave(
        filename = file.path(output_dir, filename),
        plot     = plot,
        width    = width,
        height   = height,
        units    = "in",
        dpi      = dpi,
        device   = device,
        ...
    )
}

# Export in multiple formats to meet journal requirements
save_plot(p, "pQTL_counts_sex_combined.pdf",      width = 7,   height = 5)
save_plot(p, "pQTL_counts_300dpi.png",            width = 7,   height = 5, dpi = 300)
save_plot(p, "pQTL_counts_600dpi.png",            width = 7,   height = 5, dpi = 600)
save_plot(p, "pQTL_counts.tiff",                  width = 7,   height = 5, dpi = 300, device = "tiff", compression = "lzw")
save_plot(p + theme(legend.direction = "vertical", legend.position = "right"),
          "pQTL_counts_single_column.pdf",       width = 3.5, height = 4)
save_plot(p, "pQTL_counts_double_column.pdf",     width = 7.2, height = 5)
save_plot(p, "pQTL_counts.eps",                   width = 7,   height = 5, device = "eps")

# 7. Generate & Save Summary Statistics Table -------------------------------
summary_table <- plot_df %>%
    pivot_wider(names_from = Type, values_from = Count) %>%
    mutate(
        Total            = `cis-SNPs` + `trans-SNPs`,
        `cis-SNPs (%)`   = round(`cis-SNPs`   / Total * 100, 1),
        `trans-SNPs (%)` = round(`trans-SNPs` / Total * 100, 1)
    ) %>%
    select(Dataset, `cis-SNPs`, `trans-SNPs`, Total, `cis-SNPs (%)`, `trans-SNPs (%)`)

cat("\n=== pQTL Summary Statistics by Dataset ===\n")
print(summary_table, n = Inf)

# Export summary as CSV
write_csv(summary_table, file.path(output_dir, "pQTL_summary_statistics.csv"))

cat("\nAll figures and summary table have been saved to:", output_dir, "\n")