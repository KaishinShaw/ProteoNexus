# =========================================================================
# Stage 1: Setup - Loading Required Libraries
# =========================================================================
# 'readr' is for fast and friendly reading of rectangular data (like csv)
# 'dplyr' is for powerful and intuitive data manipulation
# 'ggplot2' is the primary library for creating graphics
# 'ggsci' provides a collection of high-quality color palettes from scientific journals

# If you haven't installed these packages, uncomment and run the following line:
# install.packages(c("readr", "dplyr", "ggplot2", "ggsci"))

library(readr)
library(dplyr)
library(ggplot2)
library(ggsci)


# =========================================================================
# Stage 2: Data Loading
# =========================================================================
# Define the file paths. 
# IMPORTANT: Ensure this path is correct for your system. R uses forward slashes '/'
# even on Windows, which is generally safer than backslashes '\'.

path_all    <- "C:/Users/shaok/Desktop/ProteoNexus_Citation/GRABBING/merged_sig_summ_all2_with_geneid1_filtered_annotated.csv"
path_female <- "C:/Users/shaok/Desktop/ProteoNexus_Citation/GRABBING/merged_sig_summ_female2_with_geneid1_filtered_annotated.csv"
path_male   <- "C:/Users/shaok/Desktop/ProteoNexus_Citation/GRABBING/merged_sig_summ_male2_with_geneid1_filtered_annotated.csv"

# Read the CSV files into data frames
# We use a tryCatch block to provide a helpful error message if a file is not found.
tryCatch({
    all_pqtl    <- read_csv(path_all)
    female_pqtl <- read_csv(path_female)
    male_pqtl   <- read_csv(path_male)
}, error = function(e) {
    stop("Error loading files. Please check if the file paths are correct and the files exist.\nOriginal error: ", e$message)
})


# =========================================================================
# Stage 3: Data Processing and Transformation
# =========================================================================
# The goal is to create a single, "tidy" data frame suitable for ggplot2.
# This data frame will have three columns: Cohort, pQTL_Type, and Count.

# Calculate the sums of 'cis' and 'trans' columns for each dataset
summary_all <- all_pqtl %>% 
    summarise(cis = sum(cis, na.rm = TRUE), trans = sum(trans, na.rm = TRUE))

summary_female <- female_pqtl %>% 
    summarise(cis = sum(cis, na.rm = TRUE), trans = sum(trans, na.rm = TRUE))

summary_male <- male_pqtl %>% 
    summarise(cis = sum(cis, na.rm = TRUE), trans = sum(trans, na.rm = TRUE))

# Combine the summaries into a single data frame
# We add a 'Cohort' column to identify the source of the data
plot_data <- bind_rows(
    "All"    = summary_all,
    "Female" = summary_female,
    "Male"   = summary_male,
    .id = "Cohort"
) %>%
    # Now, we pivot the data from a "wide" format to a "long" format.
    # This is the standard and most effective way to structure data for ggplot2.
    tidyr::pivot_longer(
        cols = c(cis, trans),
        names_to = "pQTL_Type",
        values_to = "Count"
    )

# To ensure the plot order is logical (All, Female, Male), we convert
# the 'Cohort' column to a factor with a specified order.
plot_data$Cohort <- factor(plot_data$Cohort, levels = c("All", "Female", "Male"))

# We also make 'pQTL_Type' a factor to control its order and labels in the plot.
plot_data$pQTL_Type <- factor(plot_data$pQTL_Type, levels = c("cis", "trans"), labels = c("Cis-pQTL", "Trans-pQTL"))


# =========================================================================
# Stage 4: Visualization with ggplot2
# =========================================================================
# Now we build the plot layer by layer for a publication-quality result.

pQTL_bar_chart <- ggplot(plot_data, aes(x = Cohort, y = Count, fill = pQTL_Type)) +
    
    # Add the bar layer. 
    # 'stat = "identity"' means the bar heights are the values from the 'Count' column.
    # 'position = "dodge"' places the bars for cis and trans next to each other.
    geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
    
    # Add text labels on top of each bar for clarity.
    # 'position_dodge' ensures the labels are centered above their respective bars.
    # 'vjust' adjusts the vertical position to be slightly above the bar.
    geom_text(aes(label = Count), 
              position = position_dodge(width = 0.8), 
              vjust = -0.5, 
              size = 3.5,
              family = "sans") + # Using a standard sans-serif font
    
    # Apply a professional color palette.
    # 'scale_fill_npg()' is from the 'ggsci' package and emulates the Nature Publishing Group palette.
    # Other excellent choices: scale_fill_brewer(palette = "Paired"), scale_fill_lancet()
    scale_fill_npg() +
    
    # Set the labels for the axes and the legend.
    labs(
        title = "Number of Significant pQTLs by Cohort and Type",
        subtitle = "Comparison across combined, female, and male cohorts",
        x = "Cohort",
        y = "Number of Significant pQTLs",
        fill = "pQTL Type" # This sets the legend title
    ) +
    
    # Apply a clean, classic theme suitable for scientific papers.
    # 'base_size' sets a baseline font size for all text elements.
    theme_classic(base_size = 14) +
    
    # Further theme customizations for a polished look.
    theme(
        # Axis text and titles
        axis.text.x = element_text(size = 12, color = "black", angle = 0, hjust = 0.5),
        axis.text.y = element_text(size = 12, color = "black"),
        axis.title.x = element_text(size = 14, face = "bold", margin = margin(t = 10)),
        axis.title.y = element_text(size = 14, face = "bold", margin = margin(r = 10)),
        
        # Axis lines
        axis.line = element_line(colour = "black", linewidth = 0.6),
        
        # Legend customization
        legend.position = "top",
        legend.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 11),
        
        # Plot title and subtitle
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 12, hjust = 0.5, margin = margin(b = 15)),
        
        # Remove grid lines for a cleaner look
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
    ) +
    
    # Expand the y-axis slightly to ensure the text labels on top of the bars fit well.
    scale_y_continuous(expand = expansion(mult = c(0, 0.1)))

# Display the plot in the RStudio Plots pane
print(pQTL_bar_chart)


# =========================================================================
# Stage 5: Exporting the Plot for Publication
# =========================================================================
# Use ggsave() to save the plot to a file.
# PDF is a vector format, ideal for publications as it can be scaled without losing quality.
# PNG or TIFF are high-quality raster formats.

# Saving as a PDF
ggsave(
    filename = "pQTL_counts_by_cohort.pdf",
    plot = pQTL_bar_chart,
    width = 8,  # inches
    height = 6, # inches
    device = "pdf"
)

# Saving as a high-resolution PNG (e.g., 300 DPI)
ggsave(
    filename = "pQTL_counts_by_cohort_300dpi.png",
    plot = pQTL_bar_chart,
    width = 8,
    height = 6,
    dpi = 300
)

# You can also save in TIFF format, often requested by journals
ggsave(
    filename = "pQTL_counts_by_cohort.tiff",
    plot = pQTL_bar_chart,
    width = 8,
    height = 6,
    dpi = 300
)