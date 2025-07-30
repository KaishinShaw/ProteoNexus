# -------------------------------------------------------------------
# 1. Load required libraries
# -------------------------------------------------------------------
# install.packages(c("readr","dplyr","ggplot2","RColorBrewer","scales"))
library(readr)       # fast TSV import
library(dplyr)       # data manipulation
library(ggplot2)     # plotting
library(RColorBrewer)# professional palettes
library(scales)      # percent formatting

# -------------------------------------------------------------------
# 2. Define file paths and read data
# -------------------------------------------------------------------
mpd_path <- "C:/Users/shaok/Desktop/ProteoNexus_Citation/GRABBING/M-P-D_annotated.tsv"
epd_path <- "C:/Users/shaok/Desktop/ProteoNexus_Citation/GRABBING/E-P-D_annotated.tsv"

mpd <- read_tsv(mpd_path, col_types = cols()) 
epd <- read_tsv(epd_path, col_types = cols())

# -------------------------------------------------------------------
# 3. Pre‐process: set factor levels for consistent ordering
# -------------------------------------------------------------------
# Adjust the gender levels to your desired order (e.g. Male, Female, All)
gender_levels <- c("male", "female", "all")

# Automatically pick up the four Category levels in alphabetical order
mpd <- mpd %>%
    mutate(
        gender   = factor(gender,   levels = gender_levels),
        Category = factor(Category, levels = sort(unique(Category)))
    )

epd <- epd %>%
    mutate(
        gender   = factor(gender,   levels = gender_levels),
        Category = factor(Category, levels = sort(unique(Category)))
    )

# -------------------------------------------------------------------
# 4. A helper function to create a percent-stacked bar chart
# -------------------------------------------------------------------
create_percent_stacked_bar <- function(df, dataset_name, palette = "Dark2") {
    ggplot(df, aes(x = gender, fill = Category)) +
        geom_bar(position = "fill", width = 0.7, color = "black") +
        scale_y_continuous(
            labels = percent_format(accuracy = 1),
            expand = c(0, 0)
        ) +
        scale_fill_brewer(palette = palette) +
        labs(
            title = paste0(dataset_name, " Dataset"),
            x     = "Gender",
            y     = "Proportion (%)",
            fill  = "Category"
        ) +
        theme_classic(base_size = 14) +
        theme(
            plot.title        = element_text(face = "bold", hjust = 0.5),
            axis.title        = element_text(face = "bold"),
            axis.text         = element_text(color = "black"),
            legend.title      = element_text(face = "bold"),
            legend.text       = element_text(size = 12)
        )
}

# -------------------------------------------------------------------
# 5. Generate the two plots
# -------------------------------------------------------------------
plot_mpd <- create_percent_stacked_bar(mpd, "M-P-D")
plot_epd <- create_percent_stacked_bar(epd, "E-P-D")

# Display in RStudio’s plot pane (or equivalent)
print(plot_mpd)
print(plot_epd)

# -------------------------------------------------------------------
# 6. (Optional) Save high-resolution figures for journal submission
# -------------------------------------------------------------------
ggsave(
    filename = "mpd_percent_stacked_bar.png",
    plot     = plot_mpd,
    width    = 6,   # inches
    height   = 4,   # inches
    dpi      = 300  # publication quality
)

ggsave(
    filename = "epd_percent_stacked_bar.png",
    plot     = plot_epd,
    width    = 6,
    height   = 4,
    dpi      = 300
)