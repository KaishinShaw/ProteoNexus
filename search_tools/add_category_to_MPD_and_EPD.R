# Load required packages
library(readr)    # fast file I/O
library(dplyr)    # data manipulation

# Define file paths
base_dir        <- "C:/Users/shaok/Desktop/ProteoNexus_Citation/GRABBING"
mpd_file        <- file.path(base_dir, "M-P-D_merged.tsv")
measurement_file<- file.path(base_dir, "measurement.csv")
epd_file        <- file.path(base_dir, "E-P-D_merged.tsv")
exposure_file   <- file.path(base_dir, "exposure.csv")

# Read in the data
MPD         <- read_tsv(mpd_file,         col_types = cols())
measurement <- read_csv(measurement_file, col_types = cols())
EPD         <- read_tsv(epd_file,         col_types = cols())
exposure    <- read_csv(exposure_file,    col_types = cols())

# -------------------------------------------------------
# 1) Annotate MPD:
#    - add a new column 'Category'
#    - for each row in MPD, match PHM to measurement$PFID
#      and pull in measurement$Category
#    - then move 'Category' to be the first column
# -------------------------------------------------------
MPD <- MPD %>%
    mutate(
        Category = measurement$Category[ match(PHM, measurement$PFID) ]
    ) %>%
    select(Category, everything())

# -------------------------------------------------------
# 2) Annotate EPD:
#    - add a new column 'Category'
#    - for each row in EPD, match PHE to exposure$PHID
#      and pull in exposure$Category
#    - then move 'Category' to be the first column
# -------------------------------------------------------
EPD <- EPD %>%
    mutate(
        Category = exposure$Category[ match(PHE, exposure$PHID) ]
    ) %>%
    select(Category, everything())

# (Optional) Inspect the results
glimpse(MPD)
glimpse(EPD)

# (Optional) Write back to disk
write_tsv(MPD, file.path(base_dir, "M-P-D_annotated.tsv"))
write_tsv(EPD, file.path(base_dir, "E-P-D_annotated.tsv"))