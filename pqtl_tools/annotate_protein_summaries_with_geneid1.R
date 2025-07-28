# Load necessary libraries
library(dplyr)

# Define the directory containing the CSV files
input_dir <- "C:/Users/shaok/Desktop/ProteoNexus_Citation/GRABBING"

# List the target files
file_names <- c(
    "merged_sig_summ_all2.csv",
    "merged_sig_summ_female2.csv",
    "merged_sig_summ_male2.csv"
)

# Helper function to process one file
process_file <- function(file_path, gene_tbl) {
    # Read the CSV into a data frame
    df <- read.csv(file_path, stringsAsFactors = FALSE)
    
    # Perform a case‐insensitive match of protein_name to genesymbol
    # and pull in the corresponding geneid1 from gene_tbl
    idx <- match(
        tolower(df$protein_name),
        tolower(gene_tbl$genesymbol)
    )
    
    # Create the new column geneid1
    df$geneid1 <- gene_tbl$geneid1[idx]
    
    return(df)
}

# Process each file and optionally write back to disk
processed_dfs <- list()

for (fn in file_names) {
    input_path  <- file.path(input_dir, fn)
    
    # Process and store in a list
    processed_df <- process_file(input_path, gene_tbl)
    processed_dfs[[fn]] <- processed_df
    
    # Write the augmented data frame to a new CSV
    output_path <- sub("\\.csv$", "_with_geneid1.csv", input_path)
    write.csv(processed_df, output_path, row.names = FALSE)
    
    message("Processed and wrote: ", output_path)
}

# If you want to attach them to the global environment:
# names(processed_dfs) will be the original file names
list2env(processed_dfs, envir = .GlobalEnv)