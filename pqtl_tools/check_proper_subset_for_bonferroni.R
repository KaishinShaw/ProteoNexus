# Load required libraries
library(dplyr)
library(tidyverse)

# Function to check if df1 is a proper subset of df2
check_proper_subset <- function(df1, df2, file1_name, file2_name) {
    
    cat("\n", paste(rep("=", 80), collapse = ""), "\n")
    cat("Checking if", file1_name, "is a proper subset of", file2_name, "\n")
    cat(paste(rep("-", 80), collapse = ""), "\n")
    
    # Check column names consistency
    cols1 <- colnames(df1)
    cols2 <- colnames(df2)
    
    # Find common columns
    common_cols <- intersect(cols1, cols2)
    
    # Check if all columns from df1 exist in df2
    if (!all(cols1 %in% cols2)) {
        missing_cols <- setdiff(cols1, cols2)
        cat("WARNING: The following columns from", file1_name, "are not in", file2_name, ":\n")
        cat(paste(missing_cols, collapse = ", "), "\n")
        cat("Cannot determine proper subset relationship due to column mismatch.\n")
        return(FALSE)
    }
    
    # Report dataset dimensions
    cat("Dataset dimensions:\n")
    cat(sprintf("  %s: %d rows x %d columns\n", file1_name, nrow(df1), ncol(df1)))
    cat(sprintf("  %s: %d rows x %d columns\n", file2_name, nrow(df2), ncol(df2)))
    
    # Reorder columns to match for comparison
    df1_ordered <- df1[, cols1]
    df2_subset <- df2[, cols1]
    
    # Convert all columns to character for consistent comparison
    df1_char <- df1_ordered %>% mutate(across(everything(), as.character))
    df2_char <- df2_subset %>% mutate(across(everything(), as.character))
    
    # Create a unique row identifier for each row
    df1_rows <- df1_char %>%
        mutate(row_id = row_number()) %>%
        unite("row_signature", -row_id, sep = "|__|", remove = FALSE)
    
    df2_rows <- df2_char %>%
        unite("row_signature", everything(), sep = "|__|", remove = FALSE)
    
    # Check if each row from df1 exists in df2
    matches <- df1_rows$row_signature %in% df2_rows$row_signature
    
    # Calculate statistics
    n_matched <- sum(matches)
    n_unmatched <- sum(!matches)
    match_percentage <- (n_matched / nrow(df1)) * 100
    
    # Report results
    cat("\nMatching Statistics:\n")
    cat(sprintf("  Rows in %s that exist in %s: %d (%.2f%%)\n", 
                file1_name, file2_name, n_matched, match_percentage))
    cat(sprintf("  Rows in %s that DO NOT exist in %s: %d (%.2f%%)\n", 
                file1_name, file2_name, n_unmatched, 100 - match_percentage))
    
    # Determine proper subset status
    is_proper_subset <- (n_matched == nrow(df1)) && (nrow(df1) < nrow(df2))
    is_subset <- (n_matched == nrow(df1))
    
    # Report subset status
    cat("\nSubset Analysis:\n")
    if (is_proper_subset) {
        cat(sprintf("  ✓ %s IS a PROPER SUBSET of %s\n", file1_name, file2_name))
        cat("  (All rows match and the first dataset is smaller)\n")
    } else if (is_subset && nrow(df1) == nrow(df2)) {
        cat(sprintf("  %s is EQUAL to %s (same rows, not a proper subset)\n", 
                    file1_name, file2_name))
    } else if (is_subset) {
        cat(sprintf("  %s IS a SUBSET of %s (but not specified if proper)\n", 
                    file1_name, file2_name))
    } else {
        cat(sprintf("  ✗ %s is NOT a subset of %s\n", file1_name, file2_name))
        
        # Show examples of non-matching rows if any exist
        if (n_unmatched > 0 && n_unmatched <= 10) {
            cat("\nExamples of non-matching rows (row numbers from", file1_name, "):\n")
            non_matching_indices <- which(!matches)
            cat(paste(non_matching_indices, collapse = ", "), "\n")
        } else if (n_unmatched > 10) {
            cat("\nFirst 10 non-matching row numbers from", file1_name, ":\n")
            non_matching_indices <- which(!matches)[1:10]
            cat(paste(non_matching_indices, collapse = ", "), "\n")
        }
    }
    
    return(list(
        is_subset = is_subset,
        is_proper_subset = is_proper_subset,
        n_matched = n_matched,
        n_unmatched = n_unmatched,
        match_percentage = match_percentage
    ))
}

# Main execution
main <- function() {
    
    cat("Starting Proper Subset Analysis\n")
    cat("===============================\n")
    
    # Define file pairs to check
    file_pairs <- list(
        list(pip080 = "pip080_merged_sig_summ_male2_with_geneid1_filtered_annotated.csv",
             merged = "merged_sig_summ_male2_with_geneid1_filtered_annotated.csv",
             label = "male"),
        list(pip080 = "pip080_merged_sig_summ_female2_with_geneid1_filtered_annotated.csv",
             merged = "merged_sig_summ_female2_with_geneid1_filtered_annotated.csv",
             label = "female"),
        list(pip080 = "pip080_merged_sig_summ_all2_with_geneid1_filtered_annotated.csv",
             merged = "merged_sig_summ_all2_with_geneid1_filtered_annotated.csv",
             label = "all")
    )
    
    # Store results
    results <- list()
    
    # Process each file pair
    for (pair in file_pairs) {
        
        # Read the CSV files
        tryCatch({
            cat("\nReading", pair$label, "dataset files...\n")
            
            # Read pip080 file
            df_pip080 <- read.csv(pair$pip080, 
                                  stringsAsFactors = FALSE, 
                                  check.names = FALSE)
            
            # Read merged file
            df_merged <- read.csv(pair$merged, 
                                  stringsAsFactors = FALSE, 
                                  check.names = FALSE)
            
            # Perform subset check
            result <- check_proper_subset(df_pip080, df_merged, 
                                          pair$pip080, pair$merged)
            
            # Store result
            results[[pair$label]] <- result
            
        }, error = function(e) {
            cat("ERROR processing", pair$label, "files:\n")
            cat(as.character(e), "\n")
            results[[pair$label]] <- NULL
        })
    }
    
    # Summary report
    cat("\n", paste(rep("=", 80), collapse = ""), "\n")
    cat("FINAL SUMMARY\n")
    cat(paste(rep("=", 80), collapse = ""), "\n\n")
    
    for (pair in file_pairs) {
        if (!is.null(results[[pair$label]])) {
            result <- results[[pair$label]]
            status <- ifelse(result$is_proper_subset, "YES ✓", "NO ✗")
            cat(sprintf("%-10s dataset: Proper Subset = %s (%.2f%% rows matched)\n", 
                        toupper(pair$label), status, result$match_percentage))
        } else {
            cat(sprintf("%-10s dataset: ERROR - Could not process files\n", 
                        toupper(pair$label)))
        }
    }
    
    cat("\n", paste(rep("=", 80), collapse = ""), "\n")
    
    return(results)
}

# Execute the analysis
results <- main()

# Optional: Save detailed results to an R object for further analysis
# saveRDS(results, "subset_analysis_results.rds")