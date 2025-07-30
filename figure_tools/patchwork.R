# Load necessary library
library(patchwork)

# ------------------------------------------------------------------------------
# Step 1: Create 'row1' by vertically stacking two plots
#   - Top:    compare_cis_trans_pqtl_counts
#   - Bottom: plot_protein_cis_trans_count_distribution
# ------------------------------------------------------------------------------
row1 <- compare_cis_trans_pqtl_counts /
    plot_protein_cis_trans_count_distribution

# ------------------------------------------------------------------------------
# Step 2: Create 'mpd' by horizontally aligning pathway and protein plots
#   - Left:  mpd_pathway
#   - Right: mpd_protein
# ------------------------------------------------------------------------------
mpd <- mpd_pathway | mpd_protein

# ------------------------------------------------------------------------------
# Step 3: Create 'epd' by horizontally aligning pathway and protein plots
#   - Left:  epd_pathway
#   - Right: epd_protein
# ------------------------------------------------------------------------------
epd <- epd_pathway | epd_protein

# ------------------------------------------------------------------------------
# Step 4: Create 'gpd' by horizontally aligning pathway and protein plots
#   - Left:  gpd_pathway
#   - Right: gpd_protein
# ------------------------------------------------------------------------------
gpd <- gpd_pathway | gpd_protein

# ------------------------------------------------------------------------------
# Step 5: Stack the three medium‐pathway assemblies vertically to form 'mediation'
#   - Top:    mpd
#   - Middle: epd
#   - Bottom: gpd
# ------------------------------------------------------------------------------
mediation <- mpd /
    epd /
    gpd

# ------------------------------------------------------------------------------
# Step 6: Compose the final layout by placing 'row1' on the left
#         and 'mediation' on the right
# ------------------------------------------------------------------------------
final_plotA <- (row1 | mediation) +
    plot_layout(widths = c(1, 2))

final_plotB <- final_plotA /
    figure_AUC

# Render the combined plot
print(final_plotB)

final_plotB <- final_plotB + 
    plot_annotation(tag_levels = 'A')

print(final_plotB)