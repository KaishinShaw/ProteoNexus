library(ggplot2)
library(patchwork)

tag_plot <- function(p, tag, size = 18){
    p +
        labs(tag = tag) +
        theme(
            plot.tag        = element_text(face = "bold", size = size),
            plot.tag.position = c(0, 1)
        )
}

pA <- tag_plot(compare_cis_trans_pqtl_counts,      "A")
pB <- tag_plot(plot_protein_cis_trans_count_distribution, "B")

pC <- tag_plot(mpd_pathway,  "C")
pD <- tag_plot(mpd_protein,  "D")

pE <- tag_plot(epd_pathway,  "E")
pF <- tag_plot(epd_protein,  "F")

pG <- tag_plot(gpd_pathway,  "G")
pH <- tag_plot(gpd_protein,  "H")

pI <- tag_plot(figure_AUC,   "I")

row1     <- pA / pB
mediation<- (pC | pD) / (pE | pF) / (pG | pH)

final_plotA <- (row1 | mediation) +
    plot_layout(widths = c(1, 2))

final_plotB <- final_plotA / pI

print(final_plotB)
