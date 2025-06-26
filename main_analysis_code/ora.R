#############################
## Perform over representation analysis

# Load packages 
library(clusterProfiler)
library(ReactomePA)
library(DOSE)
library(org.Hs.eg.db)
library(dplyr)
library(plyr)
library(ggplot2)
library(ggprism)
library(gground)

# Set parameters
N_PATHWAY <- 5

# Function 1: perform ORA 
ora_run <- function(pro_symbol){
  
  ## transfer gene symbol
  gene_symbol <- toupper(pro_symbol)
  gene_id <- bitr(gene_symbol, fromType = "SYMBOL",
                  toType = c("ENTREZID"), OrgDb = org.Hs.eg.db)
  ## GO 
  ora_res_GO <- enrichGO(gene = gene_id[, 1], keyType = 'SYMBOL',
                         OrgDb = org.Hs.eg.db, ont = "ALL",
                         pAdjustMethod = "BH", pvalueCutoff  = 1,
                         readable = TRUE)@result
  ## KEGG
  ora_res_KEGG <- enrichKEGG(gene = gene_id[, 2], organism = "hsa",
                             keyType = "ncbi-geneid",
                             pvalueCutoff = 1)@result
  ora_res_KEGG$geneID <- aaply(ora_res_KEGG$geneID, 1, function(ss){
    
    gene_eid <- strsplit(ss, "\\/") %>% unlist 
    gene_sym <- gene_id[gene_id$ENTREZID %in% gene_eid, 1, drop = T] %>% 
      paste0(collapse ="/")
    return(gene_sym)
  })
  
  ora_res <- list(GO = ora_res_GO, KEGG = ora_res_KEGG)
  
  return(ora_res)
}

# Function 2: 
ora_plt <- function(ora_res){
  
  ## select pathways
  use_pathway <- group_by(ora_res[['GO']], ONTOLOGY) %>%
    top_n(N_PATHWAY, wt = -p.adjust) %>%
    group_by(p.adjust) %>%
    top_n(1, wt = Count) %>%
    rbind(
      top_n(ora_res[['KEGG']], N_PATHWAY, -p.adjust) %>%
        group_by(p.adjust) %>%
        top_n(1, wt = Count) %>%
        mutate(ONTOLOGY = 'KEGG')
    ) %>%
    ungroup() %>%
    mutate(ONTOLOGY = factor(ONTOLOGY, 
                             levels = rev(c('BP', 'CC', 'MF', 'KEGG')))) %>%
    dplyr::arrange(ONTOLOGY, p.adjust) %>%
    mutate(Description = factor(Description, levels = Description)) %>%
    tibble::rowid_to_column('index')

  gene_out <- strsplit(use_pathway$geneID, "\\/") %>%
    aaply(., 1, function(x){
      if(length(x) > 6) {
        out_x <- sample(x, 6) %>% paste0(collapse = "/") %>% paste0("/...(", length(x) -6, ")") 
      } else {
        out_x <- paste0(x, collapse = "/")
      }
      return(out_x)
  })
  width <- 0.5
  xaxis_max <- max(-log10(use_pathway$p.adjust)) + 1
  rect.data <- group_by(use_pathway, ONTOLOGY) %>%
    reframe(n = n()) %>%
    ungroup() %>%
    mutate(xmin = -3 * width, xmax = -2 * width,
           ymax = cumsum(n), ymin = lag(ymax, default = 0) + 0.6, ymax = ymax + 0.4)
  pal <- c( '#fbb05b', '#acd372', '#ed6ca4','#7bc4e2')
  ## plot
  plt <- ggplot(use_pathway, aes(-log10(p.adjust), y = index, fill = ONTOLOGY)) +
    geom_round_col(aes(y = Description), width = 0.6, alpha = 0.8) +
    geom_text(aes(x = 0.05, label = Description), hjust = 0, size = 5) +
    geom_text(aes(x = 0.1, label = gene_out, colour = ONTOLOGY), 
              hjust = 0, vjust = 3, size = 3.5, fontface = 'italic', 
              show.legend = FALSE) +
    geom_point(aes(x = -width, size = Count), shape = 21) +
    geom_text(aes(x = -width, label = Count)) +
    scale_size_continuous(name = 'Count', range = c(2, 8)) +
    geom_round_rect(aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = ONTOLOGY),
                    data = rect.data, radius = unit(2, 'mm'), inherit.aes = FALSE) +
    geom_text(aes(x = (xmin + xmax) / 2, y = (ymin + ymax) / 2, label = ONTOLOGY),
              data = rect.data, angle = 90, hjust = 0.5,
              vjust = 0.5, inherit.aes = FALSE) +
    geom_segment(aes(x = 0, y = 0, xend = xaxis_max, yend = 0),
                 linewidth = 1.5, inherit.aes = FALSE) +
    labs(y = NULL, x = expression("-log"[10]*"FDR")) +
    scale_fill_manual(name = 'Category', values = pal) +
    scale_colour_manual(values = pal) +
    scale_x_continuous(breaks = seq(0, xaxis_max, 2), expand = expansion(c(0, 0))) +
    theme_prism() +
    theme(axis.text.y = element_blank(), 
          axis.line = element_blank(),
          axis.ticks.y = element_blank(),
          legend.title = element_text())
  return(plt)
}
