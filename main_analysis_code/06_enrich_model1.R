#!/bin/bash

########################
## Enrichment for the significant proteins associated with exposure/measurement

# Load code
library(data.table)
source("/gpfs/chencao/ysbioinfor/project/proteohubProject/01_code/ora.R")
# enrich_plt <- ora_plt(ora_result)
# Set parameter
ALPHA <- 0.05
DATA_PATH <- "/gpfs/chencao/ysbioinfor/project/proteohubProject/02_data"
TRAIT <- read.table(paste0(DATA_PATH, "/model1_label.txt"))[, 1]
## 1:183
for (tt in 128: 150){
  
  trait_id = TRAIT[tt]
  OUT_PATH <- "/gpfs/chencao/ysbioinfor/project/proteohubProject/web_out/"
  OUT_PATH <- ifelse(grepl("PHE", trait_id), paste0(OUT_PATH, "E-P/"), paste0(OUT_PATH, "M-P/"))
  for (sex in c("all", "female", "male")){
    print(sex)
    # ORA analysis
    pro_result <- fread(paste0(OUT_PATH, trait_id, "/effect_size_", sex, ".tsv"))  %>% 
      as.data.frame()
    sig_pro <- pro_result$protein[pro_result$P_BH < ALPHA]
    if(is.na(sig_pro[1]) | length(sig_pro) < 6){
      
      go_col <- c("ONTOLOGY", "ID", "Description", "GeneRatio", "BgRatio", "pvalue", 
                  "p.adjust", "qvalue", "geneID", "Count") %>% matrix(., nrow = 1)
      kegg_col <- c("category", "subcategory", "ID", "Description", "GeneRatio",  
                    "BgRatio", "pvalue", "p.adjust", "qvalue", "geneID", "Count") %>% matrix(., nrow = 1)
      write.table(go_col, sep = "\t", quote = F, col.names = F, row.names = F,
                  file = paste0(OUT_PATH, trait_id, "/enrichment_GO_", sex, ".tsv"))
      write.table(kegg_col, sep = "\t", quote = F, col.names = F, row.names = F,
                  file = paste0(OUT_PATH, trait_id, "/enrichment_KEGG_", sex, ".tsv"))
    } else {
      
      ora_result <- ora_run(sig_pro)
      fwrite(ora_result[[1]] %>% filter(p.adjust<ALPHA), sep = "\t", 
             file = paste0(OUT_PATH, trait_id, "/enrichment_GO_", sex, ".tsv"))
      fwrite(ora_result[[2]] %>% filter(p.adjust<ALPHA), sep = "\t", 
             file = paste0(OUT_PATH, trait_id, "/enrichment_KEGG_", sex, ".tsv"))
      enrich_plt <- ora_plt(ora_result)
      ggsave(filename = paste0(OUT_PATH, trait_id, "/enrichment_plot_", sex, ".png"), 
             enrich_plt, width = 18, height = 12, units = "in", dpi = 320)
    }
  }
}


