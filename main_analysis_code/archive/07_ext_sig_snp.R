#!/bin/bash

########################
## Select significant SNPs for each protein

# Load packages
library(plyr)
library(dplyr)
library(data.table)


# Set parameter
P_BON <- 1e-8
PIP <- 0.85
PROJ_PATH <- "/gpfs/chencao/ysbioinfor/project/proteohubProject/"
DATA_PATH <- paste0(PROJ_PATH, "02_data/")
EXP_PATH <- paste0(PROJ_PATH, "/web_out/G-P/")

# Select significant pQTL data
pro_label <- read.table(paste0(DATA_PATH, "protein_label.txt"))[, 1]
for (sex in c("all", "female", "male")){

  sig_snp_dat <- alply(c(1: length(pro_label)), 1, function (ss){

    cat("Protein:", ss, "\n")
    sig_summstat <- fread(paste0(EXP_PATH, pro_label[ss], "/output/summ_", sex, "2.assoc.txt.gz")) %>%
      filter(p_wald < P_BON & pip_susie > PIP)
    fwrite(sig_summstat, file = paste0(EXP_PATH, pro_label[ss], "/output/sig_summ_", sex, "2.assoc.tsv"),
           sep = "\t")
    if (nrow(sig_summstat) != 0){

      pro_snp <- cbind(pro_label[ss], sig_summstat$rs)
      return(pro_snp)
    }
  }) %>% do.call("rbind", .)
  write.table(sig_snp_dat, row.names = F, col.names = F, quote = F,
              file = paste0(DATA_PATH, "/sig_genotype/pqtl_", sex, ".txt"))
  write.table(unique(sig_snp_dat[,2]), row.names = F, col.names = F, quote = F,
              file = paste0(DATA_PATH, "/sig_genotype/sig_snp_", sex, ".txt"))
}



# Select significant pQTL data
pro_label <- read.table(paste0(DATA_PATH, "protein_label.txt"))[, 1]
for (sex in c("all", "female", "male")){
  
  pro_count <- 0
  num_pqtl <- 0
  num_pqtl_out <- 0
  # sig_snp_dat <- alply(c(1: length(pro_label)), 1, function (ss){
  for (ss in 1: length(pro_label))  {
    
    cat("Protein:", ss, "\n")
    num_pqtl <- fread(paste0(EXP_PATH, pro_label[ss], "/output/sig_summ_", sex, "2.assoc.tsv")) %>%
      nrow
    if (num_pqtl != 0) pro_count <- pro_count + 1
    num_pqtl_out <- num_pqtl_out + num_pqtl
  }

    
    
    
    
    
    if (nrow(sig_summstat) != 0){
      
      pro_snp <- cbind(pro_label[ss], sig_summstat$rs)
      return(pro_snp)
    }
  }) %>% do.call("rbind", .)
  write.table(sig_snp_dat, row.names = F, col.names = F, quote = F,
              file = paste0(DATA_PATH, "/sig_genotype/pqtl_", sex, ".txt"))
  write.table(unique(sig_snp_dat[,2]), row.names = F, col.names = F, quote = F,
              file = paste0(DATA_PATH, "/sig_genotype/sig_snp_", sex, ".txt"))
}







