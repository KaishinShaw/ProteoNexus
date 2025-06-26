library(dplyr)
library(plyr)

setwd("/gpfs/chencao/ysbioinfor/project/proteohubProject/01_code/check")

dis_label <- read.table("d_male.txt")[, 1]
ms_label <- read.table("m_label.txt")[, 1]

PROJ_PATH <- "/gpfs/chencao/ysbioinfor/project/proteohubProject/"
DATA_PATH <- paste0(PROJ_PATH, "02_data/")


sex <- "male"
SIG_FILE <- read.table(paste0(DATA_PATH, "/sig_total/emd_", sex, ".txt"), head = F) 
fac_trait <- unique(SIG_FILE[, 1])

# for (ff in 1: length(fac_trait)){}
miss <- alply(c(1: length(fac_trait)), 1, function(ff){
  
  dis_dir <- ifelse(grepl("PHM", fac_trait[ff]), 
                    paste0(PROJ_PATH, "/web_out/M-P-D/", fac_trait[ff]), 
                    paste0(PROJ_PATH, "/web_out/E-P-D/", fac_trait[ff]))
  # cat(dis_dir, "\n")
  dis_l <- list.files(dis_dir, pattern = paste0("*_", sex, ".tsv"))
  dis_trait <- SIG_FILE[SIG_FILE[, 1] ==  fac_trait[ff], 2] %>% 
    paste0(., "_", sex, ".tsv")
  if(length(dis_l) < length(dis_trait)){
    
    miss_s <- setdiff(dis_trait, dis_l) %>% substr(., 1, 8)
    miss_d <- cbind(rep(fac_trait[ff], length(miss_s)), miss_s)
  } else {
    
    miss_d <- rep(NA, 2)
  }
  return(miss_d)
}) %>% do.call("rbind", . ) %>% na.omit
fwrite(miss, file = paste0("miss_emd_", sex, ".txt"), 
       sep = "\t", col.names = F)

