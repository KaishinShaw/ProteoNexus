library(data.table)
library(dplyr)

PATH <- "/gpfs/chencao/ysbioinfor/project/proteohubProject/02_data/GWAS_fam/"
GENO_PATH <- "/gpfs/chencao/ysbioinfor/Datasets/ukb/geno/EUR_protein/hm3/"
# SEX <- c("all", "female", "male")
SEX <- c( "female", "male")
for (sex in SEX){
  
  pro_dat <- readRDS(paste0(PATH, "pro_resid_", sex, ".rds"))
  fam <- fread(paste0(GENO_PATH, sex, "/merge.fam"))
  fam[[6]] <- NULL
  fam_pro <- cbind.data.frame(fam, pro_dat)
  fwrite(fam_pro, file = paste0(GENO_PATH, sex, "/merge.fam"), 
         sep = " ")
}