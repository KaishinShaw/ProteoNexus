
PATH <- "/gpfs/chencao/ysbioinfor/project/proteohubProject"
GENO_PATH <- "/gpfs/chencao/ysbioinfor/Datasets/ukb/geno/EUR_protein/hm3/"



sex = "female"


dis_dat <- readRDS(paste0(PATH, "/02_data/out_pheno/disease_trait_", sex, "_dat.rds"))
cov_dat <- readRDS(paste0(PATH, "/02_data/cov/cov_", sex, ".rds"))
dis_resid <- alply(c(1: ncol(dis_dat)), 1, function(dis_s){
  
  dd <- dis_dat[, dis_s]
  dis_r <- rep(NA, length(dd))
  if(!all(is.na(dd))){

    dis_df <- data.frame(dd, cov_dat[, -1])
    dis_r[!is.na(dd)] <- lm(dd ~ ., data = dis_df)$residuals 
  }
  return(dis_r)
}) %>% do.call("cbind", .)


fam <- fread(paste0(GENO_PATH, sex, "/merge1.fam")) %>% as.data.frame()
# fam[, 6] <- NULL
fam_dis <- cbind.data.frame(fam[, 1:5], dis_resid)
fwrite(fam_dis, file = paste0(GENO_PATH, sex, "/merge.fam"), na = "-9", 
       sep = "\t", col.names = F)



