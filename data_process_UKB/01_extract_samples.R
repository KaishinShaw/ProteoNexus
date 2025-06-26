########################
## Extract EUR samples 

# Load packages
library(bigreadr)
library(dplyr)
library(plyr)


# Set parameters
PATH <- "/public/home/Datasets/ukb/pheno/04_olink"
PROJ_PATH <- "/public/home/biostat03/project/proteohubProject/"
TRAIT <- "/public/home/Datasets/ukb/pheno/03_trait/participant.csv.gz"


# Load data
covar_dat <- fread2(TRAIT, 
                    select = c("eid", 
                               "p31",                          ## Sex
                               "p21000_i0",                    ## Ethnic background | Instance 0
                               "p21000_i1",                    ## Ethnic background | Instance 1
                               "p21000_i2",                    ## Ethnic background | Instance 2
                               "p21000_i3",                    ## Ethnic background | Instance 3
                               "p22001",                       ## Genetic sex
                               "p22018",                       ## Genetic relatedness exclusions
                               "p22010",                       ## Recommended genomic analysis exclusions
                               "p22020",                       ## Used in genetic principal components 
                               "p54_i0",                       ## site
                               "p21022",                       ## Age
                               "p21001_i0",                    ## BMI
                               "p22000",                       ## UKB Array
                               paste0("p22009_a", 1:18)))      ## Genetic PC
pro_dat <- fread2(paste0(PATH, "/olink_instance_0.csv.gz"))
fam_file <- fread2("/public/home/Datasets/ukb/geno/plink/chr1.fam")[, 1]

# QC 
## Step1: select individuals with protein data
covar_pro_dat <- covar_dat[match(pro_dat$eid, covar_dat$eid), ]

## Step2: select individuals with EUR genetic background
eth_EUR <- c("White", "British", "Any other white background") 
cnd_i_EUR <- !(covar_pro_dat$p21000_i0 %in% eth_EUR | covar_pro_dat$p21000_i1 %in% eth_EUR |
  covar_pro_dat$p21000_i2 %in% eth_EUR | covar_pro_dat$p21000_i3 %in% eth_EUR)
cat("Exclusion: ", sum(cnd_i_EUR), " individuals are EUR in UKB-PPP\n")
cnd_ii <- covar_pro_dat$p22018 != ""
cat("Exclusion: ", sum(cnd_ii), " individuals are relatedness.\n")
cnd_iii <- covar_pro_dat$p31 != covar_pro_dat$p22001
cat("Exclusion: ", sum(cnd_iii), " individuals are gender incongruence.\n")
cnd_iv <- covar_pro_dat$p22010 != "" | !covar_pro_dat$eid %in% fam_file
cat("Exclusion: ", sum(cnd_iv), " individuals are genomic analysis exclusions.\n")
cnd_v <- covar_pro_dat$p22020 != "Yes" 
cat("Exclusion: ", sum(cnd_v), " individuals are without genetic PC.\n")
cnd_vi <- is.na(covar_pro_dat$p21001_i0)
cat("Exclusion: ", sum(cnd_vi), " individuals are without BMI.\n")
covar_pro_dat <- covar_pro_dat[!cnd_i_EUR & !cnd_ii & !cnd_iii & !cnd_iv & !cnd_v & !cnd_vi, ]

## Step3: select individuals with high quality protein data
pro_EUR_dat <- pro_dat[pro_dat$eid %in% covar_pro_dat$eid, ]
pro_EUR_mat <- as.matrix(pro_EUR_dat[, -1])
miss_rate_indiv <- aaply(pro_EUR_mat, 1, function(a) sum(is.na(a))/ncol(pro_EUR_mat))
cat("Exclusion: ", sum(miss_rate_indiv > 0.3), "individuals are with missing data of >30% protein.\n")


# Fix eid and covariables
## Build datasets
eid_pro <- pro_EUR_dat$eid[miss_rate_indiv <= 0.3]
covar_pro <- data.frame(eid = covar_pro_dat$eid, 
                        Site = as.factor(covar_pro_dat$`p54_i0`), 
                        Age = covar_pro_dat$`p21022`, 
                        Sex = ifelse(covar_pro_dat$`p22001` == "Female", 0, 1), 
                        BMI = covar_pro_dat$`p21001_i0`, 
                        Array = ifelse(grepl("UKBiLEVEAX", covar_pro_dat$`p22000`), 0, 1), 
                        covar_pro_dat[, paste0("p22009_a", 1:18)])
covar_pro2 <- covar_pro[covar_pro$eid %in% eid_pro, ]
pro_EUR_mat <- pro_EUR_mat[miss_rate_indiv <= 0.3, ]
saveRDS(covar_pro2, paste0(PROJ_PATH, "/02_data/cov/cov_All.rds"))

## Output
write.table(cbind(eid_pro, eid_pro), 
            file = paste0(PROJ_PATH, "/02_data/sample_id/all_eid"), 
            col.names = F, row.names = F, quote = F, sep = "\t")

# Process protein data
## Scale data
miss_rate_pro <- aaply(pro_EUR_mat, 2, function(a) sum(is.na(a))/nrow(pro_EUR_mat))
cat("Exclusion: ", sum(miss_rate_pro > 0.2), "protiens are with missing data of >20% individuals.\n")
pro_fin_mat <- pro_EUR_mat[, miss_rate_pro <= 0.2]
write.table(colnames(pro_fin_mat), file = paste0(PATH, "/protein_label.txt"), 
            row.names = F, col.names = F, quote = F)
pro_fin_imp <- alply(pro_fin_mat, 2, function(a){
  a[is.na(a)] <- mean(a, na.rm = T)
  a_s <- scale(a)
  return(a_s)
}) %>% do.call("cbind", .)
saveRDS(pro_fin_imp, file = paste0(PATH, "/protein_clean_imp.rds"))

## Obtain residual for GWAS in all sex
pro_resid_all <- alply(pro_fin_imp, 2, function(pp){
  
  pro_df <- data.frame(pp, covar_pro2[, -1])
  pro_r <- lm(pp ~ ., data = pro_df)$residuals 
  return(pro_r)
}) %>% do.call("cbind", .)
saveRDS(pro_resid_all, paste0(PATH, "/02_data/GWAS_fam/pro_resid_all.rds"))

## Obtain residual for GWAS in male
pro_fin_imp_male <- pro_fin_imp[covar_pro2$Sex == 1, ]
covar_pro2_male <- covar_pro2[covar_pro2$Sex == 1, ]
write.table(cbind(covar_pro2_male$eid, covar_pro2_male$eid), 
            file = paste0(PROJ_PATH, "/02_data/sample_id/male_eid"), 
            col.names = F, row.names = F, quote = F, sep = "\t")
covar_pro2_male$Sex <- NULL
pro_resid_male <- alply(pro_fin_imp_male, 2, function(pp){
  
  pro_df <- data.frame(pp, covar_pro2_male[, -1])
  pro_r <- lm(pp ~ ., data = pro_df)$residuals 
  return(pro_r)
}) %>% do.call("cbind", .)
saveRDS(pro_resid_male, paste0(PROJ_PATH, "/02_data/GWAS_fam/pro_resid_male.rds"))
saveRDS(covar_pro2_male, paste0(PROJ_PATH, "/02_data/cov/cov_Male.rds"))


## Obtain residual for GWAS in female
pro_fin_imp_female <- pro_fin_imp[covar_pro2$Sex == 0, ] 
covar_pro2_female <- covar_pro2[covar_pro2$Sex == 0, ]
write.table(cbind(covar_pro2_female$eid, covar_pro2_female$eid), 
            file = paste0(PROJ_PATH, "/02_data/sample_id/female_eid"), 
            col.names = F, row.names = F, quote = F, sep = "\t")
covar_pro2_female$Sex <- NULL
pro_resid_female <- alply(pro_fin_imp_female, 2, function(pp){
  
  pro_df <- data.frame(pp, covar_pro2_female[, -1])
  pro_r <- lm(pp ~ ., data = pro_df)$residuals 
  return(pro_r)
}) %>% do.call("cbind", .)
saveRDS(pro_resid_female, paste0(PROJ_PATH, "/02_data/GWAS_fam/pro_resid_Female.rds"))
saveRDS(covar_pro2_female, paste0(PROJ_PATH, "/02_data/cov/cov_Female.rds"))

