#!/bin/bash

########################
## Fit the mediation model for Genotype to Disease (G-P-D)
## Model1: protein ~ genotype + age + sex + BMI + site + 18 PCs
## Model2: logit(disease) ~ sig. protein + age + sex + BMI + site + 18 PCs

# Load packages
library(plyr)
library(dplyr)
library(data.table)
library(medflex)
library(optparse)

# Input parameters
args_list = list(
  make_option("--sex", type = "character", default=NULL,
              help="INPUT: sex label", metavar = "character"), 
  make_option("--snp_num", type = "integer", default=NULL,
              help="INPUT: snp number", metavar = "character")
)
opt_parser = OptionParser(option_list=args_list)
opt = parse_args(opt_parser)

# opt = list(sex = "female", snp_num = 1)

# Set parameter
ALPHA <- 0.05
ALPHA_OUT <- 1
PROJ_PATH <- "/gpfs/chencao/ysbioinfor/project/proteohubProject/"
DATA_PATH <- paste0(PROJ_PATH, "02_data/")
pro_header <- c("Direct Effect Size", "Direct Effect S.E", "Direct Effect P",
                "Direct Effect Padj","Direct Effect LLCI ", "Direct Effect ULCI",
                "Indirect Effect Size", "Indirect Effect S.E", "Indirect Effect P",
                "Indirect Effect Padj","Indirect Effect LLCI ", "Indirect Effect ULCI")
# Load data
## Select sig. snp in genotype data
sig_snp <- read.table(paste0(DATA_PATH, "sig_genotype/sig_snp_", opt$sex, ".txt"),
                      header = F)[opt$snp_num, 1]
geno_mat <- fread(paste0(DATA_PATH, "sig_genotype/sig_geno_", opt$sex, ".raw"))
sig_idx <- which(strsplit(colnames(geno_mat), "_") %>% laply(., function(gn) gn[1]) == sig_snp)
sig_geno_mat <- geno_mat[[sig_idx]] %>% as.vector()
## Select protein data
sig_pqtl <- paste0(DATA_PATH, "sig_genotype/sig_pqtl_", opt$sex, ".txt") %>% 
  fread(, header = F) %>% 
  filter(V1 == sig_snp) %>% 
  as.data.frame()
pro_label <- paste0(PROJ_PATH, "web_out/protein_label.txt") %>% fread(header = F)
pro_dat <- paste0(PROJ_PATH, "/02_data/protein_clean_imp.rds") %>% readRDS
pro_dat <- pro_dat[, pro_label[[1]] %in% sig_pqtl[, 2], drop = F]
if (opt$sex != "all"){
  
  id_dat_all <- fread(paste0(DATA_PATH, "/sample_id/all_eid"))$V1
  id_dat_sex <- fread(paste0(DATA_PATH, "/sample_id/", opt$sex, "_eid"))$V1
  pro_dat <- pro_dat[id_dat_all %in% id_dat_sex, , drop = F]
}
cat(sig_snp, "associates to", nrow(sig_pqtl), "protein(s).\n")
## Select disease data
sig_dis <- paste0(DATA_PATH, "sig_genotype/sig_trait_", opt$sex, ".txt") %>% 
  fread(, header = F) %>% 
  filter(V1 == sig_snp) %>% 
  as.data.frame()
dis_dat <- paste0(DATA_PATH, "/out_pheno/disease_trait_", opt$sex, "_dat.rds") %>% readRDS()
dis_dat <- dis_dat[, colnames(dis_dat) %in% sig_dis[, 2], drop = F]
cat(sig_snp, "associates to", nrow(sig_dis), "diseases(s).\n")
## Process covariables
cov_dat <- readRDS(paste0(DATA_PATH, "cov/cov_", opt$sex, ".rds"))
cov_dat$eid <- NULL
dummy <- model.matrix(~Site, data = cov_dat)
cov_dat <- cbind.data.frame(cov_dat, dummy[, -1])
cov_dat$Site <- NULL
colnames(cov_dat)[grep("SiteStockport", colnames(cov_dat))] <- "SiteStockport"
cov_dat[, grep("p22009", colnames(cov_dat))] <- NULL

# Fit mediation model for each disease
for (dd in 1: ncol(dis_dat)){
  
  ## Fit mediation model for ddth disease
  dis_dat_sub = dis_dat[, dd]
  OUT_FILE <- paste0(PROJ_PATH, "/web_out/G-P-D/", sig_snp, "/", 
                     colnames(dis_dat)[dd], "_", opt$sex,  ".tsv")
  pro_res <- matrix(NA, nrow = ncol(pro_dat), ncol = length(pro_header)-2)
  for (pro_s in 1: ncol(pro_dat)){
    
    df_tmp <- data.frame(outcome = dis_dat_sub,
                         pro = pro_dat[, pro_s],
                         ms = sig_geno_mat,
                         cov_dat) %>% na.omit
    med_f <- as.formula(paste0("outcome ~ ms + pro + ", 
                               paste0(colnames(cov_dat), collapse = " + ")))
    expData <- neImpute(med_f, family = "binomial", data = df_tmp)
    tot_f <- as.formula(paste0("outcome ~ ms0 + ms1 + ", 
                               paste0(colnames(cov_dat), collapse = " + ")))
    med_fit <- try(neModel(tot_f, expData = expData, se = "robust", family = "binomial"), 
                   silent = T)
    if (inherits(med_fit, "try-error") == F)
      pro_res[pro_s, ] <- c(coef(summary(med_fit))[2, -3], confint(med_fit)["ms0", ],
                            coef(summary(med_fit))[3, -3], confint(med_fit)["ms1", ])
  }
  ## Process output
  if(nrow(na.omit(pro_res)) == 0){
    
    cat("Mediation model is failed in all proteins.\n")
    pro_out_i <- c("protein", pro_header) %>% matrix(., nrow = 1)
    write.table(pro_out_i, row.names = F, sep = "\t", quote = F, col.names = F,
                file = OUT_FILE)
  } else {
    
    dir_eff_padj <- p.adjust(pro_res[, 3], method = "BH")
    ind_eff_padj <- p.adjust(pro_res[, 8], method = "BH")
    pro_out <- cbind.data.frame(pro_res[, 1:3], dir_eff_padj, 
                                pro_res[, 4:8], ind_eff_padj, 
                                pro_res[, 9:10])
    colnames(pro_out) <- pro_header
    pro_out_i <- pro_out[pro_out[, 10] <= ALPHA, , drop = F]
    sig_pro_res <- sig_pqtl[pro_out[, 10] <= ALPHA, 2]
    cat("Indirect effect of", length(sig_pro_res), "proteins are significant.\n")
    # pro_out
    if(nrow(pro_out_i) != 0){

      if(nrow(pro_out_i) == 1){
        
        pro_out_ii <- as.data.frame(c(protein = sig_pro_res, pro_out_i))
        write.table(pro_out_ii,  sep = "\t", row.names = F, quote = F, file = OUT_FILE)
      } else{
        
        pro_out_ii <- cbind.data.frame(protein = sig_pro_res, pro_out_i)
        fwrite(pro_out_ii,  sep = "\t", file = OUT_FILE)
      }
      
    } else {
      
      cat("No significant proteins in mediation model.\n")
      pro_out_ii <- c("protein", pro_header) %>% matrix(., nrow = 1)
      write.table(pro_out_ii, row.names = F, sep = "\t", quote = F, col.names = F,
                  file = OUT_FILE)
    }
  }
}
