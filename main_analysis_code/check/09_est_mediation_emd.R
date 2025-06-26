#!/bin/bash

########################
## Fit the mediation model for Exposure factor to Trait (E-P-T)
## Model1: protein ~ exposure + age + sex + BMI + site + 18 PCs
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
  make_option("--pair", type = "integer", default=NULL,
              help="INPUT: exposure factor id", metavar = "character")
)
opt_parser = OptionParser(option_list=args_list)
opt = parse_args(opt_parser)

# opt = list(sex = "male", pair = 131)

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
## Sig. emd
sig_emd <- fread(paste0(PROJ_PATH, "01_code/check/miss_emd_", opt$sex, ".txt"), header = F)
sig_emd_s = sig_emd[opt$pair, ]
if(grepl("PHE", sig_emd_s[[1]])) {
  
  FAC_FILE <- paste0(DATA_PATH, "/out_pheno/exposure_factors_", opt$sex, "_dat.rds")
  IN_PATH <- paste0(PROJ_PATH, "/web_out/E-P/", opt$fac_id, "/")
  OUT_FILE <- paste0(PROJ_PATH, "/web_out/E-P-D/", sig_emd_s[[1]], "/", 
                     sig_emd_s[[2]], "_", opt$sex,  ".tsv")
  # field_id <- fread(paste0(DATA_PATH, "/fieldID_exposure_factor.txt"))
}  
if(grepl("PHM", sig_emd_s[[1]])) {
  
  FAC_FILE <- paste0(DATA_PATH, "/out_pheno/measurement_trait_", opt$sex, "_dat.rds")
  IN_PATH <- paste0(PROJ_PATH, "/web_out/M-P/", opt$fac_id, "/")
  OUT_FILE <- paste0(PROJ_PATH, "/web_out/M-P-D/", sig_emd_s[[1]], "/", 
                     sig_emd_s[[2]], "_", opt$sex,  ".tsv")
  # field_id <- fread(paste0(DATA_PATH, "/fieldID_measurement_trait.txt"))
}

# Load Protein
pro_label <- paste0(PROJ_PATH, "web_out/protein_label.txt") %>% 
  fread(header = F)
sig_pro <- fread(paste0(IN_PATH, sig_emd_s[[1]], "/effect_size_", opt$sex, ".tsv")) %>% 
  filter(P_Bon <= ALPHA)
# sig_pro = sig_pro[1:100,]
cat(sig_emd_s[[1]], "associates to", nrow(sig_pro), "proteins.\n")

if (nrow(sig_pro) != 0){
  
  pro_dat <- paste0(PROJ_PATH, "/02_data/protein_clean_imp.rds") %>% readRDS
  pro_dat <- pro_dat[, pro_label[[1]] %in% sig_pro[[1]], drop = F]
  if (opt$sex != "all"){
    
    id_dat_all <- fread(paste0(DATA_PATH, "/sample_id/all_eid"))$V1
    id_dat_sex <- fread(paste0(DATA_PATH, "/sample_id/", opt$sex, "_eid"))$V1
    pro_dat <- pro_dat[id_dat_all %in% id_dat_sex, , drop = F]
  }
  ## Outcome
  outcome_trait <- readRDS(paste0(DATA_PATH, "/out_pheno/disease_trait_", opt$sex, "_dat.rds"))
  outcome_trait <- outcome_trait[, colnames(outcome_trait) == sig_emd_s[[2]]]

  ## Exposure
  fac_dat <- readRDS(FAC_FILE)
  fac_dat <- fac_dat[, sig_emd_s[[1]]]
  ## Covariables
  cov_dat <- readRDS(paste0(DATA_PATH, "cov/cov_", opt$sex, ".rds"))
  cov_dat$eid <- NULL
  dummy <- model.matrix(~Site, data = cov_dat)
  cov_dat <- cbind.data.frame(cov_dat, dummy[, -1])
  cov_dat$Site <- NULL
  colnames(cov_dat)[grep("SiteStockport", colnames(cov_dat))] <- "SiteStockport"
  cov_dat[, grep("p22009", colnames(cov_dat))] <- NULL
  # Fit mediation model for each protein
  pro_res <- matrix(NA, nrow = ncol(pro_dat), ncol = length(pro_header)-2)
  # ncol(pro_dat)
  for (pro_s in 1: ncol(pro_dat)){
    
    if (pro_s%%20 == 0) cat("Protein: ", pro_s, "\n")
    df_tmp <- data.frame(outcome = outcome_trait,
                         pro = pro_dat[, pro_s],
                         ms = fac_dat %>% as.numeric,
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
    
  # Output result
  if(nrow(na.omit(pro_res)) == 0){
    
    cat("Mediation model is failed in all proteins.\n")
    pro_out_i <- c("protein", pro_header) %>% matrix(., nrow = 1)
    write.table(pro_out_i, row.names = F, sep = "\t", quote = F, col.names = F,
                file = OUT_FILE)
  } else {
    
    dir_eff_padj <- p.adjust(pro_res[, 3], method = "BH")
    ind_eff_padj <- p.adjust(pro_res[, 8], method = "BH")
    if(nrow(pro_res) == 1){
      
      pro_out <- matrix(c(pro_res[1:3], dir_eff_padj, pro_res[4:8], 
                            ind_eff_padj, pro_res[9:10]), nrow=1)
    } else {
      
      pro_out <- cbind.data.frame(pro_res[, 1:3], dir_eff_padj, 
                                  pro_res[, 4:8], ind_eff_padj, 
                                  pro_res[, 9:10])
    }
    colnames(pro_out) <- pro_header
    pro_out_i <- pro_out[pro_out[, 10] <= ALPHA, ]
    sig_pro_res <- sig_pro[[1]][pro_out[, 10] <= ALPHA]
    cat("Indirect effect of", length(sig_pro_res), "proteins are significant.\n")
    # pro_out
    if(nrow(pro_out_i) != 0){
      
      pro_out_i <- apply(pro_out_i, 2, function (a) round(a, 5)) 
      if(length(pro_out_i) == length(pro_header)){
        
        pro_out_ii <- as.data.frame(c(protein = sig_pro_res, pro_out_i)) %>% t
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
} else {
  
  cat("No significant protein in model1.\n")
  pro_out_ii <- c("protein", pro_header) %>% matrix(., nrow = 1)
  write.table(pro_out_ii, row.names = F, sep = "\t", quote = F, col.names = F,
              file = OUT_FILE)
}