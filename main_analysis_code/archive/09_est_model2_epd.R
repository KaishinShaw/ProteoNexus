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
  make_option("--exp_id", type = "integer", default=NULL,
              help="INPUT: exposure factor id", metavar = "character"),
  make_option("--out_id", type = "integer", default=NULL,
              help = "INPUT: disease id", metavar = "character")
)
opt_parser = OptionParser(option_list=args_list)
opt = parse_args(opt_parser)

# opt = list(sex = "female", exp_id = 3, out_id = 14)

# Set parameter
ALPHA <- 0.05
PROJ_PATH <- "/gpfs/chencao/ysbioinfor/project/proteohubProject/"
DATA_PATH <- paste0(PROJ_PATH, "02_data/")
EXP_PATH <- paste0(PROJ_PATH, "/web_out/E-P/")
# source("/gpfs/chencao/ysbioinfor/software/processv43/PROCESS v4.3 for R/process.R")

# Load data
## Protein
pro_label <- paste0(PROJ_PATH, "web_out/protein_label.txt") %>% 
  fread(header = F)
exp_label <- fread(paste0(DATA_PATH, "e_label.txt"), header = F)[[1]]
sig_pro <- fread(paste0(EXP_PATH, exp_label[opt$exp_id], "/effect_size_", opt$sex, ".tsv")) %>% 
  filter(P_BH <= ALPHA)

if (nrow(sig_pro) != 0){
  
  pro_dat <- paste0(PROJ_PATH, "/02_data/protein_clean_imp.rds") %>% readRDS
  pro_dat <- pro_dat[, pro_label[[1]] %in% sig_pro[[1]]]
  if (opt$sex != "all"){
    
    id_dat_all <- fread(paste0(DATA_PATH, "/sample_id/all_eid"))$V1
    id_dat_sex <- fread(paste0(DATA_PATH, "/sample_id/", opt$sex, "_eid"))$V1
    pro_dat <- pro_dat[id_dat_all %in% id_dat_sex, ]
  }
  ## Outcome
  trait_label <- fread(paste0(DATA_PATH, "fliedID_disease_trait.txt"))[[1]]
  outcome_trait <- readRDS(paste0(DATA_PATH, "/out_pheno/disease_trait_", opt$sex, "_dat.rds"))
  outcome_trait <- outcome_trait[, colnames(outcome_trait) == trait_label[opt$out_id]]
  if (!all(is.na(outcome_trait))){
    
    ## Exposure
    exp_factor <- readRDS(paste0(DATA_PATH, "out_pheno/exposure_factors_", opt$sex, "_dat.rds"), )
    exp_factor <- exp_factor[, opt$exp_id]
    ## Covaribles
    cov_dat <- readRDS(paste0(DATA_PATH, "cov/cov_", opt$sex, ".rds"))
    cov_dat$eid <- NULL
    dummy <- model.matrix(~Site, data = cov_dat)
    cov_dat <- cbind.data.frame(cov_dat, dummy[, -1])
    cov_dat$Site <- NULL
    colnames(cov_dat)[grep("SiteStockport", colnames(cov_dat))] <- "SiteStockport"
    cov_dat[, grep("p22009", colnames(cov_dat))] <- NULL
    # Fit mediation model for each protein
    pro_res <- matrix(NA, nrow = ncol(pro_dat), ncol = 8)
    # ncol(pro_dat)
    for (pro_s in 1: ncol(pro_dat)){
      
      cat("Protein: ", pro_s, "\n")
      df_tmp <- data.frame(outcome = outcome_trait,
                           pro = pro_dat[, pro_s],
                           exp = exp_factor,
                           cov_dat) %>% na.omit
      # med_res <- try(process(data = df_tmp, y = "outcome", x = "exp", m = "pro", model = 4,
      #                        cov = colnames(df_tmp)[-c(1:3)], seed = 20250525, 
      #                        modelbt = 0, outscreen = 1, boot = NBOOST, save = 3), 
      #                silent = T)
      med_f <- as.formula(paste0("outcome ~ exp + pro + ", 
                                 paste0(colnames(cov_dat), collapse = " + ")))
      expData <- neImpute(med_f, family = "binomial", data = df_tmp)
      tot_f <- as.formula(paste0("outcome ~ exp0 + exp1 + ", 
                                 paste0(colnames(cov_dat), collapse = " + ")))
      med_fit <- neModel(tot_f, expData = expData, se = "robust", family = "binomial")
      if (inherits(med_fit, "try-error") == F)
        # pro_res[pro_s, ] <- c(coef(summary(med_fit))[2, -3], confint(med_fit)["exp0", ], 
        #                       coef(summary(med_fit))[3, -3], confint(med_fit)["exp1", ])
        pro_res[pro_s, ] <- c(coef(summary(neEffdecomp(med_fit)))[3, -3], 
                              coef(summary(med_fit))[3, -3], confint(med_fit)["exp1", ])
    }

    # Output result
    pro_header <- c("Total Effect Size", "Total Effect S.E", "Total Effect P",
                    "Total Effect Padj", 
                    "Indirect Effect Size", "Indirect Effect S.E", "Indirect Effect P",
                    "Indirect Effect Padj","Indirect Effect LLCI ", "Indirect Effect ULCI")
    dir_eff_padj <- p.adjust(pro_res[, 3], method = "BH")
    ind_eff_padj <- p.adjust(pro_res[, 6], method = "BH")
    pro_out <- cbind.data.frame(pro_res[, 1:3], dir_eff_padj, pro_res[, 4:8], ind_eff_padj)
    colnames(pro_out) <- pro_header
    pro_out_i <- pro_out[pro_out[, 10] <= ALPHA, ]
    sig_pro_res <- sig_pro[[1]][pro_out[, 10] <= ALPHA]
    # pro_out
    if(nrow(pro_out_i) != 0){
      
      pro_out_i <- apply(pro_out_i, 2, function (a) round(a, 5)) 
      pro_out_i <- cbind.data.frame(protein = sig_pro_res, pro_out_i)
      fwrite(pro_out_i,  sep = "\t", 
             file = paste0(PROJ_PATH, "/web_out/E-P-D/", exp_label[opt$exp_id], "/", 
                           trait_label[opt$out_id],"_", opt$sex,  ".tsv"))
    } else {
      
      cat("No significant proteins in mediation model.\n")
      pro_out_i <- c("protein", pro_header) %>% matrix(., nrow = 1)
      write.table(pro_out_i, row.names = F, sep = "\t", quote = F, col.names = F,
                  file = paste0(PROJ_PATH, "/web_out/E-P-D/", exp_label[opt$exp_id], "/", 
                                trait_label[opt$out_id], "_", opt$sex, ".tsv"))
    }
  } else{
    
    cat("Sex-specific trait.\n")
  }
} else {
  
  cat("NO significant protein in model1.\n")
}
