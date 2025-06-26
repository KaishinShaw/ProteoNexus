# Load packages
library(bigreadr)
library(dplyr)
library(plyr)

# Set parameters
PROJ_PATH <- "/public/home/biostat03/project/proteohubProject"
TRAIT <- "/public/home/Datasets/ukb/pheno/03_trait/participant.csv.gz"
NONCANCER_DATE <- "/public/home/Datasets/ukb/pheno/03_trait/noncancer_date.rds"
NONCANCER_ID <- "/public/home/Datasets/ukb/pheno/03_trait/fieldID_noncancer_date.txt"

# Load sample id
all_id <- fread2(paste0(PROJ_PATH, "/02_data/sample_id/all_eid"))[, 1]
male_id <- fread2(paste0(PROJ_PATH, "/02_data/sample_id/male_eid"))[, 1]
female_id <- fread2(paste0(PROJ_PATH, "/02_data/sample_id/female_eid"))[, 1]

# Load disease
## Set  flied id
cancer_self_report <- c(paste0("p20001_i0_a", c(0:5)), paste0("p20001_i1_a", c(0:5)),
                        paste0("p20001_i2_a", c(0:5)), paste0("p20001_i3_a", c(0:5)))
cancer_self_report_age <- c(paste0("p20007_i0_a", c(0:5)), paste0("p20007_i1_a", c(0:5)),
                            paste0("p20007_i2_a", c(0:5)), paste0("p20007_i3_a", c(0:5)))
cancer_icd10 <- paste0("p40006_i", c(0: 21))
cancer_icd10_age <- paste0("p40008_i", c(0: 21))
noncancer_self_report <- c(paste0("p20002_i0_a", c(0:33)), paste0("p20002_i1_a", c(0:33)),
                           paste0("p20002_i2_a", c(0:33)), paste0("p20002_i3_a", c(0:33)))
noncancer_self_report_age <- c(paste0("p20009_i0_a", c(0:33)), paste0("p20009_i1_a", c(0:33)),
                                paste0("p20009_i2_a", c(0:33)), paste0("p20009_i3_a", c(0:33)))
## Load binary data
trait_binary_dat <- fread2(TRAIT, 
                           select = c("eid", 
                                      "p21022",         # Age at recruitment
                                      "p34",            # Year of birth
                                      "p41202",         # Main diagnose
                                      "p41204",         # Second diagnose
                                      cancer_self_report, cancer_self_report_age,
                                      cancer_icd10, cancer_icd10_age, 
                                      noncancer_self_report, noncancer_self_report_age))
trait_binary_dat <- trait_binary_dat[match(all_id, trait_binary_dat$eid), ]
sample_size <- nrow(trait_binary_dat)
recruitment_age <- trait_binary_dat$p21022
noncancer_date_df <- readRDS(NONCANCER_DATE)
noncancer_date_df <- noncancer_date_df[match(all_id, noncancer_date_df$eid), ]
noncancer_id_date <- fread2(NONCANCER_ID, header= F)

# Process field ID data
field_id <- fread2(paste0(PROJ_PATH, "/02_data/fliedID_disease_trait.txt"))
field_id$ICD10 <- gsub(",", "\\|", field_id$ICD10)

# Identify cancer
trait_binary_dat_cancer <- trait_binary_dat[, c(cancer_self_report, cancer_self_report_age, 
                                                cancer_icd10, cancer_icd10_age)]
field_id_cancer <- field_id[c(1:13), ]
trait_cancer_out <- alply(c(1: nrow(field_id_cancer)), 1, function(tt) {
  
  ## define cases and age
  self_dat <- alply(c(1: sample_size), 1, function(ss){
    # ss=235
    cancer_s <- trait_binary_dat_cancer[ss, c(cancer_self_report, cancer_self_report_age)]
    cnd_s <- any(cancer_s[, cancer_self_report] %in% field_id_cancer$`Self-reported Trait`[tt])
    if(cnd_s){
      # print(ss)
      cancer_idx_s <- which(cancer_s[, cancer_self_report] %in% field_id_cancer$`Self-reported Trait`[tt])
      cancer_age_s <- cancer_s[, (cancer_idx_s + length(cancer_self_report))]
      if (length(cancer_idx_s) >= 2)
        cancer_age_s <- cancer_age_s[1, 1]
      return(c(ss, cancer_age_s))
    }
  }) %>% do.call("rbind", .)
  icd10_dat <- alply(c(1: sample_size), 1, function(ss){
  
    cancer_s <- trait_binary_dat_cancer[ss, c(cancer_icd10, cancer_icd10_age)]
    cnd_s <- any(grepl(field_id_cancer$`ICD10`[tt], cancer_s[, cancer_icd10]))
    if(cnd_s){
      
      cancer_idx_s <- which(grepl(field_id_cancer$`ICD10`[tt], cancer_s[, cancer_icd10]))
      cancer_age_s <- cancer_s[, (cancer_idx_s + length(cancer_icd10))]
      if (length(cancer_idx_s) >= 2)
        cancer_age_s <- min(cancer_age_s)
      return(c(ss, cancer_age_s))
    }
  }) %>% do.call("rbind", .)
  ## define cases before recruiting 
  if (!is.null(icd10_dat) | !is.null(self_dat)){
  
    inter_idx <- intersect(self_dat[, 1], icd10_dat[, 1])
    if (length(inter_idx) > 0){
      
      event_dat <- rbind(self_dat[!self_dat[, 1] %in% inter_idx, ], icd10_dat)
    } else {
      
      event_dat <- rbind(self_dat, icd10_dat)
    }
    del_cnd <- ifelse(event_dat[, 2] < recruitment_age[event_dat[, 1]], 1, 0)
    event_dat <- cbind(event_dat, del_cnd)
    ## output 
    out_s <- rep(0, sample_size)
    out_s[event_dat[event_dat[, 3] == 0, 1]] <- 1
    out_s[event_dat[event_dat[, 3] == 1, 1]] <- NA
    cat("For", field_id_cancer$`Self-reported Trait`[tt], ", ",  
        sum(out_s, na.rm = T), "cases are selected.\n")
    return(out_s)
  }
}) %>% do.call("cbind", .) 

# Identify non-cancer
trait_binary_dat_noncancer <- trait_binary_dat[, c("p41202", "p41204", 
                                                    noncancer_self_report, noncancer_self_report_age)]
field_id_noncancer <- field_id[-c(1: 13), ]
trait_noncancer_out <- alply(c(1: nrow(field_id_noncancer)), 1, function(tt) {
  
  print(field_id_noncancer$`Self-reported Trait`[tt])
  icd10_dat <- self_dat <- NULL
  ## define cases and age
  self_trait <- field_id_noncancer$`Self-reported Trait`[tt]
  if(!is.na(self_trait)){
    self_dat <- alply(c(1: sample_size), 1, function(ss){
      
      noncancer_s <- trait_binary_dat_noncancer[ss, c(noncancer_self_report, noncancer_self_report_age)]
      cnd_s <- any(noncancer_s[, noncancer_self_report] %in% self_trait)
      if(cnd_s){
        
        noncancer_idx_s <- which(noncancer_s[, noncancer_self_report] %in% self_trait)
        noncancer_age_s <- noncancer_s[, (noncancer_idx_s + length(noncancer_self_report))]
        if (length(noncancer_idx_s) >= 2)
          noncancer_age_s <- noncancer_age_s[1, 1]
        if (noncancer_age_s != -1)
          return(c(ss, noncancer_age_s))
      }
    }) %>% do.call("rbind", .)
  } 
  noncancer_date_df_s <- strsplit(field_id_noncancer$`ICD10`[tt], "\\|")[[1]] %>% 
    substr(1, 3) %>% 
    aaply(., 1, function(ss) noncancer_id_date[grep(ss, noncancer_id_date[, 2]), 1]) %>%
    noncancer_date_df[, ., drop = F] 
  if (ncol(noncancer_date_df_s) != 0){
    
    idx_main <- strsplit(field_id_noncancer$`ICD10`[tt], "\\|")[[1]] %>% 
      alply(., 1, function(ss) grepl(ss, trait_binary_dat_noncancer[, "p41202"])) %>%
      do.call("cbind", .)
    main_dat <- alply(c(1: sample_size), 1, function(ss) {
      
      if (any(idx_main[ss, ])){
        
        date_s <- noncancer_date_df_s[ss, , drop = F] 
        date_s2 <- date_s[1, idx_main[ss, ]]
        if (length(date_s2) >= 2)
          date_s2 <- date_s2[1, 1]
        if(!grepl("Code has event date", as.character(date_s2))){
          
          age_s <- difftime(as.Date(date_s2), 
                            as.Date(as.character(trait_binary_dat$p34[ss]), "%Y"), 
                            units="days") %>% as.numeric
          return(c(ss, round(age_s/365, 1)))
        }
      }
    }) %>% do.call("rbind", .)
    idx_second <- strsplit(field_id_noncancer$`ICD10`[tt], "\\|")[[1]] %>% 
      alply(., 1, function(ss) grepl(ss, trait_binary_dat_noncancer[, "p41204"])) %>%
      do.call("cbind", .)
    second_dat <- alply(c(1: sample_size), 1, function(ss) {
      
      if (any(idx_second[ss, ])){

        date_s <- noncancer_date_df_s[ss, , drop = F] 
        date_s2 <- date_s[1, idx_second[ss, ]]
        if (length(date_s2) >= 2)
          date_s2 <- date_s2[1, 1]
        if(!grepl("Code has event date", as.character(date_s2))){
          
          age_s <- difftime(as.Date(date_s2), 
                            as.Date(as.character(trait_binary_dat$p34[ss]), "%Y"), 
                            units="days") %>% as.numeric
          return(c(ss, round(age_s/365, 1)))
        }
      }
    }) %>% do.call("rbind", .)
    inter_diag <- intersect(second_dat[, 1], main_dat[, 1])
    icd10_dat <- rbind(main_dat, second_dat[!second_dat[, 1] %in% inter_diag, ])
  }
  ## define cases before recruiting 
  if (!is.null(icd10_dat) | !is.null(self_dat)){
    
    inter_idx <- intersect(self_dat[, 1], icd10_dat[, 1])
    if (length(inter_idx) > 0){
      
      event_dat <- rbind(self_dat[!self_dat[, 1] %in% inter_idx, ], icd10_dat)
    } else {
      
      event_dat <- rbind(self_dat, icd10_dat)
    }
    del_cnd <- ifelse(event_dat[, 2] < recruitment_age[event_dat[, 1]], 1, 0)
    event_dat <- cbind(event_dat, del_cnd)
    ## output 
    out_s <- rep(0, sample_size)
    out_s[event_dat[event_dat[, 3] == 0, 1]] <- 1
    out_s[event_dat[event_dat[, 3] == 1, 1]] <- NA
    cat("For", field_id_noncancer$`Reported Trait`[tt], ", ",  
        sum(out_s, na.rm = T), "cases are selected.\n")
    return(out_s)
  }
}) %>% do.call("cbind", .) 

# Output 
trait_disease_all_dat <- trait_disease_all_tot <- cbind(trait_cancer_out, trait_noncancer_out)
trait_disease_all_dat[, which(field_id$`Sex Specificity` != "All")] <- rep(NA, nrow(trait_disease_all_dat))
trait_disease_female_dat <- trait_disease_all_tot[all_id%in%female_id, ]
trait_disease_female_dat[, which(field_id$`Sex Specificity` == "Male")] <- rep(NA, nrow(trait_disease_female_dat))
trait_disease_male_dat <- trait_disease_all_tot[all_id%in%male_id, ]
trait_disease_male_dat[, which(field_id$`Sex Specificity` == "Female")] <- rep(NA, nrow(trait_disease_male_dat))
colnames(trait_disease_all_dat) <- colnames(trait_disease_female_dat) <- colnames(trait_disease_male_dat) <- field_id$PHID
saveRDS(trait_disease_all_dat, file = paste0(PROJ_PATH, "/02_data/out_pheno/disease_trait_all_dat.rds"))
saveRDS(trait_disease_female_dat, file = paste0(PROJ_PATH, "/02_data/out_pheno/disease_trait_female_dat.rds"))
saveRDS(trait_disease_male_dat, file = paste0(PROJ_PATH, "/02_data/out_pheno/disease_trait_male_dat.rds"))

# write.csv(trait_disease_all_dat, file = paste0(PROJ_PATH, "/02_data/out_pheno/disease_trait_all_dat.csv"), 
#          quote = F, row.names = F)
# write.csv(trait_disease_female_dat, file = paste0(PROJ_PATH, "/02_data/out_pheno/disease_trait_female_dat.csv"), 
#           quote = F, row.names = F)
# write.csv(trait_disease_male_dat, file = paste0(PROJ_PATH, "/02_data/out_pheno/disease_trait_male_dat.csv"), 
#           quote = F, row.names = F)


write.csv(xx, file = paste0(PROJ_PATH, "/02_data/out_pheno/measurement_trait_all.csv"),
         quote = F, row.names = F)

# Summarize disease traits
trait_disease_all_summ <- aaply(1: ncol(trait_disease_all_dat), 1, function (ll){
  
  label <- trait_disease_all_dat[, ll] %>% na.omit()
  if (length(label) == 0){

    return(NA)
  } else {
    # paste0(length(label), " (", sum(label), ":", sum(label == 0), ")")
    return(paste0(length(label), " (", sum(label), ":", sum(label == 0), ")"))
  }
})
trait_disease_female_summ <- aaply(1: ncol(trait_disease_female_dat), 1, function (ll){
  
  label <- trait_disease_female_dat[, ll] %>% na.omit()
  if (length(label) == 0){
    
    return(NA)
  } else {
    
    return(paste0(length(label), " (", sum(label), ":", sum(label == 0), ")"))
  }
})
trait_disease_male_summ <- aaply(1: ncol(trait_disease_male_dat), 1, function (ll){
  
  label <- trait_disease_male_dat[, ll] %>% na.omit()
  if (length(label) == 0){
    
    return(NA)
  } else {
    
    return(paste0(length(label), " (", sum(label), ":", sum(label == 0), ")"))
  }
})
trait_summ <- cbind.data.frame(field_id$PHID, trait_disease_all_summ,
                               trait_disease_female_summ, trait_disease_male_summ)
write.csv(trait_summ, file = paste0(PROJ_PATH, "/02_data/out_pheno/trait_disease_summ.csv"), 
          quote = F, row.names = F)
