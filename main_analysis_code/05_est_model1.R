#!/bin/bash

########################
## Estimate the association between protein and exposure/measurement
## Model_all: protein ~ exposure/measurement + age + sex + BMI + site + 18 PCs
## Model_spec: protein ~ exposure/measurement + age  + BMI + site + 18 PCs

# Load packages
library(plyr)
library(dplyr)
library(data.table)
library(ggplot2)
library(ggrepel)
library(optparse)

# Input parameters
args_list = list(
  make_option("--sex", type="character", default=NULL,
              help="INPUT: sex label", metavar="character"), 
  make_option("--trait_id", type="character", default=NULL,
              help="INPUT: trait id", metavar="character")
)
opt_parser = OptionParser(option_list=args_list)
opt = parse_args(opt_parser)

# opt = list (sex = "all", trait_id = "PHE04005")

# Set parameter
ALPHA <- 0.05
DATA_PATH <- "/gpfs/chencao/ysbioinfor/project/proteohubProject/02_data"
OUT_PATH <- "/gpfs/chencao/ysbioinfor/project/proteohubProject/web_out/"
pro_label <- fread(paste0(OUT_PATH, "protein_label.txt"), header = F)$V1
# pro_label <- pro_label[1:50]

if(grepl("PHE", opt$trait_id)) {
  
  OUT_PATH <- paste0(OUT_PATH, "E-P/")
  DATA_FILE <- paste0(DATA_PATH, "/out_pheno/exposure_factors_", opt$sex, "_dat.rds")
  field_id <- fread(paste0(DATA_PATH, "/fieldID_exposure_factor.txt"))
}  
if(grepl("PHM", opt$trait_id)) {

  OUT_PATH <- paste0(OUT_PATH, "M-P/")
  DATA_FILE <- paste0(DATA_PATH, "/out_pheno/measurement_trait_", opt$sex, "_dat.rds")
  field_id <- fread(paste0(DATA_PATH, "/fieldID_measurement_trait.txt"))
}

# Load data
fac_dat <- readRDS(DATA_FILE)
fac_dat_use <- fac_dat[, colnames(fac_dat) == opt$trait_id]

id_dat_all <- fread(paste0(DATA_PATH, "/sample_id/all_eid"))$V1
cov_dat <- readRDS(paste0(DATA_PATH, "/cov/cov_", opt$sex, ".rds"))
cov_dat$eid <- NULL
# strata_dat <- readRDS(paste0(DATA_PATH, "/cov/strata_all.rds"))
dummy <- model.matrix(~Site, data = cov_dat)
cov_dat <- cbind.data.frame(cov_dat, dummy)
cov_dat$Site <- NULL

pro_dat <- readRDS(paste0(DATA_PATH, "/protein_clean_imp.rds"))
if (opt$sex != "all"){

  id_dat_sex <- fread(paste0(DATA_PATH, "/sample_id/", opt$sex, "_eid"))$V1
  pro_dat <- pro_dat[id_dat_all %in% id_dat_sex, ]
  # strata_dat <- strata_dat[id_dat_all %in% id_dat_sex, ]
}

# Summarize result
pro_summ <- alply(pro_dat, 2, function(pp){
# pro_summ <- alply(pro_dat[, 1:50], 2, function(pp){

  # pp=pro_dat[, 1]
  df_tmp <- data.frame(pro = pp, fac = fac_dat_use, cov_dat)
  fit_mod <- lm(pro ~ ., data = df_tmp)
  # summary(fit_mod)
  
  ci <- confint(fit_mod, names(coef(fit_mod))[2], method = "Wald") %>% round(5)
  or_ci <- paste0(ci[1], ", ", ci[2])
  res_tmp <- c(coef(fit_mod)[2] , or_ci, coef(summary(fit_mod))[2, c(3, 4)])
  return(res_tmp)
}) %>% do.call("rbind", .)

pro_result <- data.frame(protein = pro_label,
                         beta = as.numeric(pro_summ[, 1]),
                         CI =  pro_summ[, 2],
                         Z = pro_summ[, 3],
                         P = as.numeric(pro_summ[, 4]))
pro_result$P_Bon <- p.adjust(pro_result$P, method = "bonferroni")
pro_result$P_BH <- p.adjust(pro_result$P, method = "BH")
pro_result$diffexpressed <- "NO"
pro_result$diffexpressed[pro_result$P_BH < ALPHA & pro_result$beta > 0] <- "UP"
pro_result$diffexpressed[pro_result$P_BH < ALPHA & pro_result$beta < 0] <- "DOWN"
pro_show <- pro_result[order(pro_result$P_BH), ] %>% head(5)

# # volcano plot
# trait_label <- field_id$`Reported Trait`[field_id$PHID == opt$trait_id]
# volcano_plt <- ggplot(data = pro_result,
#                       aes(x = beta, y = -log10(P_BH), col = diffexpressed)) +
#   geom_vline(xintercept = 0, col = "gray", linetype = 'dashed') +
#   geom_hline(yintercept = -log10(ALPHA), col = "gray", linetype = 'dashed') +
#   geom_point(size = 1.5) +
#   geom_label_repel(data = pro_show,
#                    aes(label = protein), show.legend = F,
#                    force = 2, nudge_y = 1) +
#   labs(color = "Association", y = expression("-log"[10]*"FDR")) +
#   theme_classic(base_size = 15) +
#   theme(axis.title.y = element_text(face = "bold", margin = margin(0,20,0,0),
#                                     size = rel(1.1), color = 'black'),
#         axis.title.x = element_text(hjust = 0.5, face = "bold", margin = margin(20,0,0,0),
#                                     size = rel(1.1), color = 'black'),
#         plot.title = element_text(hjust = 0.5, size = 10),
#         legend.position="bottom") +
#   ggtitle(paste0("Protein ~ ", trait_label, " + Covs"))
# ## color dots
# if(all(pro_result$diffexpressed == "NO")){
# 
#   volcano_plt <- volcano_plt +
#     scale_color_manual(values = "grey", labels = "Not Sig.")
# } else {
# 
#   if ("DOWN" %in% pro_result$diffexpressed){
# 
#     volcano_plt <- volcano_plt +
#       scale_color_manual(values = c("#00AFBB", "grey", "#BB0C00"),
#                          labels = c("Negative", "Not Sig.", "Positive"))
#   } else {
# 
#     volcano_plt <- volcano_plt +
#       scale_color_manual(values = c("grey", "#BB0C00"),
#                          labels = c("Not Sig.", "Positive"))
#   }
# }
# ggsave(filename = paste0(OUT_PATH, opt$trait_id, "/volcano_plot_", opt$sex, ".png"),
#        volcano_plt, width = 6, height = 5, units = "in", dpi = 320)

# Output data
result_out <- pro_result %>%  filter(P < ALPHA) %>% select(!diffexpressed)
result_out$beta <- result_out$beta %>% round(5)
fwrite(result_out, sep = "\t",
       file = paste0(OUT_PATH, opt$trait_id, "/effect_size_", opt$sex, ".tsv"))
# system(paste0("gzip -f ", OUT_PATH, opt$trait_id, "/effect_size_", opt$sex, ".tsv"))
