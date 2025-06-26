
# Load packages
library(data.table)
library(dplyr)
library(plyr)
library(ggplot2)
library(ggrepel)

# Function 1: volcano plot
volcano.plt <- function(pqtl_ss, trait_label){
  
  pqtl_ss$diffexpressed <- "NO"
  pqtl_ss$diffexpressed[pqtl_ss$P_BH < ALPHA & pqtl_ss$beta > 0] <- "UP"
  pqtl_ss$diffexpressed[pqtl_ss$P_BH < ALPHA & pqtl_ss$beta < 0] <- "DOWN"
  pro_show <- pqtl_ss[order(pqtl_ss$P_Bon), ] %>% head(5)
  
  volcano_plt <- ggplot(data = pqtl_ss,
                        aes(x = beta, y = -log10(P_BH), col = diffexpressed)) +
    geom_vline(xintercept = 0, col = "gray", linetype = 'dashed') +
    geom_hline(yintercept = -log10(ALPHA), col = "gray", linetype = 'dashed') +
    geom_point(size = 1.5) +
    geom_label_repel(data = pro_show,
                     aes(label = protein), show.legend = F,
                     force = 2, nudge_y = 1) +
    labs(color = "Association", y = expression("-log"[10]*"FDR")) +
    theme_classic(base_size = 15) +
    theme(axis.title.y = element_text(face = "bold", margin = margin(0,20,0,0),
                                      size = rel(1.1), color = 'black'),
          axis.title.x = element_text(hjust = 0.5, face = "bold", margin = margin(20,0,0,0),
                                      size = rel(1.1), color = 'black'),
          plot.title = element_text(hjust = 0.5, size = 10),
          legend.position="bottom") +
    ggtitle(paste0("Protein ~ ", trait_label, " + Covs"))
  ## color dots
  if(all(pqtl_ss$diffexpressed == "NO")){
    
    volcano_plt <- volcano_plt +
      scale_color_manual(values = "grey", labels = "Not Sig.")
  } else {
    
    if ("DOWN" %in% pqtl_ss$diffexpressed){
      
      volcano_plt <- volcano_plt +
        scale_color_manual(values = c("#00AFBB", "grey", "#BB0C00"),
                           labels = c("Negative", "Not Sig.", "Positive"))
    } else {
      
      volcano_plt <- volcano_plt +
        scale_color_manual(values = c("grey", "#BB0C00"),
                           labels = c("Not Sig.", "Positive"))
    }
  }
  return(volcano_plt)
}

# Set parameters
PROJ_PATH <- "/gpfs/chencao/ysbioinfor/project/proteohubProject"
PQTL_PATH <- paste0(PROJ_PATH, "/web_out/pqtl")
GP_PATH <- paste0(PROJ_PATH, "/web_out/G-P")
PRO_LABEL <- read.table(paste0(PROJ_PATH, "/02_data/protein_label.txt"))[, 1]
ALPHA <- 0.05
# sex = "male"

# Load data
snp_list <- paste0(PROJ_PATH, "/02_data/sig_genotype/sig_snp_", sex, ".txt") %>% fread(header = F)
snp_pqtl <- alply(seq(PRO_LABEL), 1, function(pro_n){

  if (pro_n %% 100 == 0) cat("Protein:", pro_n, "\n")
  pro <- PRO_LABEL[pro_n]
  snp_pqtl_s <- fread(paste0(PQTL_PATH, "/", pro, "/output/summ_", sex, "2.assoc.txt.gz")) %>%
    left_join(snp_list, ., by = join_by(V1 == rs), keep = TRUE)
})


# pqtl_out <- ldply(seq(snp_list[[1]]), function(ss){
pqtl_out <- parallel::mclapply(seq(snp_list[[1]]), function(ss){
  
  pqtl_ss <- alply(seq(snp_pqtl), 1, function(pro_n){
    
    pqtl_snp <- snp_pqtl[[pro_n]][ss, c("beta", "se", "p_wald", "pip_susie")]
    pqtl_snp_ci_up <- pqtl_snp$beta + 1.96*pqtl_snp$se
    pqtl_snp_ci_low <- pqtl_snp$beta - 1.96*pqtl_snp$se
    pqtl_snp_out <- data.frame(protein = PRO_LABEL[pro_n], 
                               beta = pqtl_snp$beta, 
                               CI = paste0(round(pqtl_snp_ci_low, 5), "," , round(pqtl_snp_ci_up, 5)),
                               Z = pqtl_snp$beta/pqtl_snp$se, 
                               P = pqtl_snp$p_wald, 
                               PIP = pqtl_snp$pip_susie)
    return(pqtl_snp_out)
  }) %>% do.call("rbind",.)
  pqtl_ss$P_Bon <- p.adjust(pqtl_ss$P, "bonferroni")
  pqtl_ss$P_BH <- p.adjust(pqtl_ss$P, "BH")
  
  ## Output volcano plot
  vol_plt_ss <- volcano.plt(pqtl_ss, snp_list[[1]][ss])
  ggsave(file = paste0(PROJ_PATH, "/web_out/G-P/", snp_list[[1]][ss], "/volcano_plot_", sex, ".png"),
         vol_plt_ss, width = 6, height = 5, units = "in", dpi = 320)
  ## Output Sig. pQTL
  pqtl_ss_out <- pqtl_ss[pqtl_ss$P < ALPHA, ]
  fwrite(pqtl_ss_out, file = paste0(PROJ_PATH, "/web_out/G-P/", snp_list[[1]][ss], "/effect_size_", sex, ".tsv"), 
         sep = "\t", quote = F)
  return(ss)
}, mc.cores = 12, mc.allow.recursive = TRUE)




# Set parameter
sex_label = "All"

for (sex in c("all", "female", "male")){
# for (sex in c("female", "male")){
  
  if (sex == "female") sex_label = c("All", "Female")
  if (sex == "male") sex_label = c("All", "Male")
  
  # Load data
  dis_field_id <- fread(paste0(PROJ_PATH, "/02_data/fliedID_disease_trait.txt")) %>%
    filter(`Sex Specificity` %in% sex_label)
  snp <- paste0(PROJ_PATH, "snp_FDR/common/", sex, 
                "/bonferroni_005_bonferroni_0.05-pip_0.8_trait.csv") %>%
    read.csv(., header = T) %>%
    select(snp_name) %>%
    unique
  
  # Load summary statistics for each disease
  dis_summ_list <- alply(dis_field_id[[1]], 1, function(dis_code){
    
    paste0(PROJ_PATH, "web_out/GWAS-D/", dis_code, "/output/summ_", sex, ".assoc.txt.gz") %>%
      fread(., header = T) %>%
      filter(rs %in% snp[, 1])
  })
  
  # dis_out <- parallel::mclapply(seq(snp[, 1]), function(ss){
  dis_out <- alply(seq(snp[, 1]), 1, function(ss){ 
    # ss=1
    ss_x <- which(dis_summ_list[[1]]$rs == snp[ss, 1])
    dis_ss <- alply(seq(dis_summ_list), 1, function(dis_n){
      
      dis_snp <- dis_summ_list[[dis_n]][ss_x, c("rs", "beta", "se", "p_wald")]
      dis_snp <- dis_summ_list[[dis_n]][ss_x, c("beta", "se", "p_wald")]
      dis_snp_ci_up <- dis_snp$beta + 1.96*dis_snp$se
      dis_snp_ci_low <- dis_snp$beta - 1.96*dis_snp$se
      dis_snp_out <- data.frame(FieldID = dis_field_id[dis_n]$PHID,
                                Trait = dis_field_id[dis_n]$`Reported Trait`, 
                                Category = dis_field_id[dis_n]$`ICD10 category`, 
                                # rs = dis_snp$rs,
                                beta = dis_snp$beta, 
                                CI = paste0(round(dis_snp_ci_low, 5), ", " , round(dis_snp_ci_up, 5)),
                                Z = dis_snp$beta/dis_snp$se, 
                                P = dis_snp$p_wald)
      return(dis_snp_out)
    }) %>% do.call("rbind",.)
    # dis_ss$P_Bon <- p.adjust(dis_ss$P, "bonferroni")
    # dis_ss$P_BH <- p.adjust(dis_ss$P, "BH")
    
    # Figure
    total_result_top <- dis_ss %>%
      mutate(low = as.numeric(tstrsplit(CI, ", ")[[1]])) %>%
      mutate(up = as.numeric(tstrsplit(CI, ", ")[[2]])) %>%
      arrange(P) %>%
      head(n=5)
    total_result_top$Trait <- factor(total_result_top$Trait, levels = total_result_top$Trait)
    plt <- ggplot(total_result_top, aes(x = Trait, y = beta, fill = Trait)) +
      geom_bar(stat = "identity", width = 0.8) +
      geom_errorbar(aes(ymin = low, ymax = up),
                    size = 0.8, width = 0.2) +
      geom_text(aes(label = round(beta, 4)), vjust = -0.5,
                size = 4) +
      scale_fill_manual(values = c("#800000", "#2A9D8F", "#46F0F0", "#AAFFC3", "#6A3569")) +
      xlab("Top Five Diseases") +
      ylab("Effect Size") +
      ggtitle(paste0("Disease ~ ", snp[ss, 1], " + Covs"))+
      theme_bw() +
      theme(axis.text.x = element_text(size = 12, color = "black", angle = 70, hjust = 1),
            axis.title.x = element_text(size = 15, face = "bold", color = "black"),
            axis.text.y = element_text(size = 12, color = "black"),
            axis.title.y = element_text(size = 15, face = "bold", color = "black"),
            panel.grid = element_blank(),
            strip.placement = "outside",
            strip.background = element_blank(),
            strip.text = element_blank(),
            legend.position = "none")
    ggsave(paste0(PROJ_PATH, "/web_out/G-D/", snp[ss, 1], "/top_trait_", sex, ".png"),
           plt,
           width = 6, height = 6, units = "in", dpi = 320)
    
    ## Output Sig. trait
    if(sum(dis_ss$P < ALPHA) < 5){
      
      dis_ss_out <- dis_ss %>% arrange(P) %>% head(5)
    } else {
      cat(ss, "ok\n")
      dis_ss_out <- dis_ss[dis_ss$P < ALPHA, ]
    }
    
    out_file <- paste0(PROJ_PATH, "/web_out/G-D/", snp[ss, 1], "/effect_size_", sex, ".tsv")
    fwrite(dis_ss_out, file = out_file, 
           sep = "\t", quote = F)
    return(out_file)
  })
  # }, mc.cores = 4, mc.allow.recursive = TRUE)
}
