# Load packages
library(data.table)
library(dplyr)
library(stringr)
library(ggplot2)
library(optparse)

# Input parameters
args_list <- list(
  make_option("--pro", type="character", default=NULL,
            help="INPUT: protein label", metavar="character")
)
opt_parser = OptionParser(option_list=args_list)
opt = parse_args(opt_parser)
# opt=list(pro="a1bg")
# 
# sex = "all"
# PRO = fread("/gpfs/chencao/ysbioinfor/project/proteohubProject/web_out/protein_label.txt", 
#             header = F)[[1]] %>% as.character()
# 
# plyr::aaply(PRO[1:100], 1, function (pro){
#   
#   print(pro)
#   gwasResults <- fread(paste0(PATH, pro, "/output/summ_", sex, ".assoc.txt.gz"))
#   gwasResults %>% arrange(p_wald) %>% slice(1)
# })

# Set parameters
PATH <- "/gpfs/chencao/ysbioinfor/project/proteohubProject/web_out/pqtl/"
ref_EUR <- fread("/gpfs/chencao/ysbioinfor/Datasets/reference/SNP_BP_GRCh37.txt.gz")

for (sex in c("all", "male", "female")){
  
  # Load data
  gwasResults <- fread(paste0(PATH, opt$pro, "/output/summ_", sex, ".assoc.txt.gz"))
  
  # Manhattan plot
  ## calculate the length of each chromosome
  chr_len <- ref_EUR %>% 
    group_by(CHR) %>% 
    summarise(chr_len = max(BP))
  ## set the start x of each chr
  chr_pos <- chr_len  %>% 
    mutate(total = cumsum(as.numeric(chr_len)) - chr_len) %>%
    select(-chr_len) 
  ## calculate the cumulative position of each SNP
  SNP_info <- chr_pos %>%
    left_join(ref_EUR, ., by="CHR") %>%
    arrange(CHR, BP) %>%
    mutate(BPcum = BP + total)
  # set x axis
  X_axis <-  SNP_info %>% 
    group_by(CHR) %>% 
    summarize(center=(max(BPcum) + min(BPcum))/2)
  SNP_info$PVAL <- gwasResults$p_wald[match(SNP_info$SNP, gwasResults$rs)] 
  SNP_info$PVAL[SNP_info$PVAL == 0] <- min(SNP_info$PVAL[SNP_info$PVAL != 0])
  SNP_info$CHR <- as.factor(SNP_info$CHR)
  SNP_info_sub <- SNP_info[!is.na(SNP_info$PVAL), ]
  ## plot
  up_lim <- ceiling(max(-log10(SNP_info_sub$PVAL))/10)*10
  Manh_plt <- ggplot(SNP_info_sub) +
    geom_point(aes(x = BPcum, y = -log10(PVAL), color = CHR),
               alpha = 0.8, size = 0.8)+
    scale_x_continuous(label = X_axis$CHR, 
                       breaks= X_axis$center, 
                       limits = range(SNP_info$BPcum)) +
    scale_y_continuous(limits = c(0, up_lim),
                       expand = c(0, 0)) +
    scale_color_manual(values = rep(c("#DE423D", "#3D5587"), 11)) +
    xlab("Chromosome") + ylab(bquote("-"~log[10]~"(P-value)"))+
    geom_hline(yintercept = -log10(5E-8), 
               color = 'red', linewidth = 0.8, linetype = 2) + 
    theme_bw() +
    theme(legend.position="none",
          axis.text = element_text(size = 10, color = "black"),
          axis.title = element_text(size = 12, face = "bold", color = "black"),
          panel.border = element_blank(),
          axis.line.y = element_line(),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank())
  ggsave(paste0(PATH, opt$pro, "/plot/", sex, "_Manhattan.png"), Manh_plt, 
         height = 6, width = 14, dpi = 320)
  
  # QQ plot
  qq_df <- data.frame(obs = -log10(sort(SNP_info_sub$PVAL, decreasing = FALSE)),
                      exp = -log10(ppoints(length(SNP_info_sub$PVAL))))
  qq_plt <- ggplot(data = qq_df, aes(exp, obs))+
    geom_point(alpha = 0.8, color = "grey60")+
    geom_abline(color = "red", linewidth = 0.8, linetype = 2)+
    xlab(bquote("Expected -"~log[10]~"(P-value)"))+
    ylab(bquote("Observed -"~log[10]~"(P-value)"))+
    theme_bw()+
    theme(axis.title = element_text(size = 12, face = "bold"),
          axis.text = element_text(size = 10, color = "black"),
          panel.grid = element_blank(),
          panel.border = element_rect(fill = NA),
          panel.background = element_blank())
  ggsave(paste0(PATH, opt$pro, "/plot/", sex, "_qqplot.png"), qq_plt, 
         height = 6, width = 6, dpi = 320)
}
