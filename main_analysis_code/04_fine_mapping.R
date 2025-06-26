# Load packages
library(data.table)
library(plyr)
library(dplyr)
library(susieR)
library(Rfast)
library(optparse)

# Input parameters
args_list <- list(
  make_option("--pro_g", type="numeric", default=NULL,
              help="INPUT: protein group", metavar="character")
)
opt_parser = OptionParser(option_list=args_list)
opt = parse_args(opt_parser)

# opt <- list(pro_g = 1)

#################################
## Make block LD matrix for each chroms
#################################
# # Load block data
# library(bigsnpr)
# BLOCK_PATH <- "/gpfs/chencao/ysbioinfor/Datasets/LD_block/genome_block/EUR"
# block_pos <- alply(c(1: 22), 1, function(chrom){
# 
#   fread(paste0(BLOCK_PATH, "/chr", chrom, ".bed"), header = F)
# })
# # Load reference panel
# REF_PATH <- "/gpfs/chencao/ysbioinfor/Datasets/reference/2k_ukb/"
# for (chrom in 1: 3){
# 
#   bim_chr <- fread(paste0(REF_PATH, "/hm3_imp/chr", chrom, ".bim"))
#   obj.bigSNP <- snp_attach(paste0(REF_PATH, "/hm3_imp/chr", chrom, ".rds"))
#   G <- obj.bigSNP$genotypes
#   blk_out <- aaply(seq_len(nrow(block_pos[[chrom]])), 1, function(b){
# 
#     block_s <- block_pos[[chrom]][b]
#     ind_chr_b <- which(bim_chr$`V4` >= block_s$V2 & bim_chr$`V4` < block_s$V3)
#     corr_chr_b <- snp_cor(G, ind.col = ind_chr_b)
#     corr_chr_b_l <- list(rs = bim_chr$V2[ind_chr_b], corr = corr_chr_b)
#     saveRDS(corr_chr_b_l, paste0(REF_PATH, "LD_matrix/chr", chrom, "/blk_", b, ".rds"))
#     return(b)
#   })
# }
# 
# Set parameters
PATH <- "/gpfs/chencao/ysbioinfor/project/proteohubProject/web_out/P-G/"
REF_PATH <- "/gpfs/chencao/ysbioinfor/Datasets/reference/2k_ukb/"
BLOCK_PATH <- "/gpfs/chencao/ysbioinfor/Datasets/LD_block/genome_block/EUR"

# Load data
## Process protein data
PRO <- fread("/gpfs/chencao/ysbioinfor/project/proteohubProject/web_out/protein_label.txt",
            header = F)[[1]] %>% as.character()
start_g <- seq(1, 2919, 3)
end_g <- seq(3, 2919, 3)
pro_sub <- PRO[start_g[opt$pro_g]: end_g[opt$pro_g]]
## Get block number
block_num <- alply(c(1: 22), 1, function(chrom){

  fread(paste0(BLOCK_PATH, "/chr", chrom, ".bed"), header = F) %>% nrow
})

# Fine mapping
# pro = pro_sub[1]
for (pro in pro_sub){

  cat("Start: ", pro, "\n")
  fine_mapping_sex <- alply(c("all", "female", "male"), 1, function (sex){

    summstats <- fread(paste0(PATH, pro, "/output/summ_", sex, ".assoc.txt.gz"))
    res_pip <- alply(c(1: 22), 1, function(chrom){
      
      cat("chrom: ", chrom, "\n")
      summstats_chr <- summstats %>% filter(chr == chrom)
      pip_chr <- alply(seq_len(block_num[[chrom]]), 1, function(b){

        corr_chr_b_l <- readRDS(paste0(REF_PATH, "LD_matrix/chr", chrom, "/blk_", b, ".rds"))
        summstats_chr_b <- summstats_chr %>% filter(rs %in% corr_chr_b_l$rs)
        if (length(corr_chr_b_l$rs) > 1){

          ind_chr_b <- corr_chr_b_l$rs %in% summstats_chr$rs
          R_chr_b <- as.matrix(corr_chr_b_l$corr)[ind_chr_b, ind_chr_b]
          pip_chr_b <- try(susie_rss(bhat = summstats_chr_b$beta, shat = summstats_chr_b$se, 
                                     R = R_chr_b, n = 33325, L = 10)$pip,
                           silent = T)
          if (inherits(pip_chr_b, "try-error")){
            
            summstats_chr_b$pip_susie <- NA
          } else {
            
            summstats_chr_b$pip_susie <- pip_chr_b
          }
        } else {

          summstats_chr_b$pip_susie <- NA
        }
        return(summstats_chr_b)
      }) %>% do.call("rbind", .)
      return(pip_chr)
    }) %>% do.call("rbind", .)

    fwrite(res_pip, file = paste0(PATH, pro, "/output/summ_", sex, "2.assoc.txt"), sep = "\t")
    zip_cmd <- paste0("gzip -f ", PATH, pro, "/output/summ_", sex, "2.assoc.txt")
    system(zip_cmd)
    return(sex)
  })
}
