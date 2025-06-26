# Load packages
library(plyr)
library(dplyr)
library(data.table)

PROJ_PATH <- "/gpfs/chencao/ysbioinfor/project/proteohubProject/"
DATA_PATH <- paste0(PROJ_PATH, "02_data/")
field_id <- read.table(paste0(DATA_PATH, "model1_label.txt"))[, 1]


sex = "female"
for (sex in c("all", "female", "male")){
  
  sig_emd <- alply(field_id, 1, function(fac_id){
    # fac_id=field_id[113]
    if(grepl("PHE", fac_id)) {
      
      OUT_PATH <- paste0(PROJ_PATH, "/web_out/E-D/", fac_id, "/")
    }  
    if(grepl("PHM", fac_id)) {
      
      OUT_PATH <- paste0(PROJ_PATH, "/web_out/M-D/", fac_id, "/")
    }
    file_s <- paste0(OUT_PATH, "/effect_size_", sex, ".tsv")
    if(file.exists(file_s)){
      
      sig_emd_s <- fread(paste0(OUT_PATH, "/effect_size_", sex, ".tsv"))
      if(nrow(sig_emd_s) != 0){
        
        sig_emd_s <- sig_emd_s[[1]] %>% cbind(fac_id, .)
      }else {
        
        sig_emd_s <- rep(NA, 2)
      }
    } else {
      
      sig_emd_s <- rep(NA, 2)
    }
    return(sig_emd_s)
  })  %>% do.call("rbind", .)  %>% na.omit()
  print(nrow(sig_emd))
  fwrite(sig_emd, file = paste0(DATA_PATH, "sig_total/emd_", sex, "2.txt"),
         sep = "\t", col.names = F)
  sig_emd1 <- read.table(paste0(DATA_PATH, "sig_total/emd_", sex, ".txt"))
  print(nrow(sig_emd1))
  sig_emd_t <- paste0(sig_emd[, 1], ",", sig_emd[, 2])
  sig_emd1_t<- paste0(sig_emd1[, 1], ",", sig_emd1[, 2])
  
  xx=setdiff(sig_emd_t, sig_emd1_t) %>% strsplit(., ",") %>% do.call("rbind", .)
  fwrite(xx, file = paste0(DATA_PATH, "sig_total/im_emd_", sex, ".txt"),
         sep = "\t", col.names = F)

  
}
