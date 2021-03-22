source("/data/Brown_lab/dap/source.R")

args <- commandArgs(trailingOnly=T)

#read in sumstats
sumstats_raw <- read.delim(paste0("/data/Brown_lab/PAINTOR/dataset/GWAS_Meta/MEL_all_RSQ0.5_", args[1], "_for_TWAS.txt.gz"), 
                           sep = " ")

#pulling the arguments and setting the region
pos <- as.numeric(args[2])
dist <- as.numeric(args[3])
start <- (pos - dist)
end <- (pos + dist)

#sumstats cleaning
gwas <- sumstats_raw %>% 
  filter(BP > start & BP < end)

bim <- genio::read_bim(paste0(args[4], ".bim"))

bim.temp <- bim %>% 
  filter(pos > start & pos < end)

gwas.bim <- merge(bim.temp, gwas, by.x = "pos", by.y = "BP", sort=F)

#Find out how many mismatches there are
length(which(gwas.bim$A1!=gwas.bim$ref))

#Flip z-scores
gwas.bim$Z_adjusted <- gwas.bim$Z_fixed

# Subtract the z-score from 0 for the SNPs which are mismatched with the LD reference
gwas.bim$Z_adjusted[which(gwas.bim$A1!=gwas.bim$ref)] <- 0-gwas.bim$Z_adjusted[which(gwas.bim$A1!=gwas.bim$ref)]

#remove zscores that are NA
sumstats <- gwas.bim %>%
  filter(!is.na(Z_adjusted))

#Read in LD reference panel and create matrix
data <- genio::read_plink(args[4])

ld <- PLINKtoLD(X = data$X, N = 5000) #5000 for UKBB

ld_info <- cbind(bim, ld) #bim (rsid, ref, alt, pos, chr) attached to LD for later reference

nrow(sumstats) == nrow(ld_info)

#Writing for dap-g 
sumstats %>% 
  select(snp_name_i = RSID, z_i = Z_adjusted) %>% 
  write_tsv(paste0("/data/Brown_lab/dap/dap_input/sumstats/", args[4], "_sumstats_UKBB.txt"))

R <- cleanLD(ld_info, start, end) #this function aligns the summary stats and LD matrix

#Writing for dap-g 
R %>% 
  as.data.frame() %>% 
  write_tsv(paste0("/data/Brown_lab/dap/dap_input/ld/", args[4], "_LD_UKBB.txt"), col_names = F)