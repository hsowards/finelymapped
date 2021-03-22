source("/data/Brown_lab/dap/source.R") #this reads in the library and my functions I've written

args <- commandArgs(trailingOnly=T)

#1000G data processing for dap-g
sumstats_tongwu <- readr::read_tsv(paste0("/data/Brown_lab/dap/sumstats/", args[4], ".GWAS.txt"))

ld_tongwu <- readr::read_tsv(paste0("/data/Brown_lab/dap/LD_1000G/", args[4], ".LD.gz"), col_names = F)

#pulling the arguments and setting the region
pos <- as.numeric(args[2])
dist <- as.numeric(args[3])
start <- (pos - dist)
end <- (pos + dist)

#Clean sumstats
gwas <- sumstats_tongwu %>% 
  filter(pos > start & pos < end)

#bim
bim <- readr::read_tsv(paste0("/data/Brown_lab/dap/LD_1000G/", args[4], ".LD.gz"), col_names = F) %>% 
  select(X1, X2, X3, X4, X5) %>% 
  filter(X2 > start & X2 < end)

gwas.bim <- merge(bim, gwas, by.x = "X2", by.y = "pos", sort=F)

#Find out how many mismatches there are
length(which(gwas.bim$A1!=gwas.bim$ref))

#Flip z-scores
gwas.bim$Z_adjusted <- gwas.bim$Z_fixed

# Subtract the z-score from 0 for the SNPs which are mismatched with the LD reference
gwas.bim$Z_adjusted[which(gwas.bim$A1!=gwas.bim$ref)] <- 0-gwas.bim$Z_adjusted[which(gwas.bim$A1!=gwas.bim$ref)]

sumstats <- gwas.bim %>%
  filter(!is.na(zscore))

##this is the cleanLD function from my aource file  adjjusted for 1000G
##this code aligns the summary statistics and LD matrix
#CLEAN LD
ld_headers <- ld_tongwu %>% 
  select(-X1, pos = X2, -X3, -X4, -X5)

ld_headers$pos <- paste("X", ld_headers$pos, sep="")

headers <- ld_headers$pos

ld_parable <- ld_headers %>% 
  select(-pos) %>% 
  set_names(headers) %>% 
  mutate(pos = ld_tongwu$X2) %>% 
  select(pos, contains("X"))

#defining extra snps
R_extras <- anti_join(ld_tongwu, sumstats, by = c("X2")) #if there are no extras? dont proceed?
ex_snps <- R_extras$X2 #pos of extra snps
ex_snps_cols <- paste("X", ex_snps, sep="") #column of extra snps

ld_filter <- ld_parable %>%
  filter(pos %in% ex_snps) #keeping missing snp rows

#clean LD matrix to region (rows)
ld_clean <- ld_parable %>% 
  anti_join(ld_filter, by = c('pos')) %>% #antijoin filters out positions that are in R_filter from R_parable
  select(-all_of(ex_snps_cols)) %>% #selecting out missing snp columns
  filter(pos > start & pos < end)

#clean LD matrix to region (cols)
region_snps <- ld_clean$pos
region_snps_col <- paste("X", region_snps, sep = "")
R_clean <- ld_clean %>% 
  select(region_snps_col)

R <- as.matrix(R_clean)

bps <- ld_clean$pos

sumstats_filtered  <- sumstats %>% #dataset with snps in LD matrix
  filter(X2 %in% bps) %>% 
  distinct(X2, .keep_all = T) #removing one duplicated position

#Writing for dap-g 
sumstats_filtered %>% 
  select(snp_name_i = rsnum, z_i = zscore) %>% 
  write_tsv(paste0("/data/Brown_lab/dap/dap_input/sumstats/", args[4], "_sumstats_1000G.txt"))

#Writing for dap-g 
R %>% 
  as.data.frame() %>% 
  write_tsv(paste0("/data/Brown_lab/dap/dap_input/ld/", args[4], "_LD_1000G.txt"), col_names = F)