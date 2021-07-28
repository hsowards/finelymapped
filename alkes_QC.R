source("/data/Brown_lab/dap/source.R") #this reads in the library and my functions I've written

args <- commandArgs(trailingOnly=T)

sumstats<- readr::read_tsv(paste0("/data/Brown_lab/dap/sumstats/", args[4], ".GWAS.txt"))

ld <- readr::read_delim(paste0("/data/Brown_lab/UKBB_LD_Alkes/", args[4], "_LD.txt.gz"), delim = ",", col_names = F)
#ld <- read.table(paste0("/data/Brown_lab/UKBB_LD_Alkes/", args[4], "_LD.txt.gz"), header = F, sep = ",")

bim <- readr::read_tsv(paste0("/data/Brown_lab/UKBB_LD_Alkes/", args[4], "_LD.gz"))

ld_info <- cbind(bim, ld)

#pulling the arguments and setting the region
pos <- as.numeric(args[2])
dist <- as.numeric(args[3])
start <- (pos - dist)
end <- (pos + dist)

#Clean sumstats
gwas <- sumstats %>% 
  filter(pos > start & pos < end)

#bim
bim <- readr::read_tsv(paste0("/data/Brown_lab/UKBB_LD_Alkes/", args[4], "_LD.gz")) %>% 
  filter(position > start & position < end)

gwas.bim <- merge(bim, gwas, by.x = "position", by.y = "pos", sort=F)

print("number of mismatches ")

#Find out how many mismatches there are
length(which(gwas.bim$allele1!=gwas.bim$ref))

gwas.bim_mismatch <- gwas.bim[which(gwas.bim$allele1!=gwas.bim$ref),]

write_tsv(gwas.bim_mismatch, paste0("/data/Brown_lab/dap/troubleshooting/", args[4], "_mismatches.txt"))

gwas.bim_mismatch %>% 
  filter(allele1 == "A" | allele1 == "T") %>%
  filter(ref == "C" | ref == "G") %>%
  write_tsv(., paste0("/data/Brown_lab/dap/troubleshooting/", args[4], "_strandflip.txt"))
