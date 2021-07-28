library(BiocManager)
library(motifbreakR)
library(BSgenome.Hsapiens.UCSC.hg19)
library(tidyverse)

##install and load relevant datasets for RSIDs
BiocManager::install("SNPlocs.Hsapiens.dbSNP144.GRCh37")
library(SNPlocs.Hsapiens.dbSNP144.GRCh37)

df <- readxl::read_xlsx("AHR_CCVs.xlsx")

rsids <- df$rsID

snps.mb <- snps.from.rsid(rsid = rsids,
                          dbSNP = SNPlocs.Hsapiens.dbSNP144.GRCh37,
                          search.genome = BSgenome.Hsapiens.UCSC.hg19)

results <- motifbreakR(snpList = snps.mb, filterp = TRUE,
                       pwmList = motifbreakR_motif,
                       threshold = 1e-4, #might want to change this
                       method = "ic", #might want to change this as well
                       bkg = c(A=0.25, C=0.25, G=0.25, T=0.25),
                       BPPARAM = BiocParallel::bpparam())

write.csv(results, "AHR_MB_results.csv")

pdf("motifbreakR_rs117132860.pdf") #printing plot (followed with dev.off)
plotMB(results = results, rsid = "rs117132860", effect = "strong")
dev.off()