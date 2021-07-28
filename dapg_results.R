source("/data/Brown_lab/dap/source.R")

args <- commandArgs(trailingOnly=T)

results <- readr::read_tsv(paste0("/data/Brown_lab/dap/results/", args[4], "/", args[6]), col_names = F) %>%
   separate(X2, into = c("ref", "alt", "rsnum", "bp"), sep = "_")

#read in sumstats
sumstats <- readr::read_tsv(paste0("/data/Brown_lab/dap/sumstats/", args[4], ".GWAS.txt")) %>% 
  arrange(pvalue)

#calculating LLR
#for low p-value (used this one)
sumstats$chisq <- qchisq(p = sumstats$pvalue, df = 1, lower.tail = FALSE)
#calculate LLR
sumstats$likelihood <- exp(sumstats$chisq/2)
sumstats$LLR <- as.numeric(sumstats$likelihood[1]/sumstats$likelihood)

#pulling the arguments and setting the region
position <- as.numeric(args[2])
dist <- as.numeric(args[3])
start <- (position - dist)
end <- (position + dist)

#sumstats cleaning
gwas <- sumstats %>%
   filter(pos > start & pos < end)

merged <- left_join(gwas, results, by = c("rsnum")) %>% 
  select(rsnum, dapg_cs = X5, PIP = X3, chr, pos, se, zscore, pvalue, or, likelihood, LLR) %>% 
  arrange(pvalue) %>%
  filter(dapg_cs >= 1)

readr::write_tsv(merged, paste0(args[4],"_", args[5], "_alkes_summary_wLLR.txt"))

