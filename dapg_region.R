source("/data/Brown_lab/dap/source.R")

args <- commandArgs(trailingOnly=T)

pos <- as.numeric(args[2])

sumstats_raw <- read.delim(paste0("/data/Brown_lab/PAINTOR/dataset/GWAS_Meta/MEL_all_RSQ0.5_", args[1], "_for_TWAS.txt.gz"), 
                           sep = " ") %>%
   filter(BP >= (pos - 500000) & BP <= (pos + 500000)) %>% 
   arrange(P.R.)

#set region analyzed by dap-G
dist <- as.numeric(args[3])
start <- (pos - dist)
end <- (pos + dist)

#calculating LLR
#for low p-value (used this one)
sumstats_raw$chisq <- qchisq(p = sumstats_raw$P.R., df = 1, lower.tail = FALSE)
#calculate LLR
sumstats_raw$likelihood <- exp(sumstats_raw$chisq/2)
sumstats_raw$LLR <- as.numeric(sumstats_raw$likelihood[1]/sumstats_raw$likelihood)

#write_tsv(sumstats_raw, paste0(args[2], "LLR_calc.txt"))

#print("Region within 1:1000  LLR")
#sumstats_raw %>%
   #filter(LLR <= 10000) %>%
   #filter(P <= 10E-7) %>%
   #summarise(start_LLR = min(BP), end_LLR = max(BP), median_LLR = median(BP))

print("Region within 1:1000 conditional LLR")
sumstats_raw %>%
   filter(LLR <= 1000) %>%
   summarise(start_LLR = min(BP), end_LLR = max(BP), median_LLR = median(BP))

print("Region withing +/- 250 kbp")
sumstats_raw %>%
   filter(BP >= start & BP <= end) %>%
   summarise(start_250 = min(BP), end_250 = max(BP))

#print(sumstats_raw)
