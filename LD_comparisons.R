library(tidyverse)

sumstats <- readr::read_tsv("data/AHRlocus.GWAS.txt") %>% 
  mutate(rsID = rsnum, chr = as.character(chr))

data <- readxl::read_xlsx("results/Candidate Causals for Annotation.xlsx") %>% 
  mutate(chr = 7) %>% # hardcode changing this to numeric sorry
  select(chr, pos = start, rsID, `Conditional Set`)

#DAP results

#raw <- readr::read_tsv("results/results_1000G_250_3.txt", skip_empty_rows = T)
#this sucks, might have to use some command line magic for these files

ukbb <- readxl::read_xlsx("results/ukbb_4_250.xlsx") %>% 
  left_join(., sumstats, by = "rsID") %>% 
  select(chr, pos, rsID, `UKBB DAP-G Set` = cluster)
onekg <- readxl::read_xlsx("results/1000G_4_250.xlsx") %>% 
  left_join(., sumstats, by = "rsID") %>% 
  select(chr, pos, rsID, `1000G DAP-G Set` = cluster)


dap <- full_join(ukbb, onekg, by = c("chr", "pos", "rsID")) %>% 
  mutate(chr = as.double(chr))

#LLR
llr_by_rs117_rs104 <- readxl::read_xlsx("Matthew Conditional Results, 3OM_with_LLR_KB.v.2.xlsx", sheet = "c.7_2_by_rs11713286_rs10487582") %>% 
  filter(LLR <= 100) %>%
  select(rsID = rsID_dbSNP142, LLR_by_rs117_rs104 = LLR, pos = bp, chr = Chr)

llr_by_rs117_rs730 <- readxl::read_xlsx("Matthew Conditional Results, 3OM_with_LLR_KB.v.2.xlsx", sheet = "c.7_3_by_rs117132860_rs73069846") %>% 
  filter(LLR <= 100) %>% 
  select(rsID = rsID_dbSNP142, LLR_by_rs117_rs730 = LLR, pos = bp, chr = Chr) %>% 
  full_join(., llr_by_rs117_rs104, by = c("chr", "pos", "rsID"))

llr_by_all <- readxl::read_xlsx("Matthew Conditional Results, 3OM_with_LLR_KB.v.2.xlsx", 
                                sheet = "c.7_by_all_3_SNPs") %>% 
  filter(LLR <= 100) %>% 
  select(rsID = rsID_dbSNP142, LLR_by_all = LLR, pos = bp, chr = Chr) %>% 
  full_join(., llr_by_rs117_rs730, by = c("chr", "pos", "rsID")) %>% 
  select(rsID, chr, pos, LLR_by_all, LLR_by_rs117_rs730, LLR_by_rs117_rs104)

new <- full_join(dap, llr_by_all, by = c("chr","pos","rsID"))
summary <- full_join(data, new, by = c("chr","pos","rsID"))

write_csv(summary, "results/fine-mapping_comparison.csv")