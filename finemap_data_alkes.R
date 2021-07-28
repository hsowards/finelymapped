source("/data/Brown_lab/dap/source.R") #this reads in the library and my functions I've written

args <- commandArgs(trailingOnly=T) #args: $chr $pos $dist $prefix $rsid $altstart $altend

#tongwu's pruned sumstats
#sumstats<- readr::read_tsv(paste0("/data/Brown_lab/dap/sumstats/", args[4], ".GWAS.txt"))

#sumstats for full chr
sumstats_raw <- read.delim(paste0("/data/Brown_lab/PAINTOR/dataset/GWAS_Meta/MEL_all_RSQ0.5_", args[1], "_for_TWAS.txt.gz"), 
                           sep = " ") %>% 
  rename(pos = BP, ref = A1, alt = A2, rsnum = RSID)

ld <- readr::read_delim(paste0("/data/Brown_lab/UKBB_LD_Alkes/", args[4], "_LD.txt.gz"), delim = ",", col_names = F)
#ld <- read.table(paste0("/data/Brown_lab/UKBB_LD_Alkes/", args[4], "_LD.txt.gz"), header = F, sep = ",")

bim <- readr::read_tsv(paste0("/data/Brown_lab/UKBB_LD_Alkes/", args[4], "_LD.gz"))

ld_info <- cbind(bim, ld)

#pulling the arguments and setting the region
pos <- as.numeric(args[2])
dist <- as.numeric(args[3])
start <- (pos - dist)
end <- (pos + dist)

altstart <- as.numeric(args[6])
altend <- as.numeric(args[7])

#Clean sumstats
gwas <- sumstats_raw %>% 
  filter(pos > start & pos < end)

#bim
bim <- readr::read_tsv(paste0("/data/Brown_lab/UKBB_LD_Alkes/", args[4], "_LD.gz")) %>% 
  filter(position > start & position < end)

gwas.bim <- merge(bim, gwas, by.x = "position", by.y = "pos", sort=F)

print("number of mismatches ")
#Find out how many mismatches there are
length(which(gwas.bim$allele1!=gwas.bim$ref))

#Flip z-scores
gwas.bim$Z_adjusted <- gwas.bim$Z_fixed

#Subtract the z-score from 0 for the SNPs which are mismatched with the LD reference
gwas.bim$Z_adjusted[which(gwas.bim$allele1!=gwas.bim$ref)] <- 0-gwas.bim$Z_adjusted[which(gwas.bim$allele1!=gwas.bim$ref)]

sumstats <- gwas.bim %>%
  filter(!is.na(Z_adjusted))

if(altstart > 0) {
  print(paste0("Clean sumstats for alt region:", altstart, "-", altend))
  gwas_alt <- sumstats_raw %>% 
    filter(pos > altstart & pos < altend) %>%
    distinct(pos, .keep_all=T)
  
  head(gwas_alt)
  
  #bim
  bim_alt <- readr::read_tsv(paste0("/data/Brown_lab/UKBB_LD_Alkes/", args[4], "_LD.gz")) %>% 
    filter(position > altstart & position < altend)
  
  head(bim_alt)
  
  gwas.bim_alt <- merge(bim_alt, gwas_alt, by.x = "position", by.y = "pos", sort=F)
  
  print("number of mismatches in alt region")
  #Find out how many mismatches there are
  length(which(gwas.bim_alt$allele1!=gwas.bim_alt$ref))
  
  #Flip z-scores
  gwas.bim_alt$Z_adjusted <- gwas.bim_alt$Z_fixed
  
  #Subtract the z-score from 0 for the SNPs which are mismatched with the LD reference
  gwas.bim_alt$Z_adjusted[which(gwas.bim_alt$allele1!=gwas.bim_alt$ref)] <- 0-gwas.bim_alt$Z_adjusted[which(gwas.bim_alt$allele1!=gwas.bim_alt$ref)]
  
  sumstats_alt <- gwas.bim_alt %>%
    filter(!is.na(Z_adjusted))
}

print("Cleaning the LD matrix")
##this is the cleanLD function from my aource file  adjusted for 1000G
##this code aligns the summary statistics and LD matrix
#CLEAN LD
print("Dataset to make positions for headers")
ld_headers <- ld_info %>%
  unite(id, c('allele1', 'allele2', 'rsid', 'position'), sep = "_") %>% 
  #distinct(position, .keep_all=T) %>%
  select(-chromosome)

#duplicate troubleshooting steps
dup_rsids <- ld_info[which(duplicated(ld_info$rsid)),"rsid"]

dup_pos <- ld_info[which(duplicated(ld_info$position)),"position"]

#sumstats %>% 
  #filter(rsnum %in% dup_rsids | pos %in% dup_pos) %>% 
  #write_tsv(paste0("/data/Brown_lab/dap/troubleshooting/", args[4], "_sumstats_dups.txt"))

#ld <- LDlinkR::LDproxy(args[5], pop = "EUR", token = "73c915b7a4da")

#duplicated <- ld_info %>%
   #filter(rsid %in% dup_rsids | position %in% dup_pos) %>%
   #write_tsv(paste0("/data/Brown_lab/dap/troubleshooting/", args[4], "dups.txt"))

#duplicated %>%
   #select(rsid, position) %>%
   #left_join(., ld, by = c("rsid" = "RS_Number"))%>%
   #write_tsv(paste0("/data/Brown_lab/dap/troubleshooting/", args[4], "_dups_LDlink.txt"))

print("pulling headers into a vector")

headers <- ld_headers$id

print("creating parable dataset")
ld_parable <- ld_headers %>% 
  select(-id) %>% 
  set_names(headers) %>% 
  mutate(id = ld_headers$id) %>% 
  select(id, contains("_"))

print("defining extra SNPs")
#defining extra snps
R_extras <- sumstats %>%
   unite(id, c('ref', 'alt', 'rsnum', 'position'), sep = "_") %>%
   anti_join(ld_headers, ., by = c("id")) #if there are no extras? dont proceed?

ex_snps <- R_extras$id #id of extra snps for columns and rows

ld_filter <- ld_parable %>%
  filter(id %in% ex_snps) #keeping missing snp rows

print("clean LD matrix rows")
#clean LD matrix to region (rows)
ld_clean <- ld_parable %>% 
  anti_join(ld_filter, by = c('id')) %>% #antijoin filters out positions that are in R_filter from R_parable
  select(-all_of(ex_snps)) #selecting out missing snp columns

print("clean LD matrix columns")
#clean LD matrix to region (cols)
region_snps <- ld_clean$id
#region_snps_col <- paste("X", region_snps, sep = "")
R_clean <- ld_clean %>% 
  select(all_of(region_snps))

print("setting cleaned matrix as matrix")
R <- as.matrix(R_clean)

ids <- ld_clean$id

sumstats_filtered  <- sumstats %>% #dataset with snps in LD matrix
  unite(id, c('ref', 'alt', 'rsnum', 'position'), sep = "_") %>% 
  filter(id %in% ids) %>%
  distinct(id, .keep_all = T) #removing one duplicated position

if (altstart > 0) {
  #CLEAN LD
  print("Dataset to make positions for headers (alt region)")
  ld_headers <- ld_info %>%
    unite(id, c('allele1', 'allele2', 'rsid', 'position'), sep = "_") %>% 
    #distinct(position, .keep_all=T) %>%
    select(-chromosome)
  
  #duplicate troubleshooting steps
  dup_rsids <- ld_info[which(duplicated(ld_info$rsid)),"rsid"]
  
  dup_pos <- ld_info[which(duplicated(ld_info$position)),"position"]
  
  print("pulling headers into a vector (alt region)")
  
  headers <- ld_headers$id
  
  print("creating parable dataset (alt region)")
  ld_parable <- ld_headers %>% 
    select(-id) %>% 
    set_names(headers) %>% 
    mutate(id = ld_headers$id) %>% 
    select(id, contains("_"))
  
  print("defining extra SNPs (alt region)")
  #defining extra snps
  R_extras <- sumstats_alt %>%
    unite(id, c('ref', 'alt', 'rsnum', 'position'), sep = "_") %>%
    anti_join(ld_headers, ., by = c("id")) #if there are no extras? dont proceed?
  
  ex_snps <- R_extras$id #id of extra snps for columns and rows
  
  ld_filter <- ld_parable %>%
    filter(id %in% ex_snps) #keeping missing snp rows
  
  print("clean LD matrix rows")
  #clean LD matrix to region (rows)
  ld_clean <- ld_parable %>% 
    anti_join(ld_filter, by = c('id')) %>% #antijoin filters out positions that are in R_filter from R_parable
    select(-all_of(ex_snps)) #selecting out missing snp columns
  
  print("clean LD matrix columns")
  #clean LD matrix to region (cols)
  region_snps <- ld_clean$id
  #region_snps_col <- paste("X", region_snps, sep = "")
  R_clean <- ld_clean %>% 
    select(all_of(region_snps))
  
  print("setting cleaned matrix as matrix")
  R_alt <- as.matrix(R_clean)
  
  ids <- ld_clean$id
  
  sumstats_filtered_alt  <- sumstats_alt %>% #dataset with snps in LD matrix
    unite(id, c('ref', 'alt', 'rsnum', 'position'), sep = "_") %>% 
    filter(id %in% ids) %>%
    distinct(id, .keep_all = T) #removing one duplicated position
}

#Writing for dap-g 
sumstats_filtered %>% 
  select(snp_name_i = id, z_i = Z_adjusted) %>% 
  write_tsv(paste0("/data/Brown_lab/dap/dap_input/sumstats/", args[4], "_sumstats_alkes.txt"))

#Writing for dap-g 
R %>% 
  as.data.frame() %>% 
  write_tsv(paste0("/data/Brown_lab/dap/dap_input/ld/", args[4], "_LD_alkes.txt"), col_names = F)

if (altstart > 0) {
  #Writing for dap-g 
  sumstats_filtered_alt %>% 
    select(snp_name_i = id, z_i = Z_adjusted) %>% 
    write_tsv(paste0("/data/Brown_lab/dap/dap_input/sumstats/", args[4], "_altregion_sumstats_alkes.txt"))
  
  #Writing for dap-g 
  R_alt %>% 
    as.data.frame() %>% 
    write_tsv(paste0("/data/Brown_lab/dap/dap_input/ld/", args[4], "_altregion_LD_alkes.txt"), col_names = F)
  
}
