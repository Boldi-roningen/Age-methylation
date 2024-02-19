#Merges sample data with info about individuals
library(readr)
Covn <- readRDS("")
info <- read_csv("")
merge_Covn <- merge(Covn, info, by = "Sample")


saveRDS(merged_CovN, "Cov_ZF_SampleInfo_Feb2024.rds")
#Merge reduced sample data with info about individuals
Covn_extract <- readRDS
merge_Covn_extract <- merge(Covn_extract, info, by = "Sample")
saveRDS(merged_CovN_extract, "Extract_Cov_ZF_SampleInfo_Feb2024.rds")