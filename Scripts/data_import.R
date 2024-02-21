#Merges sample data with info about individuals
library(readr)
Covn <- readRDS("/scratch/s5678684/Results/Cov_ZF_2024_MalesnoW.rds")
Covn_extract <- readRDS("/scratch/s5678684/Results/Cov_ZF_2024_MalesnoW_Subset.rds")
info <- read_csv("/scratch/s5678684/Raw_data/Zf_SampleInfo_July2023.csv")
merge_Covn <- merge(Covn, info, by = "Sample")
merge_Covn_extract <- merge(Covn_extract, info, by = "Sample")

saveRDS(merged_CovN, "/scratch/s5678684/Results/Cov_ZF_SampleInfo_Feb2024.rds")
saveRDS(merged_CovN_extract, "/scratch/s5678684/Results/Cov_ZF_SampleInfo_Feb2024_Subset.rds")


