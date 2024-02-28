#Merges sample data with info about individuals, and creates Cpg Column
library(readr)
library(dplyr)
Covn <- readRDS("/scratch/s5678684/Results/Cov_ZF_2024_MalesnoW.rds")%>%
  mutate(Sample = as.numeric(Sample))

Covn_extract <- tibble(readRDS("./Results/Cov_ZF_2024_MalesnoW_Subset.rds"))%>%
  mutate(Sample = as.numeric(Sample))

info <- readr::read_csv("./Raw_data/methyl_male_no_w/Zf_SampleInfo_July2023.csv")
print("Data loaded")

merge_Covn <- merge(Covn, info, by = "Sample")
merge_Covn_extract <- merge(Covn_extract, info, by = "Sample")
print("Data merged")

merge_Covn$cpg_site <- paste(merge_Covn$Chromosome, merge_Covn$StartPosition, sep = "_")
merge_Covn_extract$cpg_site <- paste(merge_Covn_extract$Chromosome, merge_Covn_extract$StartPosition, sep = "_")
print("Cpg site column created")

saveRDS(merge_Covn, "/scratch/s5678684/Results/Cov_ZF_SampleInfo_Feb2024.rds")
saveRDS(merge_Covn_extract, "./Results/Cov_ZF_SampleInfo_Feb2024_Subset.rds")
print("Data saved, process finished")
