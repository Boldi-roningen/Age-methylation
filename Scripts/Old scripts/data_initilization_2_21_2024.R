library(readxl)

#peregrin directory path
directory_path <- "/scratch/s5678684/Data_analysis/Raw_data/methyl_male_no_w/"

# Initialize an empty list to store the data frames
cov_data_list <- list()

# Loop through files, read data, add Sample column, and set column names
for (i in 1:100) {
  file_path <- paste0(directory_path, "S", i, ".deduplicated.bismark.cov.gz")
  cov_data <- read.delim(file_path)
  cov_data$Sample <- as.character(i)
  
  # Set column names
  colnames(cov_data) <- c("Chromosome", "StartPosition", "EndPosition", "MethylationPercentage", "CountMethylated", "CountNonMethylated", "Sample")
  
  cov_data_list[[paste0("Cov", i)]] <- cov_data
}

# Combine the data frames into a single data frame
CovN <- do.call(rbind, cov_data_list)

# Replace Chromosome names
CovN[CovN == "NC_063213.1"] <- "1"
CovN[CovN == "NC_063214.1"] <- "2"
CovN[CovN == "NC_063215.1"] <- "3"
CovN[CovN == "NC_063216.1"] <- "4"
CovN[CovN == "NC_063217.1"] <- "5"
CovN[CovN == "NC_063218.1"] <- "6"
CovN[CovN == "NC_063219.1"] <- "7"
CovN[CovN == "NC_063220.1"] <- "8"
CovN[CovN == "NC_063221.1"] <- "9"
CovN[CovN == "NC_063222.1"] <- "10"
CovN[CovN == "NC_063223.1"] <- "11"
CovN[CovN == "NC_063224.1"] <- "12"
CovN[CovN == "NC_063225.1"] <- "13"
CovN[CovN == "NC_063226.1"] <- "14"
CovN[CovN == "NC_063227.1"] <- "15"
CovN[CovN == "NC_063228.1"] <- "16"
CovN[CovN == "NC_063229.1"] <- "17"
CovN[CovN == "NC_063230.1"] <- "18"
CovN[CovN == "NC_063231.1"] <- "19"
CovN[CovN == "NC_063232.1"] <- "20"
CovN[CovN == "NC_063233.1"] <- "21"
CovN[CovN == "NC_063234.1"] <- "22"
CovN[CovN == "NC_063235.1"] <- "23"
CovN[CovN == "NC_063236.1"] <- "24"
CovN[CovN == "NC_063237.1"] <- "25"
CovN[CovN == "NC_063238.1"] <- "26"
CovN[CovN == "NC_063239.1"] <- "27"
CovN[CovN == "NC_063240.1"] <- "28"
CovN[CovN == "NC_063241.1"] <- "29"
CovN[CovN == "NC_063242.1"] <- "30"
CovN[CovN == "NC_063243.1"] <- "31"
CovN[CovN == "NC_063244.1"] <- "32"
CovN[CovN == "NC_063245.1"] <- "33"
CovN[CovN == "NC_063246.1"] <- "34"
CovN[CovN == "NC_063247.1"] <- "35"
CovN[CovN == "NC_063248.1"] <- "36"
CovN[CovN == "NC_063249.1"] <- "37"
CovN[CovN == "NC_063250.1"] <- "38"
CovN[CovN == "NC_063251.1"] <- "39"
CovN[CovN == "NC_063252.1"] <- "40"
CovN[CovN == "NC_063253.1"] <- "41"
CovN[CovN == "NC_063255.1"] <- "Z"
CovN[CovN == "NC_063254.1"] <- "W"
CovN[CovN == "NC_026783.1"] <- "MT"

# Make coverage column
CovN$Coverage<- CovN$CountMethylated+CovN$CountNonMethylated

saveRDS(CovN, "/scratch/s5678684/Results/Cov_ZF_2024_MalesnoW.rds")

#Reduced table with random extracts
set.seed(123)
number <- nrow(CovN)
random_numbers <- sample(number, 10000, replace = FALSE)
CovN_extract <- CovN[random_numbers, ]
saveRDS(CovN_extract, "/scratch/s5678684/Results/Cov_ZF_2024_MalesnoW_Subset.rds")
