#R script for grabbing the dataset
library(GEOquery)


# load counts table from GEO
urld <- "https://www.ncbi.nlm.nih.gov/geo/download/?format=file&type=rnaseq_counts"
path <- paste(urld, "acc=GSE86884", "file=GSE86884_raw_counts_GRCh38.p13_NCBI.tsv.gz", sep="&");
tbl <- as.matrix(data.table::fread(path, header=T, colClasses="integer"), rownames=1)

#Looks into package TiSA, and other time series packages to do in this analysis
#DEG using time series. 
BiocManager::install("GEOmeta")
X <- getGEO("GSE86884")
library(tidyverse)

# Assuming 'tbl' is your existing data frame
# Create a vector of column names
column_names <- colnames(tbl)

# Create a data frame with the desired structure
new_df1 <- tibble(
  sample = rep(column_names),  # Repeat each column name 4 times
  group = rep("transplant", times = length(column_names)),  # Fill 'group' column with "transplant"
  replicate = rep(paste0("transplant_", 1:96), each = 4),
  time_point = rep(paste0('preop', '3wk', '3mo', '7mo'), rep=24))
  # Setting up the replicate column
  # Repeat the string "transplant_1" 4 times, then "transplant_2" 4 times, and so on 
  
  # Combine columns of total_col into a single column
  

  
# View the new data frame
