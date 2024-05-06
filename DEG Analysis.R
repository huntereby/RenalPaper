#R script for grabbing the dataset

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
bio_pkgs <- c('DESeq2','GOSemSim','GO.db','limma','ComplexHeatmap',
'AnnotationDbi','GenomeInfoDb','GenomicRanges','Biobase','S4Vectors',
'BiocGenerics','MatrixGenerics','IRanges','BiocFileCache',
'SummarizedExperiment','org.Hs.eg.db','org.Mm.eg.db','org.Ce.eg.db',
'tximport','tximportData', 'GEOquery')
BiocManager::install(bio_pkgs)
install.packages("devtools")
devtools::install_github("Ylefol/TimeSeriesAnalysis@master")
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
  # Replicate "transplant_1" 4 times, then "transplant_2" 4 times, and so on until there's 96 rows total
  replicate = rep(paste0("transplant_", rep(1:24, each = 4)), times = 1),
  time_point = rep(c('preop', '3wk', '3mo', '7mo'), 24))

# Place this command in the terminal:
# Rscript -e "rmarkdown::render('rmarkdown_method/TS_analysis.Rmd',output_file='TS_analysis.html')"
