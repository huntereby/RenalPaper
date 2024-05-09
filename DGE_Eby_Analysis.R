#This is Eby attempt at creating a DGE using time series analysis
#We also want to create a trajactory for the Genes at each time point

# load counts table from GEO
urld <- "https://www.ncbi.nlm.nih.gov/geo/download/?format=file&type=rnaseq_counts"
path <- paste(urld, "acc=GSE86884", "file=GSE86884_raw_counts_GRCh38.p13_NCBI.tsv.gz", sep="&");
tbl <- as.matrix(data.table::fread(path, header=T, colClasses="integer"))
                 


#Normalization of RNASeq tbl
library(edgeR)
#Set Genes as row names
row.names(tbl)=tbl[,1]
tbl<- tbl[,-1]
#Normalize with cpm from edgeR
tbl_norm <- cpm(tbl)


#Going to be using maSigPro for time series analysis
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("maSigPro")
library(maSigPro)
maSigProUsersGuide()
data(data.abiotic)
