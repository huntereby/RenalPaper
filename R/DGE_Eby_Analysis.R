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

#Example data structure
data("edesign.abiotic")



SampDes <- read.csv("SampleDesign.csv")
SampDes <- SampDes[,-1]

install.packages("dplyr")
library(dplyr)
# Assuming SampDes is your data frame and time_point is the column you want to modify
SampDes <- SampDes %>%
  mutate(time_point = recode(time_point,
                             "preop" = "0",
                             "3wk" = "1",
                             "3mo" = "2",
                             "7mo" = "3"))




rownames(ss.edsign) <- paste("Array", c(1:12), sep = "")

data(data.abiotic)


ss.GENE <- function(n, r, var11 = 0.01, var12 = 0.02, var13 = 0.02,
                    var14 = 0.02, a1 = 0, a2 = 0, a3 = 0, a4 = 0) {
  tc.dat <- NULL
  for (i in 1:n) {
    gene <- c(rnorm(r, a1, var11), rnorm(r, a1, var12),
              rnorm(r, a3, var13), rnorm(r, a4, var14))
    tc.dat <- rbind(tc.dat, gene)
  }
  tc.dat }
flat <-ss.GENE(n = 85, r = 3) # flat
induc <- ss.GENE(n = 5, r = 3, a1 = 0, a2 = 0.2, a3 = 0.6, a4 = 1) # induction
sat <- ss.GENE(n = 5, r = 3, a1 = 0, a2 = 1, a3 = 1.1, a4 = 1.2) # saturation
ord <- ss.GENE(n = 5, r = 3, a1 = -0.8, a2 = -1, a3 = -1.3, a4 =-0.9) # intercept
ss.DATA <- rbind(flat, induc,sat,ord)
rownames(ss.DATA) <- paste("feature", c(1:100), sep = "")
colnames(ss.DATA) <- paste("Array", c(1:12), sep = "")

# This function does not eexsist AND is in there documentation
ss.example <- maSigPro(ss.DATA, ss.edsign, vars="each")
