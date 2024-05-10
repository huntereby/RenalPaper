# This is a file that will be used to create a DGE using time series analysis

# Load count data
urld <- "https://www.ncbi.nlm.nih.gov/geo/download/?format=file&type=rnaseq_counts"
path <- paste(urld, "acc=GSE86884", "file=GSE86884_raw_counts_GRCh38.p13_NCBI.tsv.gz", sep="&");
tbl <- as.matrix(data.table::fread(path, header=T, colClasses="integer"))
                 

# Normalize count data
library(edgeR)
#Set Genes as row names
row.names(tbl)=tbl[,1]
tbl<- tbl[,-1]
#Normalize with cpm from edgeR
tbl_norm <- cpm(tbl)

# Use limma to create a DGE
library(limma)

print(head(tbl_norm))

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

# Create a DGE
f <- factor(SampDes$time_point)
design <- model.matrix(~0+f)
colnames(design) <- levels(f)
fit <- lmFit(tbl_norm, design)

# Fit the DGE
fit <- eBayes(fit)
top <- topTable(fit)
print(top)
