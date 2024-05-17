# This is a file that will be used to create a DGE using time series analysis
# Install packages
install.packages(c("tidyverse", "R.utils"))
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("edgeR", "sva"))
library(tidyverse)
# Load count data
urld <- "https://www.ncbi.nlm.nih.gov/geo/download/?format=file&type=rnaseq_counts"
path <- paste(urld, "acc=GSE86884", "file=GSE86884_raw_counts_GRCh38.p13_NCBI.tsv.gz", sep="&");
tbl <- as.matrix(data.table::fread(path, header=T, colClasses="integer"))

set.seed(123)

# Normalize count data
library(edgeR)
#Set Genes as row names



print(head(tbl_norm))

SampDes <- read.csv("SampleDesign.csv")
SampDes <- SampDes[,-1]

install.packages("dplyr")
library(dplyr)
# Assuming SampDes is your data frame and time_point is the column you want to modify
SampDes <- SampDes %>%
  mutate(time_point = recode(time_point,
                             "0" = "preop",
                             "3wk" = "wk3",
                             "3mo" = "mo3",
                             "7mo" = "mo7"))

# Get factor levels of the time poitns
time_point  <- recode(levels(as.factor(SampDes$time_point)), "0" = "preop", "1" = "wk3", "2" = "mo3", "3" = "mo7")

design <- model.matrix(~0 + time_point, data=SampDes)

comparisons <- limma::makeContrasts(
  preopv3wk = "time_pointmo3 - time_pointpreop",
  preopv3mo = "time_pointmo3 - time_pointpreop",
  preopv7mo = "time_pointmo7 - time_pointpreop",
  o3wkv3mo = "time_pointmo3 - time_pointwk3",
  o3wkv7mo = "time_pointmo7 - time_pointwk3",
  o3mov7mo = "time_pointmo7 - time_pointmo3",
   levels = design)

contrasts <- comparisons
row.names(tbl)=tbl[,1]
tbl<- tbl[,-1]
x <- tbl

y <- DGEList(x, group = SampDes$time_point)

keep <- filterByExpr(y)
y <- y[keep,, keep.lib.sizes=FALSE]
tbl_norm <- normLibSizes(y)

print(SampDes)
library(edgeR)
library(sva)

model <- edgeR::voomLmFit(tbl_norm, 
                          design = design,
                          plot = F, 
                          block = SampDes$sample,
                          sample.weights = T)

contrasts <- comparisons
results <- list() 
for (i in colnames(contrasts)){
  fit <- contrasts.fit(model, contrasts = contrasts[,i])  
  fit <- eBayes(fit, robust = T)
  results[[i]] <- topTable(fit, coef = colnames(contrasts[,i]), n = Inf) %>%
    rownames_to_column("replicate")
}

# Get only those genes which were significant (below 0.05 adj p value)
results <- results %>%
  map(~ filter(.x, P.Value < 0.05))

results %>%
  imap(~ write_csv(.x, paste0("Results/", .y, ".csv")))


