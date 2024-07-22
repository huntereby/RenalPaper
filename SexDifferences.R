#Identify the male and female GSM in the GEO dataset
library(GEOquery)
library(limma)
library(Biobase)

# Define the GEO accession number
geo_acc <- "GSE50084"  # Replace with your GEO dataset ID

gset <- getGEO("GSE50084", GSEMatrix =TRUE, AnnotGPL=TRUE)
if (length(gset) > 1) idx <- grep("GPL6244", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]
data <- gset[[1]]
exprs_data <- exprs(gset)


#table from Geo to R 
tT <- read.csv("ebytT.csv")

# Define Y-linked genes
x_genes <- c("FOXL2", "ESR1", "ESR2", "XIST","AMH")
y_genes <- c("RPS4Y1", "KDM5D", "UTY", "ZFY", "USP9Y", "DDX3Y", "EIF1AY", "TMSB4Y","SRY")

# Match row names (probe IDs) with gene symbols
id_to_gene <- setNames(tT$Gene.symbol, tT$ID)
new_row_names <- id_to_gene[rownames(exprs_data)]
rownames(exprs_data) <- new_row_names

# Define Y-linked genes (you might need to customize this list based on your dataset)
y_probes <- rownames(exprs_data)[rownames(exprs_data) %in% y_genes]
y_exprs <- exprs_data[y_probes, ]

# Define X-linked genes (you might need to customize this list based on your dataset)
x_probes <- rownames(exprs_data)[rownames(exprs_data) %in% x_genes]
x_exprs <- exprs_data[x_probes, ]

