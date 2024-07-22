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
pheno_data <- pData(data)

exprs_data <- exprs(gset)

# Download the platform annotation data
platform_id <- annotation(gset)
platform_data <- getGEO(platform_id)
platform_table <- Table(platform_data)

# Define Y-linked genes
y_genes <- c("FOXL2", "ESR1", "ESR2", "XIST","AMH")

# Match row names (probe IDs) with gene symbols
id_to_gene <- setNames(tT$Gene.symbol, tT$ID)
new_row_names <- id_to_gene[rownames(exprs_data)]
rownames(exprs_data) <- new_row_names

# Find probes corresponding to Y-linked genes
y_probes <- rownames(exprs_data)[gene_symbols %in% y_genes]

# Define Y-linked genes (you might need to customize this list based on your dataset)
y_genes <- c("RPS4Y1", "KDM5D", "UTY", "ZFY", "USP9Y", "DDX3Y", "EIF1AY", "TMSB4Y")
y_probes <- rownames(exprs_data)[rownames(exprs_data) %in% y_genes]

y_exprs <- exprs_data[y_probes, ]
