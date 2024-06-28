#New version
#Going to do comparisons of pre transplant versus the 3 timepoints
#Then we are going to run 3 pod on the 3 time point comparisons

# Version info: R 4.2.2, Biobase 2.58.0, GEOquery 2.66.0, limma 3.54.0
################################################################
#   Differential expression analysis with DESeq2
#This is my test with positron
install.packages("BiocManager")
library(BiocManager)
BiocManager::install("DESeq2")
library(DESeq2)

# load counts table from GEO
urld <- "https://www.ncbi.nlm.nih.gov/geo/download/?format=file&type=rnaseq_counts"
path <- paste(urld, "acc=GSE86884", "file=GSE86884_raw_counts_GRCh38.p13_NCBI.tsv.gz", sep="&");
tbl <- as.matrix(data.table::fread(path, header=T, colClasses="integer"), rownames="GeneID")

# load gene annotations 
apath <- paste(urld, "type=rnaseq_counts", "file=Human.GRCh38.p13.annot.tsv.gz", sep="&")
annot <- data.table::fread(apath, header=T, quote="", stringsAsFactors=F, data.table=F)
rownames(annot) <- annot$GeneID

# sample selection
gsms <- paste0("02301201230101010123012012301230101301201230101231",
               "2013012301301013013012130120230120120201010101")
sml <- strsplit(gsms, split="")[[1]]

# group membership for samples
gs <- factor(sml)
groups <- make.names(c("PreTransplant","1week","3 months","6 months"))
levels(gs) <- groups
sample_info <- data.frame(Group = gs, row.names = colnames(tbl))

# pre-filter low count genes
# keep genes with at least N counts > 10, where N = size of smallest group
keep <- rowSums( tbl >= 10 ) >= min(table(gs))
tbl <- tbl[keep, ]

ds <- DESeqDataSetFromMatrix(countData=tbl, colData=sample_info, design= ~Group)

#1 Week comparison
res_1_week <- results(ds, contrast=c("Group", groups[1], groups[2]), alpha=0.05, pAdjustMethod = "fdr")
res_1_week <- merge(as.data.frame(res_1_week), annot, by=0, sort=F)
res_1_week <- subset(res_1_week, select=c("Symbol","log2FoldChange","pvalue"))
res_1_week
write.csv(res_1_week, "1weekComp.csv")
#3 month comparison
res_3_week <- results(ds, contrast=c("Group", groups[1], groups[3]), alpha=0.05, pAdjustMethod = "fdr")
res_3_week <- merge(as.data.frame(res_3_week), annot, by=0, sort=F)
res_3_week <- subset(res_3_week, select=c("Symbol","log2FoldChange","pvalue"))
res_3_week
write.csv(res_3_week, "3MonthComp.csv")
#6 month comparison
res_6_week <- results(ds, contrast=c("Group", groups[1], groups[4]), alpha=0.05, pAdjustMethod = "fdr")
res_6_week <- merge(as.data.frame(res_6_week), annot, by=0, sort=F)
res_6_week <- subset(res_3_week, select=c("Symbol","log2FoldChange","pvalue"))
res_6_week
write.csv(res_6_week, "6MonthComp.csv")

#Time series analysis
#Pre transplant to 1 week analysis
time1 <- results(ds, contrast=c("Group", groups[1], groups[2]), alpha=0.05, pAdjustMethod = "fdr")
time1 <- merge(as.data.frame(time1), annot, by=0, sort=F)
time1 <- subset(time1, select=c("Symbol","log2FoldChange","pvalue"))
time1
write.csv(time1, "time1.csv")
#1 week to 3 month analysis
time2 <- results(ds, contrast=c("Group", groups[2], groups[3]), alpha=0.05, pAdjustMethod = "fdr")
time2 <- merge(as.data.frame(time2), annot, by=0, sort=F)
time2 <- subset(time2, select=c("Symbol","log2FoldChange","pvalue"))
time2
write.csv(time2, "time2.csv")
#3 months to 6 month analysis
time3 <- results(ds, contrast=c("Group", groups[3], groups[4]), alpha=0.05, pAdjustMethod = "fdr")
time3 <- merge(as.data.frame(time3), annot, by=0, sort=F)
time3 <- subset(time3, select=c("Symbol","log2FoldChange","pvalue"))
time3
write.csv(time3, "time3.csv")





