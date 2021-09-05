

export PATH=/software/biosoft/software/Rsoft/R3.6/bin:$PATH
R
library(DESeq2)

counts <- read.csv("ciGdal-scr.gene_count_matrix.csv",head=1)

row.names(counts) <- counts$gene_id

counts <- counts[,c(4,5,2,3)]

# condition=factor(c("scr","scr","ci","ci"))

colData <- data.frame(rep=factor(c(1,2,1,2));condition=factor(c("scr","scr","ci","ci")))
rownames(colData) <- colnames(counts)
counts_cds <- DESeqDataSetFromMatrix(countData=counts,colData=colData,design=~condition)

keep <- rowSums(counts(counts_cds) >= 10) >= 2
counts_cds_keep <- counts_cds[keep,]

counts_cds_keep <- estimateSizeFactors(counts_cds_keep)

counts_cds_keep <- estimateDispersions(counts_cds_keep)

res <- nbinomWaldTest(counts_cds_keep)
res <- results(res)
write.table(res,file="ciGdal1_DESeq2Res.txt",sep="\t",row.names=T,col.names=T,quote=FALSE)

