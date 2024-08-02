setwd("DEseq2/")

library("R.utils")

# Installing packages
if(!requireNamespace("BiocManager", quietly = TRUE))
{
    install.packages("BiocManager")
}
BiocManager::install("Rsubread")
BiocManager::install("DESeq2")
BiocManager::install("EnhancedVolcano")


library("Rsubread")
library(xlsx)
library(tidyverse)
library(data.table)
library(EnhancedVolcano)
library("DESeq2")

directory <- "./subreads_counts/"

#read table genes for the gene names
genes <- read.table("FB_genes_symbols.tsv", header = TRUE)

#create sample table for DEseq
sampleFiles <- grep(".txt", list.files(directory),value=TRUE)
sampleFiles

sampleCondition <- sub("*[0-9]*.txt","",sampleFiles)
sampleCondition <- sub("*[0-9]*_merged","",sampleCondition)
sampleCondition
sampleTable <- data.frame(sampleName = sampleFiles,
                          fileName = sampleFiles,
                          condition = sampleCondition)
sampleTable

####### running the Deseq2 object dds
ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                       directory = directory,
                                       design= ~ condition)
ddsHTSeq

#pre-filtering of genes with reads > 10
keep <- which(rowSums(counts(ddsHTSeq)) >= 10)
ddsHTSeq <- ddsHTSeq[keep,]

##pca clustering with vst transformed counts
vsd <- vst(ddsHTSeq, blind=FALSE)
pdf("PCA.pdf")
plotPCA(vsd)
dev.off()

#### differential expression analysis
dds <- DESeq(ddsHTSeq)
res <- results(dds)
res

res_DvsWT <- results(dds, contrast = c("condition", "D", "WT"))

ddsHTSeq_normalized <- estimateSizeFactors(ddsHTSeq)
ddsHTSeq_normalized_counts <- counts(ddsHTSeq_normalized, normalized=TRUE)
ddsHTSeq_normalized_counts <- as.data.frame(ddsHTSeq_normalized_counts)   
ddsHTSeq_normalized_counts<- add_rownames(ddsHTSeq_normalized_counts, var = "GENEID")

final_D_vs_WT <- data.frame(
  GENEID = row.names(res_DvsWT),
  log2BaseMean = log2(res_DvsWT$baseMean),
  log2Ratio = res_DvsWT$log2FoldChange,
  Stderr_log2Ratio = res_DvsWT$lfcSE,
  pvalue = res_DvsWT$pvalue,
  padjust = res_DvsWT$padj
  # assays(dds)$counts
)

final_D_vs_WT <- merge(final_D_vs_WT, ddsHTSeq_normalized_counts, by.x = "GENEID", by.y = "GENEID", all.x = TRUE)
final_D_vs_WT <- merge(final_D_vs_WT, genes, by.x = "GENEID", by.y = "FBgn", all.x = TRUE)
write.table(final_D_vs_WT, file="final_D_vs_WT_normCounts.tsv", sep="\t", row.names=FALSE, quote=FALSE)

final_D_vs_WT_padj0.01_FC2 <- final_D_vs_WT[which(final_D_vs_WT$padjust < 0.01 & abs(final_D_vs_WT$log2Ratio) > 1), ]
write.table(final_D_vs_WT_padj0.01_FC2, file="final_D_vs_WT_padj0.01_FC2.tsv", sep="\t", row.names=FALSE, quote=FALSE)



pdf(file = "./Volcano_plot.pdf")   # The directory you want to save the file in
EnhancedVolcano(res_DvsWT,
                lab = rownames(res),
                x = 'log2FoldChange',
                y = 'padj',
                title = 'D vs WT',
                selectLab = c('FBgn0005677','FBgn0028523','FBgn0032586','FBgn0020416','FBgn0020415','FBgn0020414','FBgn0000307','FBgn0028506','FBgn0025678','FBgn0265578'),
                labCol = 'black',
                labFace = 'bold',
                drawConnectors = TRUE,
                lengthConnectors = unit(0.02, 'npc'),
                widthConnectors = 1.0,
                colConnectors = 'black',
                boxedLabels = TRUE,
                pCutoff = 0.01,
                FCcutoff = 1,
                pointSize = 2.0,
                labSize = 4.0)
dev.off() 

pdf(file = "./Volcano_plot_2.pdf")   # The directory you want to save the file in
EnhancedVolcano(res_DvsWT,
                lab = rownames(res),
                x = 'log2FoldChange',
                y = 'padj',
                title = 'D vs WT',
                selectLab = c('FBgn0005677','FBgn0028523','FBgn0032586','FBgn0020416','FBgn0020415','FBgn0020414','FBgn0000307','FBgn0028506','FBgn0025678','FBgn0265578'),
                labCol = 'black',
                labFace = 'bold',
                drawConnectors = TRUE,
                # lengthConnectors = unit(0.02, 'npc'),
                maxoverlapsConnectors = Inf,
                # widthConnectors = 1.0,
                # colConnectors = 'black',
                boxedLabels = TRUE,
                pCutoff = 0.01,
                FCcutoff = 1,
                pointSize = 2.0,
                col=c('black', 'black', 'black', 'red3'),
                labSize = 4.0)
dev.off() 
