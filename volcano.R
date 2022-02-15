#Mapping differential gene volcano plot
library("labeling")
library('ggplot2')
library('ggpubr')
library("ggrepel")
library("ggthemes")
res <- read.table("F:/major/All_results.csv",sep = ",", header = T)
res$logP <- -log10(res$padj)
head(res)
res$Group = "not-significant"
res$Group[which((res$padj < 0.05) & (res$log2FoldChange > 2))] = "up-regulated"
res$Group[which((res$padj < 0.05) & (res$log2FoldChange < 2))] = "down-regulated"
ggscatter(res, x = "log2FoldChange", y = "logP", color = "Group", palette = c("#2f5688", "#BBBBBB", "#990000"), size = 2, label = res$X, font.label = 8, repel = T, xlab = "log2FoldChange", ylab = "-log10(P.Value)") + theme_base() + geom_hline(yintercept = -log10(0.05), linetype= "dashed") + geom_vline(xintercept = c(-1,1), linetype = "dashed")
