#绘制差异基因火山图
library("labeling")
library('ggplot2')
library('ggpubr')
library("ggrepel")
library("ggthemes")
res <- read.table("F:/major/GSE111168/All_results.csv",sep = ",", header = T)
res$logP <- -log10(res$padj)
head(res)
res$Group = "not-significant"
res$Group[which((res$padj < 0.05) & (res$log2FoldChange > 2))] = "up-regulated"
res$Group[which((res$padj < 0.05) & (res$log2FoldChange < 2))] = "down-regulated"
ggscatter(res, x = "log2FoldChange", y = "logP", color = "Group", palette = c("#2f5688", "#BBBBBB", "#990000"), size = 2, label = res$X, font.label = 8, repel = T, xlab = "log2FoldChange", ylab = "-log10(P.Value)") + theme_base() + geom_hline(yintercept = -log10(0.05), linetype= "dashed") + geom_vline(xintercept = c(-1,1), linetype = "dashed")


#GENE_NAME转换为GENE_ID
library(org.Hs.eg.db)
output.gene_id <- read.csv("gene_name.txt")
DEG.gene_symbol = as.character(output.gene_id$Gene)
DEG.entrez_id = mapIds(x = org.Hs.eg.db,
                       keys = DEG.gene_symbol,
                       keytype = "SYMBOL",
                       column = "ENTREZID")
DEG.entrez_id = na.omit(DEG.entrez_id)


#功能富集分析
library(clusterProfiler)
library(topGO)
library(Rgraphviz)
library(pathview)
library(org.Hs.eg.db)
library(tidyverse)
library("clusterProfiler.dplyr")
library(DOSE)

erich.go.BP = enrichGO(gene = DEG.entrez_id,
                       OrgDb = org.Hs.eg.db,
                       keyType = "ENTREZID",
                       ont = "BP",
                       pvalueCutoff = 0.01,
                       qvalueCutoff = 0.05)
BP <- mutate(erich.go.BP, RichFactor = parse_ratio(GeneRatio) / parse_ratio(BgRatio))
dotplot(BP, x = "RichFactor", showCategory = 20, orderBy = "RichFactor") + ggtitle("Gene Ontology Biological Progress Terms Enrichment")

erich.go.CC = enrichGO(gene = DEG.entrez_id,
                        OrgDb = org.Hs.eg.db,
                        keyType = "ENTREZID",
                        ont = "CC",
                        pvalueCutoff = 0.01,
                        qvalueCutoff = 0.05)
CC <- mutate(erich.go.CC, RichFactor = parse_ratio(GeneRatio) / parse_ratio(BgRatio))
dotplot(CC, x = "RichFactor", showCategory = 20, orderBy = "RichFactor") + ggtitle("Gene Ontology Cell Component Terms Enrichment")

erich.go.MF = enrichGO(gene = DEG.entrez_id,
                        OrgDb = org.Hs.eg.db,
                        keyType = "ENTREZID",
                        ont = "MF",
                        pvalueCutoff = 0.01,
                        qvalueCutoff = 0.05)
MF <- mutate(erich.go.MF, RichFactor = parse_ratio(GeneRatio) / parse_ratio(BgRatio))
dotplot(MF, x = "RichFactor", showCategory = 20, orderBy = "RichFactor") + ggtitle("Gene Ontology Molecular Function Terms Enrichment")

enrich.KEGG = enrichKEGG(gene = DEG.entrez_id, organism = 'hsa', pvalueCutoff = 0.01, qvalueCutoff = 0.05)
KEGG <- mutate(enrich.KEGG, RichFactor = parse_ratio(GeneRatio) / parse_ratio(BgRatio))
dotplot(KEGG, x = "RichFactor", showCategory = 20, orderBy = "RichFactor") + ggtitle("KEGG Enrichment Analysis")

enrich.DO = enrichDO(gene = DEG.entrez_id,ont = "DO",pvalueCutoff = 0.01,qvalueCutoff = 0.05)
DO <- mutate(enrich.DO, RichFactor = parse_ratio(GeneRatio) / parse_ratio(BgRatio))
dotplot(DO, x = "RichFactor", showCategory = 20, orderBy = "RichFactor") + ggtitle("Disease Ontology (DO) Enrichment Analysis")
