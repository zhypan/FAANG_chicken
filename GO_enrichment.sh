#!/usr/bin/Rscript
args = commandArgs(trailingOnly=TRUE)
setwd(args[1]) 
library("BiocGenerics",lib.loc = "/home/zhypan/R/library")
library("Biobase",lib.loc = "/home/zhypan/R/library")
library("S4Vectors",lib.loc = "/home/zhypan/R/library")
library("IRanges",lib.loc = "/home/zhypan/R/library")
library("AnnotationDbi",lib.loc = "/home/zhypan/R/library")
library("org.Hs.eg.db",lib.loc = "/home/zhypan/R/library")
library("clusterProfiler", lib.loc = "/home/zhypan/R/library")
library("pathview", lib.loc = "/home/zhypan/R/library")


gs <- read.csv(args[2],header = F)
gs <- gs$V1
go <- enrichGO(gene          = gs,
               keyType       = 'ENSEMBL',
               OrgDb         = org.Hs.eg.db,
               ont           = "BP",
               pAdjustMethod = "fdr",
               pvalueCutoff  = 0.05,
               qvalueCutoff  = 0.05,
               readable      = TRUE)

goo <- as.data.frame(go)
write.table(goo, file=paste("GO_", args[2], sep=""), sep="\t", quote=F, row.names = FALSE)


pdf(file=paste("GO_", args[2], ".pdf", sep=""), width=5, height=5)
cnetplot(go, categorySize = "pvalue")
dev.off() 


pdf(file=paste("GO_", args[2], "heatplot.pdf", sep=""), width=5, height=5)
dotplot(go, showCategory=10)
dev.off() 

#!/bin/bash
module load R/4.1.0
cd /group/zhougrp/zhangyuan/dataanalysis/ChromHMM/OUTPUTSAMPLE_chicken_all_16_modify/state_variability/AAGs/Target_gene_pair
ls output_*_human.txt| while read id;
do
Rscript /group/zhougrp/zhangyuan/script/100_GO_enrichemnt_1.R /group/zhougrp/zhangyuan/dataanalysis/ChromHMM/OUTPUTSAMPLE_chicken_all_16_modify/state_variability/AAGs/Target_gene_pair ${id}
done