# Title     : Identify DEGs for bio repeat samples (n>4) by wilcoxson rank test
# Created by: stl23
# Created on: 2021/8/5

#library(degeR)
library(tibble)
args <- commandArgs(trailingOnly = TRUE)
groups <- args[1]
input <- args[2]
path <- args[3]
prefix <- args[4]
pvalue_limit <- as.numeric(args[5])
FC_limit <- as.numeric(args[6])  ## log2(FC=1.5)~0.585     log2(FC=2) = 1
control_group <- args[7]     # control
treatment_group <- args[8]   # case


if (!dir.exists(paths = path)) {
      stop("Directory provided does not exist")
    }
#input_file <- file.path(path,'./before.all.wt.mut.combined.changeid.txt')
input_file <- file.path(input)
#group_name <- unlist(strsplit(groups,','))
group_info <- read.table(file = groups, sep="\t", header=T,row.names = 1)
out_gene_tsv <- paste0(path,'/',prefix,'.RNAseq_gene_results.xls')
out_DE_gene <- paste0(path,'/',prefix,'.RNAseq_different_expression_genes_results.xls')
out_up_gene <- paste0(path,'/',prefix,'.RNAseq_UP_genes_results.xls')
out_down_gene <- paste0(path,'/',prefix,'.RNAseq_DOWN_genes_results.xls')

# Load gene data(Geneid file_name1 file_name2...)
genecountData <- read.table(file = input_file, sep="\t", row.names = 1, header=T,check.names = F)
sample_names <- colnames(genecountData)

#conditions <- factor(group_name$group,levels = unique(group_name$group))  ## order limited: c("Control","Treatment") --> Treatment vs. Control
# 按样本名称匹配分组
conditions <- factor(group_info[sample_names, "group"], 
                    levels = c(control_group, treatment_group))  # 明确水平顺序
conditions
# Get target matrix
#genecolData <- data.frame(row.names=colnames(genecountData)[1:length(colnames(genecountData))], conditions)
genecolData <- data.frame(condition = conditions, row.names = sample_names)
genecountData <- genecountData[, rownames(genecolData)]

# pvalue and DEGs
pvalues <- sapply(1:nrow(genecountData),function(i){
    data<-cbind.data.frame(gene=as.numeric(t(genecountData[i,])),conditions)
    p=wilcox.test(gene ~ conditions, data)$p.value
     return(p)
})

fdr<-p.adjust(pvalues,method = "fdr")
conditionsLevel<-levels(conditions)
dataCon1 <- genecountData[,c(which(conditions==conditionsLevel[1]))]
dataCon2 <- genecountData[,c(which(conditions==conditionsLevel[2]))]
foldChanges=rowMeans(dataCon2) - rowMeans(dataCon1)
outRst<-data.frame(log2FoldChange=foldChanges, pValues=pvalues, FDR=fdr)

rownames(outRst)=rownames(genecountData)
#colnames(outRst)[1] <- 'GeneID'
outRst <- na.omit(outRst)

# Sort by q value
generesordered <- outRst[order(outRst$pValues),]

generesordered <- rownames_to_column(generesordered)
colnames(generesordered )[1] <- 'GeneID'
write.table(as.data.frame(generesordered),file = out_gene_tsv,row.names = F,sep="\t", quote=FALSE)

diff_gene_deseq2 <-subset(generesordered, pValues < pvalue_limit & abs(log2FoldChange) > log2(FC_limit))

up_regulated_genes <- subset(diff_gene_deseq2,log2FoldChange > 0)
down_regulated_genes <- subset(diff_gene_deseq2,log2FoldChange < 0)

write.table(diff_gene_deseq2,file = out_DE_gene,sep="\t", quote=FALSE, row.names=FALSE)
write.table(up_regulated_genes,file = out_up_gene,sep="\t", quote=FALSE, row.names=FALSE)
write.table(down_regulated_genes, file = out_down_gene,sep="\t", quote=FALSE, row.names=FALSE)
