# Title     : Identify DEGs for bio repeat samples (n<4) by DEGseq
# Created by: stl23
# Created on: 2021/8/5

library(DESeq2)
args <- commandArgs(trailingOnly = TRUE)
groups <- args[1]
input <- args[2]
path <- args[3]
prefix <- args[4]
p <- 0.05
fc <- 2 # log2(1.5) FC~1.5; log2(2) FC = 2

if (!dir.exists(paths = path)) {
      stop("Directory provided does not exist")
    }
#input_file <- file.path(path,'./before.all.wt.mut.combined.changeid.txt')
input_file <- file.path(input)
#group_name <- unlist(strsplit(groups,','))
group_name <- read.table(file = groups, sep="\t", header=T)
out_gene_tsv <- paste0(path,'/',prefix,'.RNAseq_gene_results.xls')
out_DE_gene <- paste0(path,'/',prefix,'.RNAseq_DEG_results.xls')
out_up_gene <- paste0(path,'/',prefix,'.RNAseq_UP_genes_results.xls')
out_down_gene <- paste0(path,'/',prefix,'.RNAseq_DOWN_genes_results.xls')

# Load gene data(Geneid file_name1 file_name2...)
genecountData <- read.table(file = input_file, sep="\t", row.names = 1, header=T,check.names = F)
## Data format (Geneid sample1 sample2 ...)
#genecountData <- genecountData[,-2:-6]
#names(genecountData[,-1]) <- sample_name ## change file names to sample names
#genecountData <- as.matrix(genecountData)
#rownames(genecountData) <- genecountData$Geneid
#colnames(genecountData)[1] <- 'GeneID'
condition <- factor(group_name$group,levels = unique(group_name$group))  ## order limited: c("Control","Treatment") --> Treatment vs. Control
id <- factor(group_name$patient,levels = unique(group_name$patient))
# Get target matrix
genecolData <- data.frame(row.names=colnames(genecountData)[1:length(colnames(genecountData))], condition,id)
genecountData <- genecountData[, rownames(genecolData)]
genedds <- DESeqDataSetFromMatrix(countData = as.matrix(round(genecountData)), colData = genecolData, design = ~ id+condition)
# Identify signficant differently expressed Transcripts/genes
genedataset <- DESeq(genedds)
generes <- as.data.frame(results(genedataset))
generes <- cbind(rownames(generes),generes)
colnames(generes)[1] <- 'GeneID'

# Sort by q value
generesordered <- generes[order(generes$pvalue),]
write.table(as.data.frame(generesordered),file = out_gene_tsv,row.names = F,sep="\t", quote=FALSE)

#generesordered$sig <- ifelse(generesordered$log2FoldChange>log2(fc)&generesordered$pvalue<p,
#                              "Up",ifelse(generesordered$log2FoldChange<(-log2(fc))&generesordered$pvalue<p,"Down","None"))
#write.table(as.data.frame(generesordered),file = out_sig,row.names = F,sep="\t", quote=FALSE)
# DEG output
diff_gene_deseq2 <-subset(generesordered, pvalue < p & abs(log2FoldChange) > log2(fc))
up_regulated_genes <- subset(diff_gene_deseq2,log2FoldChange > 0)
down_regulated_genes <- subset(diff_gene_deseq2,log2FoldChange < 0)


write.table(diff_gene_deseq2,file = out_DE_gene,sep="\t", quote=FALSE, row.names=FALSE)
write.table(up_regulated_genes,file = out_up_gene,sep="\t", quote=FALSE, row.names=FALSE)
write.table(down_regulated_genes, file = out_down_gene,sep="\t", quote=FALSE, row.names=FALSE)
