suppressPackageStartupMessages({
library(immunedeconv)
library(tinyarray)
library(tidyverse)
library(estimate)

library(preprocessCore)
library(e1071)
source("/Users/stl/Documents/AmoyDX/scripts/CIBERSORT.R")
library(MCPcounter)

library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(ggpubr)
library(pheatmap)
library(RColorBrewer)
library(ggsci)

library(openxlsx)
library(GSVA)
library(egg)
#library(xCell)
library(xCell2)
library(tidyplots)
})
## 画图部分写成函数，并且设定控制，可以跳过GSVA或ssGSEA，直接用前面生成的文件画图

draw_boxplot_v1 <- function(tool,group_color,group_comp,ncol,filename,ht,wd){
  p <- tool %>%
    tidyplot(x = group, y = Proportion, color = group) %>%
    add_boxplot(alpha = 0.5) %>%
    #add_data_points_beeswarm(size = 2, shape = 21, color = "black") %>%
    add_data_points_beeswarm(size = 1.5) %>%
    #adjust_x_axis(title = "Group") %>%
    remove_x_axis_labels() %>%
    remove_x_axis_title() %>%
    remove_x_axis_ticks() %>%
    adjust_y_axis_title(title = "") %>%
    adjust_size(width = NA,height = NA) %>%
    adjust_colors(new_colors = group_color) %>%
    add_test_pvalue(comparisons = group_comp, hide_info = TRUE, label = "p.format",method = "wilcox_test") %>%
    split_plot(by = Cell_type, ncol = ncol) 
  ggsave(p,filename = filename,height = ht,width = wd,family = "ArialMT")
}

draw_boxplot_v2 <- function(tool,group_color,method,filename,ht,wd,font_size){
  origin <- unique(tool$Origin)
  if(origin %in% c("Cibersort","Mcpcounter","xCell")){
    my_lab_y <- "Estimated Proportion"
  }else{
    my_lab_y <- "Estimated Score"
  }
  
  p <- ggplot(tool,aes(Cell_type,Proportion,fill = group)) +
    geom_boxplot(outlier.shape = 21,color = c("black")) +
    theme_article() +
    #coord_flip()+
    #ylim(-0.1,0.8)+
    labs(x = "Cell Type", y = my_lab_y) +
    theme(legend.position = "top") +
    theme(axis.text.x = element_text(size=8,face = 'bold',angle=45,vjust = 1, hjust = 1),
          axis.text.y = element_text(size=8,face = 'bold')
          #axis.line.x = element_line(linewidth=1,colour = 'black'),
          #axis.line.y = element_line(linewidth=1.5,colour = 'black')
    )+
    #scale_fill_manual(values = mypalette(22)[c(1,6)],
    scale_fill_manual(values = as.character(group_color))+
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))+
    stat_compare_means(aes(group = group,label = ..p.signif..),
                       method = method,size=font_size)
  p$layers[[2]]$aes_params$textsize <- font_size
  ggsave(p,filename = filename,height = ht,width = wd,family = "ArialMT")
}


run_all_immune_v2_file <- function(tpm_file,norm_file,
                                   group_file,
                                   xcell=FALSE){
  
  exprMatrix <- read.table(file=norm_file, header = T, sep = "\t",check.names=F,row.names = 1)
  groups <- read.table(file=group_file,header=T)
  my_comparisons <- list(unique(groups$group))
  exprMatrix <- exprMatrix[,groups$sample]
  expr_vst <- as.matrix(exprMatrix)
  
  font_size <- 3
  ######estimate####
  filterCommonGenes(input.f= norm_file,
                    output.f="commonGenes.gct",
                    id="GeneSymbol")
  
  estimateScore(input.ds = "commonGenes.gct",
                output.ds="estimateScore.gct")
  scores <- read.table("estimateScore.gct",skip = 2,header = T)
  rownames(scores)<-scores[,1]
  scores <- t(scores[,3:ncol(scores)])
  rownames(scores)<-gsub("\\.","\\-",rownames(scores))
  immune.score<-rbind(ID=colnames(scores),scores)
  write.table(immune.score,file="ESTIMATE_scores.txt",sep="\t",quote=F,col.names=F)
  
  estimate_data <- cbind(as.data.frame(scores), group = groups$group)
  estimate_data_trans <- estimate_data %>% as.data.frame() %>%
    rownames_to_column("ID") %>%
    gather(key = Cell_type,value = Proportion,-c(ID, group))
  ## let group factor keep control:case
  estimate_data_trans$group <- factor(estimate_data_trans$group,levels = unique(groups$group))
  
  ####cibersort####
  ciber_results <- CIBERSORT("/Users/stl/Documents/AmoyDX/scripts/LM22.txt", mixture_file = tpm_file,
                             perm=100, QN=FALSE)
  #QN如果是芯片数据设置为TRUE，如果是测序数据就设置为FALSE(测序数据最好是TPM)
  re <- ciber_results[,-(23:25)]  # remove last three rows of scores
  k <- apply(re,2,function(x) {sum(x == 0) < nrow(ciber_results)/2})
  re2 <- as.data.frame(t(re[,k]))
  colnames(re2)<-gsub("\\.","\\-",colnames(re2))
  anno <- data.frame(group = groups$group,
                     row.names = groups$sample)
  mypalette <- colorRampPalette(brewer.pal(8,"Set1"))
  ciber_dat <- re %>% as.data.frame() %>%
    mutate(group = groups$group) %>%
    rownames_to_column("Sample") %>%
    gather(key = Cell_type,value = Proportion,-c(Sample, group))
  ciber_dat$Sample<-gsub("\\.","\\-",ciber_dat$Sample)
  ciber_dat$Origin <- "Cibersort"
  
  ###mcpcounter###
  genes <- read.table("/Users/stl/Documents/AmoyDX/scripts/genes.txt",
                      header = T,sep = "\t", check.names = FALSE)
  
  mcp_results<- MCPcounter.estimate(exprMatrix,
                                    featuresType= "HUGO_symbols",
                                    #probesets=probesets,
                                    genes=genes)
  results2 <- mcp_results
  colnames(results2) <- gsub("\\.", "-", colnames(results2))
  mcp_dat <- results2  %>% t() %>% as.data.frame() %>%
    mutate(group = groups$group) %>%
    rownames_to_column("Sample") %>%
    gather(key = Cell_type,value = Proportion,-c(Sample, group))
  mcp_dat$Origin <- "Mcpcounter"
  
  ## AO
  AO <- read.xlsx("/Users/stl/Documents/AmoyDX/scripts/Signatures/Signature/AO/ao.pathway.xlsx",sheet=5)
  AO_datalist <- lapply((AO), function(x) x[!is.na(x)])
  #ao.exp <- gsva(expr_vst, AO_datalist, method = "gsva")
  #gsvaPar <- gsvaParam(expr_vst, AO_datalist)
  #ao.exp <- gsva(gsvaPar)
  ao.gpar <- gsvaParam(exprData = expr_vst,geneSets = AO_datalist,kcdf = "Gaussian")
  ao.exp <- gsva(ao.gpar,verbose=TRUE)
  
  ao_dat <- ao.exp %>% t() %>% as.data.frame() %>%
    rownames_to_column("Sample") %>% mutate(group = groups$group) %>%
    gather(key = Cell_type,value = Proportion,-c(Sample, group))
  ao_dat$Origin <- "AO"
  
  #ao.exp.ss <- gsva(expr_vst, AO_datalist, method = "ssgsea")
  ao.gpar.ss <- ssgseaParam(exprData = expr_vst,geneSets = AO_datalist)
  ao.exp.ss <- gsva(ao.gpar.ss,verbose=TRUE)
  
  ao_dat.ss <- ao.exp.ss %>% t() %>% as.data.frame() %>%
    rownames_to_column("Sample") %>% mutate(group = groups$group) %>%
    gather(key = Cell_type,value = Proportion,-c(Sample, group))
  ao_dat.ss$Origin <- "AO"
  
  ## Danaher
  danaher = read.xlsx("/Users/stl/Documents/AmoyDX/scripts/Signatures/Signature/14immcell_sign/Danaher_Gene_expression_markers_of_Tumor_Infiltrating_Leukocytes.xlsx")
  dan_list = split(danaher[,1], danaher[,2])
  #dan.exp <- gsva(expr_vst, dan_list, method = "gsva")
  dan.gpar <- gsvaParam(exprData = expr_vst,geneSets = dan_list,kcdf = "Gaussian")
  dan.exp <- gsva(dan.gpar,verbose=TRUE)
  
  dan_dat <- dan.exp %>% t() %>% as.data.frame() %>%
    rownames_to_column("Sample") %>% mutate(group = groups$group) %>%
    gather(key = Cell_type,value = Proportion,-c(Sample, group))
  dan_dat$Origin <- "Danaher"
  
  #dan.exp.ss <- gsva(expr_vst, dan_list, method = "ssgsea")
  dan.gpar.ss <- ssgseaParam(exprData = expr_vst,geneSets = dan_list)
  dan.exp.ss <- gsva(dan.gpar.ss,verbose=TRUE)
  dan_dat.ss <- dan.exp.ss %>% t() %>% as.data.frame() %>%
    rownames_to_column("Sample") %>% mutate(group = groups$group) %>%
    gather(key = Cell_type,value = Proportion,-c(Sample, group))
  dan_dat.ss$Origin <- "Danaher"
  
  ## immune 29 genes ##
  immune29 <- read.xlsx("/Users/stl/Documents/AmoyDX/scripts/immune.signature.xlsx")
  immune29_datalist <- lapply((immune29), function(x) x[!is.na(x)])
  #gs.exp <- gsva(expr_vst, immune29_datalist, method = "gsva")
  gs.gpar <- gsvaParam(exprData = expr_vst,geneSets = immune29_datalist,kcdf = "Gaussian")
  gs.exp <- gsva(gs.gpar,verbose=TRUE)
  
  immune_dat29 <- gs.exp %>% t() %>% as.data.frame() %>%
    rownames_to_column("Sample") %>% mutate(group = groups$group) %>%
    gather(key = Cell_type,value = Proportion,-c(Sample, group))
  immune_dat29$Origin <- "Immune29"
  
  #gs.exp.ss <- gsva(expr_vst, immune29_datalist, method = "ssgsea")
  gs.gpar.ss <- ssgseaParam(exprData = expr_vst,geneSets = immune29_datalist)
  gs.exp.ss <- gsva(gs.gpar.ss,verbose=TRUE)
  
  immune_dat29.ss <- gs.exp.ss %>% t() %>% as.data.frame() %>%
    rownames_to_column("Sample") %>% mutate(group = groups$group) %>%
    gather(key = Cell_type,value = Proportion,-c(Sample, group))
  immune_dat29.ss$Origin <- "Immune29"
  
  ##### immune28
  immune28 <- read.csv("/Users/stl/Documents/AmoyDX/scripts/28immcells.csv",header=TRUE,sep=",",na.strings = "")
  immune28_datalist <- lapply((immune28), function(x) x[!is.na(x)])
  #gs.exp3 <- gsva(expr_vst, immune28_datalist, method = "gsva")
  gs.gpar3 <- gsvaParam(exprData = expr_vst,geneSets = immune28_datalist,kcdf = "Gaussian")
  gs.exp3 <- gsva(gs.gpar3,verbose=TRUE)
  
  immune_dat28 <- gs.exp3 %>% t() %>% as.data.frame() %>%
    rownames_to_column("Sample") %>% mutate(group = groups$group) %>%
    gather(key = Cell_type,value = Proportion,-c(Sample, group))
  immune_dat28$Origin <- "Immune28"
  
  #gs.exp3.ss <- gsva(expr_vst, immune28_datalist, method = "ssgsea")
  gs.gpar3.ss <- ssgseaParam(exprData = expr_vst,geneSets = immune28_datalist)
  gs.exp3.ss <- gsva(gs.gpar3.ss,verbose=TRUE)
  immune_dat28.ss <- gs.exp3.ss %>% t() %>% as.data.frame() %>%
    rownames_to_column("Sample") %>% mutate(group = groups$group) %>%
    gather(key = Cell_type,value = Proportion,-c(Sample, group))
  immune_dat28.ss$Origin <- "Immune28"
  
  #### IFN
  IFN <- read.xlsx("/Users/stl/Documents/AmoyDX/scripts/IFN.xlsx")
  IFN_datalist <- lapply((IFN), function(x) x[!is.na(x)])
  #gs.exp2 <- gsva(expr_vst, IFN_datalist, method = "gsva")
  gs.gpar2 <- gsvaParam(exprData = expr_vst,geneSets = IFN_datalist,kcdf = "Gaussian")
  gs.exp2 <- gsva(gs.gpar2,verbose=TRUE)
  IFN_dat <- gs.exp2 %>% t() %>% as.data.frame() %>%
    rownames_to_column("Sample") %>% mutate(group = groups$group) %>%
    gather(key = Cell_type,value = Proportion,-c(Sample, group))
  IFN_dat$Origin <- "IFN"
  
  #gs.exp2.ss <- gsva(expr_vst, IFN_datalist, method = "ssgsea")
  gs.gpar2.ss <- ssgseaParam(exprData = expr_vst,geneSets = IFN_datalist)
  gs.exp2.ss <- gsva(gs.gpar2.ss,verbose=TRUE)
  IFN_dat.ss <- gs.exp2.ss %>% t() %>% as.data.frame() %>%
    rownames_to_column("Sample") %>% mutate(group = groups$group) %>%
    gather(key = Cell_type,value = Proportion,-c(Sample, group))
  IFN_dat.ss$Origin <- "IFN"
  
  ### CAFs
  cafs <- read.csv("/Users/stl/Documents/AmoyDX/scripts/Signatures/Signature/CAF_subtype_sign/CAF_subtype_signature.csv")
  cafs_list <- lapply((cafs), function(x) x[!is.na(x) & x != ""])
  
  #caf.exp <- gsva(expr_vst, cafs_list, method = "gsva")
  caf.gpar <- gsvaParam(exprData = expr_vst,geneSets = cafs_list,kcdf = "Gaussian")
  caf.exp <- gsva(caf.gpar,verbose=TRUE)
  cafs_data <- caf.exp %>% t() %>% as.data.frame() %>%
    rownames_to_column("Sample") %>% mutate(group = groups$group) %>%
    gather(key = Cell_type,value = Proportion,-c(Sample, group))
  cafs_data$Origin <- "CAFs"
  
  #caf.exp2 <- gsva(expr_vst, cafs_list, method = "ssgsea")
  caf.gpar.ss <- ssgseaParam(exprData = expr_vst,geneSets = cafs_list)
  caf.exp.ss <- gsva(caf.gpar.ss,verbose=TRUE)
  cafs_data.ss <- caf.exp.ss %>% t() %>% as.data.frame() %>%
    rownames_to_column("Sample") %>% mutate(group = groups$group) %>%
    gather(key = Cell_type,value = Proportion,-c(Sample, group))
  cafs_data.ss$Origin <- "CAFs"
  
  all_dat <- rbind(ciber_dat,mcp_dat,ao_dat,dan_dat,immune_dat29,immune_dat28,IFN_dat,cafs_data)
  
  #### xCell
  if(xcell){   ## xcell use ssGSEA
    #xcell.result <- xCellAnalysis(exprMatrix,rnaseq = T)
    # 更新GSVA后，xcell不能用了，替换为xCell2
    data("BlueprintEncode.xCell2Ref") # https://aran-lab.com/xcell2-vignette
    xcell.result <- xCell2::xCell2Analysis(mix = exprMatrix,xcell2object = BlueprintEncode.xCell2Ref,
                                           minSharedGenes = 0.7)
    xcell_dat <- xcell.result  %>% t() %>% as.data.frame() %>%
      mutate(group = groups$group) %>%
      rownames_to_column("Sample") %>%
      gather(key = Cell_type,value = Proportion,-c(Sample, group))
    xcell_dat$Origin <- "xCell"
    all_dat.ss <- rbind(ao_dat.ss,dan_dat.ss,immune_dat29.ss,immune_dat28.ss,IFN_dat.ss,cafs_data.ss,xcell_dat)
  }else{
    all_dat.ss <- rbind(ao_dat.ss,dan_dat.ss,immune_dat29.ss,immune_dat28.ss,IFN_dat.ss,cafs_data.ss)
  }
  
  all_dat$group <- factor(all_dat$group,levels = unique(groups$group))
  all_dat.ss$group <- factor(all_dat.ss$group,levels = unique(groups$group))
  write.table(all_dat,file="all_gsva_scores.txt",sep="\t",quote=F,row.names=F)
  write.table(all_dat.ss,file="all_ssgsea_scores.txt",sep="\t",quote=F,row.names=F)
  
  # 保存
  if(xcell){ 
    res_list <- list(
      # 基础输入变量
      input = list(tpm_file = tpm_file, norm_file = norm_file, group_file = group_file,
                   exprMatrix = exprMatrix, groups = groups, my_comparisons = my_comparisons, expr_vst = expr_vst),
      # ESTIMATE分析结果
      estimate = list(scores = scores, immune.score = immune.score,
                      estimate_data = estimate_data, estimate_data_trans = estimate_data_trans),
      # CIBERSORT分析结果
      cibersort = list(ciber_results = ciber_results, re = re, re2 = re2, ciber_dat = ciber_dat),
      # MCPcounter分析结果
      mcpcounter = list(genes = genes, mcp_results = mcp_results, results2 = results2, mcp_dat = mcp_dat),
      # 各类GSVA/ssGSEA分析结果（AO/Danaher/免疫基因/IFN/CAFs）
      gsva_ssgsea = list(AO = AO, ao.exp = ao.exp, ao_dat = ao_dat, ao.exp.ss = ao.exp.ss, ao_dat.ss = ao_dat.ss,
                         danaher = danaher, dan.exp = dan.exp, dan_dat = dan_dat, dan.exp.ss = dan.exp.ss, dan_dat.ss = dan_dat.ss,
                         immune29 = immune29, immune_dat29 = immune_dat29, immune_dat29.ss = immune_dat29.ss,
                         immune28 = immune28, immune_dat28 = immune_dat28, immune_dat28.ss = immune_dat28.ss,
                         IFN = IFN, IFN_dat = IFN_dat, IFN_dat.ss = IFN_dat.ss,
                         cafs = cafs, cafs_data = cafs_data, cafs_data.ss = cafs_data.ss),
      # xCell分析结果（可选）
      xcell = list(xcell.result = xcell.result, xcell_dat = xcell_dat),
      # 合并后的最终结果
      final = list(all_dat = all_dat, all_dat.ss = all_dat.ss)
    )}else{
      res_list <- list(
        # 基础输入变量
        input = list(tpm_file = tpm_file, norm_file = norm_file, group_file = group_file,
                     exprMatrix = exprMatrix, groups = groups, my_comparisons = my_comparisons, expr_vst = expr_vst),
        # ESTIMATE分析结果
        estimate = list(scores = scores, immune.score = immune.score,
                        estimate_data = estimate_data, estimate_data_trans = estimate_data_trans),
        # CIBERSORT分析结果
        cibersort = list(ciber_results = ciber_results, re = re, re2 = re2, ciber_dat = ciber_dat),
        # MCPcounter分析结果
        mcpcounter = list(genes = genes, mcp_results = mcp_results, results2 = results2, mcp_dat = mcp_dat),
        # 各类GSVA/ssGSEA分析结果（AO/Danaher/免疫基因/IFN/CAFs）
        gsva_ssgsea = list(AO = AO, ao.exp = ao.exp, ao_dat = ao_dat, ao.exp.ss = ao.exp.ss, ao_dat.ss = ao_dat.ss,
                           danaher = danaher, dan.exp = dan.exp, dan_dat = dan_dat, dan.exp.ss = dan.exp.ss, dan_dat.ss = dan_dat.ss,
                           immune29 = immune29, immune_dat29 = immune_dat29, immune_dat29.ss = immune_dat29.ss,
                           immune28 = immune28, immune_dat28 = immune_dat28, immune_dat28.ss = immune_dat28.ss,
                           IFN = IFN, IFN_dat = IFN_dat, IFN_dat.ss = IFN_dat.ss,
                           cafs = cafs, cafs_data = cafs_data, cafs_data.ss = cafs_data.ss),
        # 合并后的最终结果
        final = list(all_dat = all_dat, all_dat.ss = all_dat.ss)
      )
    }
  
  ### 函数返回该列表
  return(res_list)
}

## Olink等panel 测序做不了xCell，以及没有TPM值做不了cibersort等反卷积
run_part_immune_file <- function(norm_file,
                                 group_file,
                                 xcell=FALSE){
  # 判断是文件还是矩阵
  if (is.matrix(norm_file) || is.data.frame(norm_file)){
    #warning("data may not be a file,treat it as matrix/data.frame")
    exprMatrix <- norm_file # row：基因ID； col： 样本
    tmp <- tempfile(fileext = ".txt")
    write.table(norm_file, tmp, sep="\t", row.names=T, quote=F)
    input_path <- tmp
  }else if(is.character(norm_file) && length(norm_file) == 1 && file.exists(norm_file)){
    exprMatrix <- read.table(file=norm_file, header = T, sep = "\t",check.names=F,row.names = 1)
    input_path <- norm_file
  }else{
    stop("norm_file 必须是：存在的文件路径 或 矩阵/数据框")
  }
  
  groups <- read.table(file=group_file,header=T)
  my_comparisons <- list(unique(groups$group))
  exprMatrix <- exprMatrix[,groups$sample]
  expr_vst <- as.matrix(exprMatrix)
  
  font_size <- 3
  
  ## AO
  AO <- read.xlsx("/Users/stl/Documents/AmoyDX/scripts/Signatures/Signature/AO/ao.pathway.xlsx",sheet=5)
  AO_datalist <- lapply((AO), function(x) x[!is.na(x)])
  #ao.exp <- gsva(expr_vst, AO_datalist, method = "gsva")
  #gsvaPar <- gsvaParam(expr_vst, AO_datalist)
  #ao.exp <- gsva(gsvaPar)
  ao.gpar <- gsvaParam(exprData = expr_vst,geneSets = AO_datalist,kcdf = "Gaussian")
  ao.exp <- gsva(ao.gpar,verbose=TRUE)
  
  ao_dat <- ao.exp %>% t() %>% as.data.frame() %>%
    rownames_to_column("Sample") %>% mutate(group = groups$group) %>%
    gather(key = Cell_type,value = Proportion,-c(Sample, group))
  ao_dat$Origin <- "AO"
  
  #ao.exp.ss <- gsva(expr_vst, AO_datalist, method = "ssgsea")
  ao.gpar.ss <- ssgseaParam(exprData = expr_vst,geneSets = AO_datalist)
  ao.exp.ss <- gsva(ao.gpar.ss,verbose=TRUE)
  
  ao_dat.ss <- ao.exp.ss %>% t() %>% as.data.frame() %>%
    rownames_to_column("Sample") %>% mutate(group = groups$group) %>%
    gather(key = Cell_type,value = Proportion,-c(Sample, group))
  ao_dat.ss$Origin <- "AO"
  
  ## Danaher
  danaher = read.xlsx("/Users/stl/Documents/AmoyDX/scripts/Signatures/Signature/14immcell_sign/Danaher_Gene_expression_markers_of_Tumor_Infiltrating_Leukocytes.xlsx")
  dan_list = split(danaher[,1], danaher[,2])
  #dan.exp <- gsva(expr_vst, dan_list, method = "gsva")
  dan.gpar <- gsvaParam(exprData = expr_vst,geneSets = dan_list,kcdf = "Gaussian")
  dan.exp <- gsva(dan.gpar,verbose=TRUE)
  
  dan_dat <- dan.exp %>% t() %>% as.data.frame() %>%
    rownames_to_column("Sample") %>% mutate(group = groups$group) %>%
    gather(key = Cell_type,value = Proportion,-c(Sample, group))
  dan_dat$Origin <- "Danaher"
  
  #dan.exp.ss <- gsva(expr_vst, dan_list, method = "ssgsea")
  dan.gpar.ss <- ssgseaParam(exprData = expr_vst,geneSets = dan_list)
  dan.exp.ss <- gsva(dan.gpar.ss,verbose=TRUE)
  dan_dat.ss <- dan.exp.ss %>% t() %>% as.data.frame() %>%
    rownames_to_column("Sample") %>% mutate(group = groups$group) %>%
    gather(key = Cell_type,value = Proportion,-c(Sample, group))
  dan_dat.ss$Origin <- "Danaher"
  
  ## immune 29 genes ##
  immune29 <- read.xlsx("/Users/stl/Documents/AmoyDX/scripts/immune.signature.xlsx")
  immune29_datalist <- lapply((immune29), function(x) x[!is.na(x)])
  #gs.exp <- gsva(expr_vst, immune29_datalist, method = "gsva")
  gs.gpar <- gsvaParam(exprData = expr_vst,geneSets = immune29_datalist,kcdf = "Gaussian")
  gs.exp <- gsva(gs.gpar,verbose=TRUE)
  
  immune_dat29 <- gs.exp %>% t() %>% as.data.frame() %>%
    rownames_to_column("Sample") %>% mutate(group = groups$group) %>%
    gather(key = Cell_type,value = Proportion,-c(Sample, group))
  immune_dat29$Origin <- "Immune29"
  
  #gs.exp.ss <- gsva(expr_vst, immune29_datalist, method = "ssgsea")
  gs.gpar.ss <- ssgseaParam(exprData = expr_vst,geneSets = immune29_datalist)
  gs.exp.ss <- gsva(gs.gpar.ss,verbose=TRUE)
  
  immune_dat29.ss <- gs.exp.ss %>% t() %>% as.data.frame() %>%
    rownames_to_column("Sample") %>% mutate(group = groups$group) %>%
    gather(key = Cell_type,value = Proportion,-c(Sample, group))
  immune_dat29.ss$Origin <- "Immune29"
  
  ##### immune28
  immune28 <- read.csv("/Users/stl/Documents/AmoyDX/scripts/28immcells.csv",header=TRUE,sep=",",na.strings = "")
  immune28_datalist <- lapply((immune28), function(x) x[!is.na(x)])
  #gs.exp3 <- gsva(expr_vst, immune28_datalist, method = "gsva")
  gs.gpar3 <- gsvaParam(exprData = expr_vst,geneSets = immune28_datalist,kcdf = "Gaussian")
  gs.exp3 <- gsva(gs.gpar3,verbose=TRUE)
  
  immune_dat28 <- gs.exp3 %>% t() %>% as.data.frame() %>%
    rownames_to_column("Sample") %>% mutate(group = groups$group) %>%
    gather(key = Cell_type,value = Proportion,-c(Sample, group))
  immune_dat28$Origin <- "Immune28"
  
  #gs.exp3.ss <- gsva(expr_vst, immune28_datalist, method = "ssgsea")
  gs.gpar3.ss <- ssgseaParam(exprData = expr_vst,geneSets = immune28_datalist)
  gs.exp3.ss <- gsva(gs.gpar3.ss,verbose=TRUE)
  immune_dat28.ss <- gs.exp3.ss %>% t() %>% as.data.frame() %>%
    rownames_to_column("Sample") %>% mutate(group = groups$group) %>%
    gather(key = Cell_type,value = Proportion,-c(Sample, group))
  immune_dat28.ss$Origin <- "Immune28"
  
  #### IFN
  IFN <- read.xlsx("/Users/stl/Documents/AmoyDX/scripts/IFN.xlsx")
  IFN_datalist <- lapply((IFN), function(x) x[!is.na(x)])
  #gs.exp2 <- gsva(expr_vst, IFN_datalist, method = "gsva")
  gs.gpar2 <- gsvaParam(exprData = expr_vst,geneSets = IFN_datalist,kcdf = "Gaussian")
  gs.exp2 <- gsva(gs.gpar2,verbose=TRUE)
  IFN_dat <- gs.exp2 %>% t() %>% as.data.frame() %>%
    rownames_to_column("Sample") %>% mutate(group = groups$group) %>%
    gather(key = Cell_type,value = Proportion,-c(Sample, group))
  IFN_dat$Origin <- "IFN"
  
  #gs.exp2.ss <- gsva(expr_vst, IFN_datalist, method = "ssgsea")
  gs.gpar2.ss <- ssgseaParam(exprData = expr_vst,geneSets = IFN_datalist)
  gs.exp2.ss <- gsva(gs.gpar2.ss,verbose=TRUE)
  IFN_dat.ss <- gs.exp2.ss %>% t() %>% as.data.frame() %>%
    rownames_to_column("Sample") %>% mutate(group = groups$group) %>%
    gather(key = Cell_type,value = Proportion,-c(Sample, group))
  IFN_dat.ss$Origin <- "IFN"
  
  ### CAFs
  cafs <- read.csv("/Users/stl/Documents/AmoyDX/scripts/Signatures/Signature/CAF_subtype_sign/CAF_subtype_signature.csv")
  cafs_list <- lapply((cafs), function(x) x[!is.na(x) & x != ""])
  
  #caf.exp <- gsva(expr_vst, cafs_list, method = "gsva")
  caf.gpar <- gsvaParam(exprData = expr_vst,geneSets = cafs_list,kcdf = "Gaussian")
  caf.exp <- gsva(caf.gpar,verbose=TRUE)
  cafs_data <- caf.exp %>% t() %>% as.data.frame() %>%
    rownames_to_column("Sample") %>% mutate(group = groups$group) %>%
    gather(key = Cell_type,value = Proportion,-c(Sample, group))
  cafs_data$Origin <- "CAFs"
  
  #caf.exp2 <- gsva(expr_vst, cafs_list, method = "ssgsea")
  caf.gpar.ss <- ssgseaParam(exprData = expr_vst,geneSets = cafs_list)
  caf.exp.ss <- gsva(caf.gpar.ss,verbose=TRUE)
  cafs_data.ss <- caf.exp.ss %>% t() %>% as.data.frame() %>%
    rownames_to_column("Sample") %>% mutate(group = groups$group) %>%
    gather(key = Cell_type,value = Proportion,-c(Sample, group))
  cafs_data.ss$Origin <- "CAFs"
  
  all_dat <- rbind(ao_dat,dan_dat,immune_dat29,immune_dat28,IFN_dat,cafs_data)
  
  #### xCell
  if(xcell){   ## xcell use ssGSEA
    #xcell.result <- xCellAnalysis(exprMatrix,rnaseq = T)
    # 更新GSVA后，xcell不能用了，替换为xCell2
    data("BlueprintEncode.xCell2Ref") # https://aran-lab.com/xcell2-vignette
    xcell.result <- xCell2::xCell2Analysis(mix = exprMatrix,xcell2object = BlueprintEncode.xCell2Ref,
                                           minSharedGenes = 0.7)
    xcell_dat <- xcell.result  %>% t() %>% as.data.frame() %>%
      mutate(group = groups$group) %>%
      rownames_to_column("Sample") %>%
      gather(key = Cell_type,value = Proportion,-c(Sample, group))
    xcell_dat$Origin <- "xCell"
    all_dat.ss <- rbind(ao_dat.ss,dan_dat.ss,immune_dat29.ss,immune_dat28.ss,IFN_dat.ss,cafs_data.ss,xcell_dat)
  }else{
    all_dat.ss <- rbind(ao_dat.ss,dan_dat.ss,immune_dat29.ss,immune_dat28.ss,IFN_dat.ss,cafs_data.ss)
  }
  
  all_dat$group <- factor(all_dat$group,levels = unique(groups$group))
  all_dat.ss$group <- factor(all_dat.ss$group,levels = unique(groups$group))
  write.table(all_dat,file="all_gsva_scores.txt",sep="\t",quote=F,row.names=F)
  write.table(all_dat.ss,file="all_ssgsea_scores.txt",sep="\t",quote=F,row.names=F)
  
  # 保存
  if(xcell){ 
    res_list <- list(
      # 基础输入变量
      input = list(tpm_file = tpm_file, norm_file = norm_file, group_file = group_file,
                   exprMatrix = exprMatrix, groups = groups, my_comparisons = my_comparisons, expr_vst = expr_vst),
      # # ESTIMATE分析结果
      # estimate = list(scores = scores, immune.score = immune.score,
      #                 estimate_data = estimate_data, estimate_data_trans = estimate_data_trans),
      
      # MCPcounter分析结果
      #mcpcounter = list(genes = genes, mcp_results = mcp_results, results2 = results2, mcp_dat = mcp_dat),
      # 各类GSVA/ssGSEA分析结果（AO/Danaher/免疫基因/IFN/CAFs）
      gsva_ssgsea = list(AO = AO, ao.exp = ao.exp, ao_dat = ao_dat, ao.exp.ss = ao.exp.ss, ao_dat.ss = ao_dat.ss,
                         danaher = danaher, dan.exp = dan.exp, dan_dat = dan_dat, dan.exp.ss = dan.exp.ss, dan_dat.ss = dan_dat.ss,
                         immune29 = immune29, immune_dat29 = immune_dat29, immune_dat29.ss = immune_dat29.ss,
                         immune28 = immune28, immune_dat28 = immune_dat28, immune_dat28.ss = immune_dat28.ss,
                         IFN = IFN, IFN_dat = IFN_dat, IFN_dat.ss = IFN_dat.ss,
                         cafs = cafs, cafs_data = cafs_data, cafs_data.ss = cafs_data.ss),
      # xCell分析结果（可选）
      xcell = list(xcell.result = xcell.result, xcell_dat = xcell_dat),
      # 合并后的最终结果
      final = list(all_dat = all_dat, all_dat.ss = all_dat.ss)
    )}else{
      res_list <- list(
        # 基础输入变量
        input = list(norm_file = norm_file, group_file = group_file,
                     exprMatrix = exprMatrix, groups = groups, my_comparisons = my_comparisons, expr_vst = expr_vst),
        # # ESTIMATE分析结果
        # estimate = list(scores = scores, immune.score = immune.score,
        #                 estimate_data = estimate_data, estimate_data_trans = estimate_data_trans),
        # MCPcounter分析结果
        #mcpcounter = list(genes = genes, mcp_results = mcp_results, results2 = results2, mcp_dat = mcp_dat),
        # 各类GSVA/ssGSEA分析结果（AO/Danaher/免疫基因/IFN/CAFs）
        gsva_ssgsea = list(AO = AO, ao.exp = ao.exp, ao_dat = ao_dat, ao.exp.ss = ao.exp.ss, ao_dat.ss = ao_dat.ss,
                           danaher = danaher, dan.exp = dan.exp, dan_dat = dan_dat, dan.exp.ss = dan.exp.ss, dan_dat.ss = dan_dat.ss,
                           immune29 = immune29, immune_dat29 = immune_dat29, immune_dat29.ss = immune_dat29.ss,
                           immune28 = immune28, immune_dat28 = immune_dat28, immune_dat28.ss = immune_dat28.ss,
                           IFN = IFN, IFN_dat = IFN_dat, IFN_dat.ss = IFN_dat.ss,
                           cafs = cafs, cafs_data = cafs_data, cafs_data.ss = cafs_data.ss),
        # 合并后的最终结果
        final = list(all_dat = all_dat, all_dat.ss = all_dat.ss)
      )
    }
  
  ### 函数返回该列表
  return(res_list)
}

run_all_immune_v2_draw <- function(res_list,group_comp=list(c(1,2)),group_color=colors_discrete_friendly,
                                   xcell=FALSE,draw_ssgsea=FALSE){
  mypalette <- colorRampPalette(brewer.pal(8,"Set1"))
  # tme_file <- all_file
  # if (!file.exists(tme_file)){
  #   stop(paste("文件未找到，程序停止:", tme_file))
  # }
  # load(tme_file)
  
  font_size = 3
  ####### draw
  # estimate
  if("estimate" %in% names(res_list)){
    estimate_data_trans <- res_list$estimate$estimate_data_trans
    draw_boxplot_v1(estimate_data_trans,group_color = group_color,group_comp = group_comp,ncol = 2,filename = "1.Estimate.boxplot.pdf",
                    ht=5,wd = 7)
  }
  all_dat <- res_list$final$all_dat
  # cibersort
  for(i in unique(all_dat$Origin)){
    tool <- all_dat[all_dat$Origin == i,]
    
    # if(length(na.omit(tool$Proportion)) < 3){
    #   stop("样本量小于3，无法进行正态性检验")
    # }else if(length(na.omit(tool$Proportion)) < 5000){
    #   pvalue <- shapiro.test(tool$Proportion)$p.value
    # }else{
    #   pvalue <- ad.test(tool$Proportion)$p.value ## Anderson-Darling检
    # }
    # 
    # if(pvalue > 0.05){  ## 符合正态分布
    #   method = "t_test"
    # }else{
    #   method = "wilcox.test"
    # }
    
    method = "wilcox.test"
    
    if(i == "Cibersort"){
      # p2.1 <- pheatmap(re2,scale = "row",
      #                  show_colnames = F,
      #                  cluster_cols = F,
      #                  annotation_col = anno,
      #                  color = colorRampPalette(c("navy", "white", "firebrick3"))(50))
      # pdf("2.cibersort.pheatmap.pdf",height = 5,width = 5,family = "ArialMT")
      # print(p2.1)
      # dev.off()
      
      ps <- ggplot(tool,aes(Sample,Proportion,fill = Cell_type)) +
        geom_bar(stat = "identity") +
        labs(fill = "Cell Type",x = "",y = "Estiamted Proportion") +
        theme_bw() +
        theme(axis.text.x = element_blank(),
              axis.ticks.x = element_blank(),
              legend.position = "bottom") +
        scale_y_continuous(expand = c(0.01,0)) +
        scale_fill_manual(values = mypalette(22))
      pdf("2.cibersort.sample.proportion.pdf",family = "ArialMT")
      print(ps)
      dev.off()
      
      pc <- ggplot(tool,aes(Cell_type,Proportion,fill = Cell_type)) +
        geom_boxplot(outlier.shape = 21,color = "black") +
        theme_bw() +
        labs(x = "Cell Type", y = "Estimated Proportion") +
        theme(axis.text.x = element_blank(),
              axis.ticks.x = element_blank(),
              legend.position = "bottom") +
        scale_fill_manual(values = mypalette(22))
      pdf("2.cibersort.celltype.proportion.pdf",family = "ArialMT")
      print(pc)
      dev.off()
      
      
      draw_boxplot_v1(tool,group_color = group_color,group_comp = group_comp,ncol = 5,filename = "2.cibersort_boxplot_v1.pdf",
                      ht=14,wd = 14)
      draw_boxplot_v2(tool,group_color = group_color,method = method,filename = "2.cibersort_boxplot_v2.pdf",
                      ht=3.5,wd = 12,font_size = font_size)
      
      
    }else if (i == "Mcpcounter"){
      
      draw_boxplot_v1(tool,group_color = group_color,group_comp = group_comp,ncol = 3,filename = "3.mcpcounter_boxplot_v1.pdf",
                      ht=9,wd = 9)
      draw_boxplot_v2(tool,group_color = group_color,method = method,filename = "3.mcpcounter_boxplot_v2.pdf",
                      ht=3.5,wd = 7,font_size = font_size)
      
      
    }else if (i == "AO"){
      
      draw_boxplot_v1(tool,group_color = group_color,group_comp = group_comp,ncol = 6,filename = "4.ao_boxplot_v1.pdf",
                      ht=12,wd = 15)
      draw_boxplot_v2(tool,group_color = group_color,method = method,filename = "4.ao_boxplot_v2.pdf",
                      ht=3.5,wd = 12,font_size = font_size)
      
      
    }else if(i == "Danaher"){
      draw_boxplot_v1(tool,group_color = group_color,group_comp = group_comp,ncol = 4,filename = "5.danaher_boxplot_v1.pdf",
                      ht=9,wd = 12)
      draw_boxplot_v2(tool,group_color = group_color,method = method,filename = "5.danaher_boxplot_v2.pdf",
                      ht=3,wd = 6.5,font_size = font_size)
      
    }else if(i == "Immune29"){
      draw_boxplot_v1(tool,group_color = group_color,group_comp = group_comp,ncol = 6,filename = "6.immune29_boxplot_v1.pdf",
                      ht=14,wd = 16)
      draw_boxplot_v2(tool,group_color = group_color,method = method,filename = "6.immune29_boxplot_v2.pdf",
                      ht=3.5,wd = 14,font_size = font_size)
      
    }else if(i == "Immune28"){
      draw_boxplot_v1(tool,group_color = group_color,group_comp = group_comp,ncol = 6,filename = "7.immune28_boxplot_v1.pdf",
                      ht=14,wd = 16)
      draw_boxplot_v2(tool,group_color = group_color,method = method,filename = "7.immune28_boxplot_v2.pdf",
                      ht=3.5,wd = 12,font_size = font_size)
      
      
    }else if(i == "IFN"){
      draw_boxplot_v1(tool,group_color = group_color,group_comp = group_comp,ncol = 3,filename = "8.IFN_boxplot_v1.pdf",
                      ht=8,wd = 9)
      draw_boxplot_v2(tool,group_color = group_color,method = method,filename = "8.IFN_boxplot_v2.pdf",
                      ht=3.5,wd = 6,font_size = font_size)
      
      
    }else if(i == "xCell"){
      draw_boxplot_v1(tool,group_color = group_color,group_comp = group_comp,ncol = 6,filename = "9.xCell_boxplot_v1.pdf",
                      ht=18,wd = 20)
      draw_boxplot_v2(tool,group_color = group_color,method = method,filename = "9.xCell_boxplot_v2.pdf",
                      ht=3.5,wd = 16,font_size = font_size)
      
    }else if(i == "CAFs"){
      draw_boxplot_v1(tool,group_color = group_color,group_comp = group_comp,ncol = 3,filename = "10.CAFs_boxplot_v1.pdf",
                      ht=5,wd = 9)
      draw_boxplot_v2(tool,group_color = group_color,method = method,filename = "10.CAFs_boxplot_v2.pdf",
                      ht=3.5,wd = 5,font_size = font_size)
      
      
    }
  }
  
  
  if(draw_ssgsea){
    all_dat.ss <- res_list$final$all_dat.ss
    for(i in unique(all_dat.ss$Origin)){
      tool <- all_dat.ss[all_dat.ss$Origin == i,]
      
      # if(length(na.omit(tool$Proportion)) < 3){
      #   stop("样本量小于3，无法进行正态性检验")
      # }else if(length(na.omit(tool$Proportion)) < 5000){
      #   pvalue <- shapiro.test(tool$Proportion)$p.value
      # }else{
      #   pvalue <- ad.test(tool$Proportion)$p.value ## Anderson-Darling检
      # }
      # if(pvalue > 0.05){  ## 符合正太分布
      #   method = "anova"
      # }else{
      #   method = "kruskal.test"
      # }
      
      method = "wilcox.test"
      if (i == "AO"){
        draw_boxplot_v1(tool,group_color = group_color,group_comp = group_comp,ncol = 6,filename = "ao_ssgsea_boxplot_v1.pdf",
                        ht=12,wd = 15)
        draw_boxplot_v2(tool,group_color = group_color,method = method,filename = "ao_ssgsea_boxplot_v2.pdf",
                        ht=3.5,wd = 12,font_size = font_size)
        
      }else if(i == "Danaher"){
        draw_boxplot_v1(tool,group_color = group_color,group_comp = group_comp,ncol = 4,filename = "danaher_ssgsea_boxplot_v1.pdf",
                        ht=9,wd = 12)
        draw_boxplot_v2(tool,group_color = group_color,method = method,filename = "danaher_ssgsea_boxplot_v2.pdf",
                        ht=3.5,wd = 6,font_size = font_size)
        
        
      }else if(i == "Immune29"){
        draw_boxplot_v1(tool,group_color = group_color,group_comp = group_comp,ncol = 6,filename = "immune29_ssgsea_boxplot_v1.pdf",
                        ht=14,wd = 16)
        draw_boxplot_v2(tool,group_color = group_color,method = method,filename = "immune29_ssgsea_boxplot_v2.pdf",
                        ht=3.5,wd = 12,font_size = font_size)
        
      }else if(i == "Immune28"){
        draw_boxplot_v1(tool,group_color = group_color,group_comp = group_comp,ncol = 6,filename = "immune28_ssgsea_boxplot_v1.pdf",
                        ht=14,wd = 16)
        draw_boxplot_v2(tool,group_color = group_color,method = method,filename = "immune28_ssgsea_boxplot_v2.pdf",
                        ht=3.5,wd = 12,font_size = font_size)
        
      }else if(i == "IFN"){
        draw_boxplot_v1(tool,group_color = group_color,group_comp = group_comp,ncol = 3,filename = "IFN_ssgsea_boxplot_v1.pdf",
                        ht=8,wd = 9)
        draw_boxplot_v2(tool,group_color = group_color,method = method,filename = "IFN_ssgsea_boxplot_v2.pdf",
                        ht=3.5,wd = 6,font_size = font_size)
        
        
      }else if(i == "xCell"){
        draw_boxplot_v1(tool,group_color = group_color,group_comp = group_comp,ncol = 6,filename = "xCell_ssgsea_boxplot_v1.pdf",
                        ht=18,wd = 20)
        draw_boxplot_v2(tool,group_color = group_color,method = method,filename = "xCell_ssgsea_boxplot_v2.pdf",
                        ht=3.5,wd = 16,font_size = font_size)
        
        
      }else if(i == "CAFs"){
        draw_boxplot_v1(tool,group_color = group_color,group_comp = group_comp,ncol = 3,filename = "CAFs_ssgsea_boxplot_v1.pdf",
                        ht=5,wd = 9)
        draw_boxplot_v2(tool,group_color = group_color,method = method,filename = "CAFs_ssgsea_boxplot_v2.pdf",
                        ht=3.5,wd = 5,font_size = font_size)
        
        
      }
    }
  }
}


draw_heatmap <- function(all_dat_file,groupA,groupB,pval=0.05,low=-1.5,high=0.8,width=5,height=4,
                         filename="all_sig_pheatmap.pdf"){
  all_dat <- read.table(all_dat_file,header = T,sep = "\t")
  tool_cal_all <- data.frame(A=NA,B=NA,Origin=NA)
  for(i in unique(all_dat$Origin)){
    tool <- all_dat[all_dat$Origin == i,]
    tool_cal <- tool %>% group_by(Cell_type,group) %>% 
      dplyr::summarise(values = list(Proportion)) %>%
      group_by(Cell_type) %>%
      dplyr::summarise(A=mean(values[[1]]),B=mean(values[[2]]),
                       pvalue=wilcox.test(values[[1]],values[[2]])$p.value)
    tool_cal <- as.data.frame(tool_cal)
    #tool_cal$sig <- ifelse(tool_cal$pvalue<0.001,"(***)",ifelse(tool_cal$pvalue<0.01,"(**)",ifelse(
    #  tool_cal$pvalue<0.05,"(*)","")))
    tool_cal$pvalue.new <- ifelse(tool_cal$pvalue < 0.001,"<0.001",round(tool_cal$pvalue,3))
    tool_cal$Cell_type_sig <- paste0(tool_cal$Cell_type," (",tool_cal$pvalue.new,")")
    #tool_cal$Cell_type_sig <- paste0(tool_cal$Cell_type,tool_cal$sig)
    tool_cal <- arrange(tool_cal,A)
    rownames(tool_cal) <- tool_cal$Cell_type_sig
    tool_cal$Origin <- i
    tool_cal_select <- tool_cal[tool_cal$pvalue < pval,]
    tool_cal_select <- tool_cal_select[,c("A","B","Origin")]
    tool_cal_select <- na.omit(tool_cal_select)
    tool_cal_all <- rbind(tool_cal_all,tool_cal_select)
    
  }
  tool_cal_all <- tool_cal_all[!is.na(tool_cal_all$Origin),]
  if(nrow(tool_cal_all) >0){
  if("Mcpcounter" %in% unique(tool_cal_all$Origin)){  ## mcpcounter值与其他方法得到的值相差较大
    mcps <- tool_cal_all[tool_cal_all$Origin == "Mcpcounter",]
    if("xCell" %in% unique(tool_cal_all$Origin)){
      xcell <- tool_cal_all[tool_cal_all$Origin == "xCell",]
      nonmcps <- tool_cal_all[tool_cal_all$Origin != "Mcpcounter" & tool_cal_all$Origin != "xCell",]
      a_sd <- sd(c(nonmcps$A,nonmcps$B))*sqrt((length(c(nonmcps$A,nonmcps$B))-1)/(length(c(nonmcps$A,nonmcps$B))))
      b_sd <- sd(c(mcps$A,mcps$B))*sqrt((length(c(mcps$A,mcps$B))-1)/(length(c(mcps$A,mcps$B))))
      c_sd <- sd(c(xcell$A,xcell$B))*sqrt((length(c(xcell$A,xcell$B))-1)/(length(c(xcell$A,xcell$B))))
      a_mean <- mean(c(nonmcps$A,nonmcps$B))
      b_mean <- mean(c(mcps$A,mcps$B))
      c_mean <- mean(c(xcell$A,xcell$B))
      nonmcps$A  <- (nonmcps$A - a_mean) /a_sd
      nonmcps$B <- (nonmcps$B - a_mean) /a_sd
      mcps$A <- (mcps$A - b_mean) /b_sd
      mcps$B <- (mcps$B -b_mean) / b_sd
      xcell$A <- (xcell$A - c_mean) /c_sd
      xcell$B <- (xcell$B -c_mean) / c_sd
      all_select <- rbind(nonmcps,mcps,xcell)
    }else{
      nonmcps <- tool_cal_all[tool_cal_all$Origin != "Mcpcounter",]
      a_sd <- sd(c(nonmcps$A,nonmcps$B))*sqrt((length(c(nonmcps$A,nonmcps$B))-1)/(length(c(nonmcps$A,nonmcps$B))))
      b_sd <- sd(c(mcps$A,mcps$B))*sqrt((length(c(mcps$A,mcps$B))-1)/(length(c(mcps$A,mcps$B))))
      a_mean <- mean(c(nonmcps$A,nonmcps$B))
      b_mean <- mean(c(mcps$A,mcps$B))
      nonmcps$A  <- (nonmcps$A - a_mean) /a_sd
      nonmcps$B <- (nonmcps$B - a_mean) /a_sd
      mcps$A <- (mcps$A - b_mean) /b_sd
      mcps$B <- (mcps$B -b_mean) / b_sd
      all_select <- rbind(nonmcps,mcps)
    }
  }else{
    all_select <- tool_cal_all
  }
    colnames(all_select) <- c(groupA,groupB,"Origin")
    print(all_select)
    ## draw
    library(circlize)
    library(ComplexHeatmap)
    
    data_for_barplot <- as.matrix(all_select[,c(1,2)])
    data_for_barplot[,1] <- as.numeric(data_for_barplot[,1])
    data_for_barplot[,2] <- as.numeric(data_for_barplot[,2])
    
    if(low == 0 & high == 0){
      pl <- Heatmap(as.matrix(all_select[,c(1,2)]),
                    name = "NES",
                    #col =  colorRampPalette(c("#00FFFF","#000000","#FD0000"))(100),
                    col =  colorRampPalette(c("#00FFFF","white","#FD0000"))(100),
                    cluster_rows = FALSE,
                    cluster_columns = FALSE,
                    column_names_side = "top",
                    column_names_rot = 45,
                    column_names_gp = gpar(fontsize = 10,fontface="bold"),
                    row_names_gp = gpar(fontsize = 10,fontface="bold"),
                    row_names_side = "left",
                    #width = ncol(all_select)*unit(5, "mm"),
                    #height = nrow(all_select)*unit(8, "mm"),
                    row_split = all_select$Origin,
                    row_title = NULL,
                    right_annotation = rowAnnotation(Origin=all_select$Origin,
                                                     # foo=anno_barplot(data_for_barplot,
                                                     #                  gp = gpar(fill = 2:3, col = 2:3)),
                                                     col=list(Origin=c("Immune29"="#9ca8b8",
                                                                       "Immune28"="#e0cdcf",
                                                                       "Cibersort"="#eee5f8",
                                                                       "Danaher"="#96a48b",
                                                                       "IFN" = "lightblue",
                                                                       "AO" = "pink",
                                                                       "Mcpcounter"="#faead3",
                                                                       "CAFs" = "#CBC3E3", # lightpurple
                                                                       "xCell" = "grey")),
                                                     show_legend = TRUE,
                                                     show_annotation_name=FALSE)
                    #right_annotation = rowAnnotation(foo = anno_barplot(as.matrix(tissue_select_new[,c(2,3)]), width = unit(4, "cm")))
                    #top_annotation = columnAnnotation(Tissue=rep("Tissue",nrow(tissue_select_new)),col=list(Tissue=c("Tissue"="#e0cdcf"))
      )
      pdf(filename,width = width,height = height,family = "ArialMT")
      print(pl)
      dev.off()
    }else{
    #col_fun = colorRamp2(c(-1, 0,1), c("#7A70B5", "#fdf9ee", "#FCBB44"))
    #col_fun = colorRamp2(c(low, 0, high), c("#00FFFF","#000000","#FD0000")) ## colorRamp2 小于-2 用#00FFFF，大于2 用#FD0000
      col_fun = colorRamp2(c(low, 0, high), c("#00FFFF","white","#FD0000"))
      pl <- Heatmap(as.matrix(all_select[,c(1,2)]),
                  name = "NES",
                  col = col_fun,
                  cluster_rows = FALSE,
                  cluster_columns = FALSE,
                  column_names_side = "top",
                  column_names_rot = 45,
                  column_names_gp = gpar(fontsize = 10,fontface="bold"),
                  row_names_gp = gpar(fontsize = 10,fontface="bold"),
                  row_names_side = "left",
                  #width = ncol(all_select)*unit(5, "mm"),
                  #height = nrow(all_select)*unit(8, "mm"),
                  row_split = all_select$Origin,
                  row_title = NULL,
                  right_annotation = rowAnnotation(Origin=all_select$Origin,
                                                  # foo=anno_barplot(data_for_barplot,
                                                  #                  gp = gpar(fill = 2:3, col = 2:3)),
                                                   col=list(Origin=c("Immune29"="#9ca8b8",
                                                                     "Immune28"="#e0cdcf",
                                                                     "Cibersort"="#eee5f8",
                                                                     "Danaher"="#96a48b",
                                                                     "IFN" = "lightblue",
                                                                     "AO" = "pink",
                                                                     "Mcpcounter"="#faead3",
                                                                     "CAFs" = "#CBC3E3", # lightpurple
                                                                     "xCell" = "grey")),
                    
                                                   show_legend = TRUE,
                                                   show_annotation_name=FALSE)
                  #right_annotation = rowAnnotation(foo = anno_barplot(as.matrix(tissue_select_new[,c(2,3)]), width = unit(4, "cm")))
                  #top_annotation = columnAnnotation(Tissue=rep("Tissue",nrow(tissue_select_new)),col=list(Tissue=c("Tissue"="#e0cdcf"))
    )
    pdf(filename,width = width,height = height,family = "ArialMT")
    print(pl)
    dev.off()
    }
  }else{
    print("No sigs")
  }
}


draw_heatmap_v2 <- function(all_dat_file,group_file,groupA,groupB,low=-2,high=2,pval=0.05,width=5,height=4,
                         filename="all_sig_pheatmap.pdf"){
  groups <- read.table(group_file,header = T,sep = "\t")
  all_dat <- read.table(all_dat_file,header = T,sep = "\t")
  tool_cal_all <- data.frame(A=NA,B=NA,Origin=NA,Cell_type=NA)
  for(i in unique(all_dat$Origin)){
    tool <- all_dat[all_dat$Origin == i,]
    tool_cal <- tool %>% group_by(Cell_type,group) %>% 
      dplyr::summarise(values = list(Proportion)) %>%
      group_by(Cell_type) %>%
      dplyr::summarise(A=mean(values[[1]]),B=mean(values[[2]]),
                       pvalue=wilcox.test(values[[1]],values[[2]])$p.value)
    tool_cal <- as.data.frame(tool_cal)
    #tool_cal$sig <- ifelse(tool_cal$pvalue<0.001,"(***)",ifelse(tool_cal$pvalue<0.01,"(**)",ifelse(
    #  tool_cal$pvalue<0.05,"(*)","")))
    tool_cal$pvalue.new <- ifelse(tool_cal$pvalue < 0.001,"<0.001",round(tool_cal$pvalue,3))
    tool_cal$Cell_type_sig <- paste0(tool_cal$Cell_type," (",tool_cal$pvalue.new,")")
    #tool_cal$Cell_type_sig <- paste0(tool_cal$Cell_type,tool_cal$sig)
    tool_cal <- arrange(tool_cal,A)
    rownames(tool_cal) <- tool_cal$Cell_type_sig
    tool_cal$Origin <- i
    tool_cal_select <- tool_cal[tool_cal$pvalue < pval,]
    #tool_cal_select <- tool_cal_select[,c("A","B","Origin")]
    tool_cal_select <- tool_cal_select[,c("A","B","Origin","Cell_type")]
    tool_cal_select <- na.omit(tool_cal_select)
    tool_cal_all <- rbind(tool_cal_all,tool_cal_select)
    
  }
  if(nrow(tool_cal_all) >0){
    tool_cal_all <- tool_cal_all[!is.na(tool_cal_all$Origin),]
    tool_cal_all$Sig <-  paste0(tool_cal_all$Origin,"_",tool_cal_all$Cell_type)
  ### select celltypes from matrix
    all_dat$Sig <- paste0(all_dat$Origin,"_",all_dat$Cell_type)
    all_dat_trans <- all_dat[,c("Sig","Sample","Proportion")] %>% tidyr::pivot_wider(names_from = Sample,
                                                                                   values_from = Proportion) %>% as.data.frame()
    rownames(all_dat_trans) <- all_dat_trans$Sig
    all_dat_trans <- all_dat_trans[,-1]
    all_dat_trans_select <- all_dat_trans[tool_cal_all$Sig,groups$sample]
    rownames(all_dat_trans_select) <- rownames(tool_cal_all)
  
    all_dat_trans_select_scale <- t(scale(t(all_dat_trans_select)))
  
    group1 <- groups[groups$group == groupA,]  ## control
    group2 <- groups[groups$group == groupB,]  ## case
  
    Delta <- rowMeans(all_dat_trans_select_scale[,group2$sample]) - rowMeans(all_dat_trans_select_scale[,group1$sample])
    Delta <- as.data.frame(Delta)
  
    ## draw
    library(circlize)
    library(ComplexHeatmap)
    group_cols <- c("#7570B3","#D95F02")
    names(group_cols) <- unique(groups$group)
    
    top.anno = HeatmapAnnotation(
      Group = groups$group,
      simple_anno_size = unit(0.6, "cm"),
      annotation_name_gp = gpar(fontsize = 8,fontface = "bold"),
      col = list(
        Group = group_cols
        ),
        annotation_legend_param = list(labels_gp = gpar(fontsize = 8),
                                       title_gp = gpar(fontsize = 8),
                                       #legend_gp = gpar(fill = 1:6),
                                       grid_width = unit(0.3, "cm"),
                                       grid_height = unit(0.5, "cm")
        )
    )
    
    row.anno1 <- rowAnnotation(`Delta` = anno_barplot(round(Delta$Delta,2),
                                                      which = c("row"),border = T,
                                                      gp = gpar(fill = ifelse(round(Delta$Delta,2) >0,"orange","blue")),
                                                      #ylim =,
                                                      width = unit(1.0, "cm")
    ),
    annotation_name_rot = 0)
    
    library(scales)
    n_path <- length(unique(tool_cal_all$Origin))
    path_col <- ggsci::pal_lancet(alpha = 0.7)(n_path)
    names(path_col) <- unique(tool_cal_all$Origin)
    rowsplits <- tool_cal_all$Origin
    rowsplits <- factor(rowsplits,levels = unique(rowsplits))
    #row_colors <- path_col[as.character(rowsplits)]
    #group_index <- as.integer(rowsplits)
    
    #old_pal <- palette()
    #palette(ggsci::pal_lancet(alpha = 0.7)(n_path))
    
    row.anno <- rowAnnotation(Origin=rowsplits,
                             col = list(Origin = path_col),
                             show_annotation_name=FALSE,
                             show_legend = FALSE,
                             #width = unit(5, "cm")
                             simple_anno_size = unit(3, "cm")
    )
    # row.anno <- rowAnnotation(Origin = anno_text(as.character(rowsplits), 
    #                                              location = 0.5,
    #                                              just = "center",
    #                                           gp = gpar(fill = group_index,col = "white",border = NULL)
    #                                           ))
    library(circlize)
    if(low ==0 | high == 0){
      bk1 <- seq(min(all_dat_trans_select_scale),0,by=0.01)
      bk2 = seq(0,max(all_dat_trans_select_scale),by=0.01)  
      col_fun = c(colorRampPalette(colors = c("#0e7f84","white"))(length(bk1)),
                  colorRampPalette(colors = c("white","#EE372C"))(length(bk2))
      )
      
    }else{
    col_fun = colorRamp2(c(low, 0, high), c("#2873A0","white","#E08130"))                          
    }
    #pdf("1.HER2_high_TME_heatmap.pdf",width = 15,height = 9,family = "ArialMT")
    pdf(filename,width = width,height = height,family = "ArialMT")
    print(Heatmap(all_dat_trans_select_scale,
                  #col = colorRampPalette(c("#293890", "white", "#BF1D2D"))(60),
                  #col =  colorRampPalette(c("#2A4C65","black","gold"))(80),
                  #col =  colorRampPalette(c("#00FFFF","white","#FD0000"))(100),
                  col = col_fun,
                  right_annotation = row.anno1,
                  left_annotation = row.anno,
                  top_annotation = top.anno,
                  cluster_rows = F,
                  cluster_columns = F,
                  column_dend_height = unit(10, "mm"),
                  row_dend_width = unit(10, "mm"), 
                  show_column_names = F,
                  column_labels = groups$group,
                  show_row_names = T,
                  row_split = rowsplits,
                  row_gap = unit(2,"mm"),
                  column_gap = unit(2,"mm"),
                  # column_km = 6,
                  # row_km = cl,
                  row_names_gp = gpar(fontsize = 8),
                  row_title_gp = gpar(fontsize = 8),
                  row_names_side = "right",
                  column_title_gp = gpar(fontsize = 8),
                  #column_split = groups$,
                  column_split = groups$group,
                  #row_split = term$g2,
                  row_title_rot = 0,
                  heatmap_legend_param = list(
                    title = "GSVA(scaled)",title_gp = gpar(fontsize = 8),
                    labels_gp = gpar(fontsize = 8)),
                  #row_names_max_width = max_text_width(
                  #  rownames(re), 
                  #  gp = gpar(fontsize = 12)
                  #)
    ))
    dev.off()
    
    
    
  }else{
    print("No sigs")
  }
}

library(ComplexHeatmap)
library(ggplot2)
library(ConsensusClusterPlus)
library(circlize)

draw_29_heatmap <- function(gs.exp,groups,cluster="kmeans",method="wilcox"){
  ## gs.exp row celltype column sample
  ## group 1st column: sample 2ed column: group
  if(cluster == "kmeans"){
  set.seed(1)
  c = kmeans(t(gs.exp), centers = 4,iter.max = 1000)
  c2 = c$cluster
  table(c2)
  identical(groups[,1], names(c2))   ## 1: sample 2: group
  groups$cluster = paste0("cluster.", c2)
  }else if(cluster == "consensus"){
    mat <- t(gs.exp)
    results <-  ConsensusClusterPlus(mat,  #mat: matrix: row: samples column: features
                                     maxK = 4,
                                     reps = 1000,              # 抽样次数(一般1000或更多)
                                     pItem = 0.8,              # 抽样比例
                                     pFeature = 1,
                                     clusterAlg = 'hc',       # 聚类方法:'hc'：层次聚类，'pam'：围绕中心点分区，'km'：K-means
                                     distance = 'pearson',       # 距离计算方法：'pearson', 'spearman', 'euclidean', 'binary', 'maximum', 'canberra', 'minkowski" or custom distance function.
                                     title="", # 结果保存路径
                                     innerLinkage="complete",  # 这里不建议使用默认的方法"average"
                                     plot=NULL, # png;pdf
                                     writeTable=FALSE)               # 结果保存形式
    
  }
  term = openxlsx::read.xlsx("/Users/stl/Documents/AmoyDX/scripts/Fge_pheatmap/Fge_genelist_group(1).xlsx")
  
  #cols <- c("#4DBBD5","#E64B35","#168675","#B88640")
  cols <- c("#4DBBD5","#E64B35")
  names(cols) <- unique(groups[,2])
  top.anno = HeatmapAnnotation(
    #treatment = clinical$新辅助化疗方案,
    Group = groups[,2],
    simple_anno_size = unit(0.6, "cm"),
    annotation_name_gp = gpar(fontsize = 12,fontface = "bold"),
    col = list(
      Group =  cols),
    annotation_legend_param = list(labels_gp = gpar(fontsize = 15),
                                   title_gp = gpar(fontsize = 15),
                                   #legend_gp = gpar(fill = 1:6),
                                   grid_width = unit(0.5, "cm"),
                                   grid_height = unit(0.5, "cm")
    )
    
  )
  rownames(gs.exp) <- gsub(" ", ".", rownames(gs.exp))
  rownames(gs.exp) <- gsub("-", ".", rownames(gs.exp))
  
  term$Sign = gsub(" ", ".", term$Sign)
  term$Sign = gsub("-", ".", term$Sign)
  gs.exp <- gs.exp[term$Sign,]
  
  identical(term$Sign, rownames(gs.exp))
  
  identical(groups[,1],colnames(gs.exp))
  
  #bk1 <- seq(min(gs.exp),0,by=0.01)
  #bk2 = seq(0,max(gs.exp),by=0.01)
  
  gs.exp.scale = t(scale(t(gs.exp)))
  
  # Define breakpoints and colors for gradient
  data_sd <- sd(min(gs.exp.scale):max(gs.exp.scale))
  breaks <- c(ceiling(min(gs.exp.scale)+data_sd),ceiling(max(gs.exp.scale)-data_sd))
  colors <- c("#002962","gold")
  col_fun <- colorRamp2(breaks, colors)
  
  p <- ComplexHeatmap::Heatmap(gs.exp.scale,
                               top_annotation = top.anno,
                               cluster_rows = F,
                               cluster_columns = F,
                               column_dend_height = unit(30, "mm"),
                               row_dend_width = unit(50, "mm"), 
                               show_column_names = F,
                               row_gap = unit(5,"mm"),
                               column_gap = unit(5,"mm"),
                               #col = c(colorRampPalette(colors = c("#0e7f84","white"))(length(bk1)),
                               #        colorRampPalette(colors = c("white","#EE372C"))(length(bk2))
                               #),
                               #col = colorRampPalette(colors = c("#002962","gold"))(length(bk1)+length(bk2)),
                               col = col_fun,
                               # column_km = 6,
                               # row_km = cl,
                               row_names_gp = gpar(fontsize = 20),
                               row_title_gp = gpar(fontsize = 20),
                               column_title_gp = gpar(fontsize = 30),
                               column_split = groups[,2],
                               row_split = term$g2,
                               row_title_rot = 0,
                               heatmap_legend_param = list(
                                 title = "Zscore",title_gp = gpar(fontsize = 15),
                                 labels_gp = gpar(fontsize = 15)),
                               row_names_max_width = max_text_width(
                                 rownames(gs.exp), 
                                 gp = gpar(fontsize = 12)
                               )
  )
  pdf("heatmap.fge.pdf",width = 17,height = 10,family = "ArialMT")
  print(p)
  dev.off()
  
  
  save(groups,cols,gs.exp,gs.exp.scale,col_fun,term,file = "draw.fge.rdata")
  
  term$g2 = factor(term$g2, levels = c("Malignant Cell Properties",   
                                      "Pro-tumor Microenvirment",
                                      "Anti-tumor Microenvirment",
                                      "Angiogenesis Fibrosis"
  ))
  term <- term[order(match(term$g2,c("Malignant Cell Properties",   
                                     "Pro-tumor Microenvirment",
                                     "Anti-tumor Microenvirment",
                                     "Angiogenesis Fibrosis"))),]
  
  gs.exp = t(gs.exp)
  gs.exp = as.data.frame(gs.exp)
  print(identical(rownames(gs.exp), groups[,1]))
  gs.exp$group = groups[,2]
  gs.exp <- gs.exp[,c(term$Sign,"group")]
  
  gene = names(gs.exp)[-ncol(gs.exp)]
  
  ###### wilcox
  if(method == "wilcox"){
    pval = sapply(gene, function(x){
      #x = "EMT.signature"
      #dat = kruskal.test(gs.exp[,x], gs.exp$group)
      #return(dat$p.value)
      formula = as.formula(paste0(x,"~ group"))
      dat = wilcox.test(formula = formula, data = gs.exp)
      
      return(dat$p.value)
      
    })}else if (method == "t.test"){
      ###### kruskal
      pval = sapply(gene, function(x){
        #x = "EMT.signature"
        #dat = kruskal.test(gs.exp[,x], gs.exp$group)
        #return(dat$p.value)
        formula = as.formula(paste0(x,"~ group"))
        dat = t.test(formula = formula, data = gs.exp)
        return(dat$p.value)
      })}
  #print(pval)
  pval.sig = round(pval,3)
  #print(pval.sig)
  gs.exp.mean = aggregate(.~group,gs.exp, mean)
  
  
  rownames(gs.exp.mean) = gs.exp.mean$group
  gs.exp.mean = gs.exp.mean[,-1]
  gs.exp.mean = t(gs.exp.mean)
  #library(GenAnalysis)
  #gs.exp.mean = scale_mat(gs.exp.mean, scale = "row")
 # bk1 <- seq(min(gs.exp.mean),-0.1,by=0.01)
  bk1 <- seq(min(gs.exp.mean),0,by=0.01)
  bk2 = seq(0,max(gs.exp.mean),by=0.01)
  
  
  
  #print(gs.exp.mean)
  #print(term)
  rownames(gs.exp.mean) = ifelse(as.numeric(pval.sig) <= 1 & as.numeric(pval.sig) > 0, 
                                 paste0(rownames(gs.exp.mean), "(",pval.sig,")"),
                                 ifelse(as.numeric(pval.sig) == 0,paste0(rownames(gs.exp.mean),"(<0.001)"),rownames(gs.exp.mean))
  )
  #print(gs.exp.mean)
  row.anno = rowAnnotation(Subgroup = term$g2,
                           col = list(Subgroup = c(
                             "Angiogenesis Fibrosis" = "#bc3c29",
                             "Pro-tumor Microenvirment" = "#20854e",
                             "Anti-tumor Microenvirment" = "#0072b5",
                             "Malignant Cell Properties" = "#e18727"
                           )
                           ),
                           show_annotation_name=FALSE
  )
  
  gs.exp.scale = t(scale(t(gs.exp.mean)))
  rownames(gs.exp.scale) <- factor(rownames(gs.exp.scale),levels = rownames(gs.exp.scale))
 # print(rownames(gs.exp.scale))
  p2 <- Heatmap(gs.exp.scale, cluster_rows = F,cluster_columns = F,
                show_column_names = T,
                row_names_gp = gpar(fontsize = 8,col = c("#e18727","#20854e","#0072b5","#bc3c29")),
                #row_names_gp = gpar(fontsize = 8,col = c("#0072b5","#20854e","#e18727","#bc3c29")),
                #row_names_gp = gpar(fontsize = 8,col = c("#bc3c29","#20854e","#0072b5","#e18727")),
                #row_names_gp = gpar(fontsize = 8,col = c("#e18727","#bc3c29","#0072b5","#20854e")),
                #row_names_gp = gpar(fontsize = 8,col = c("#bc3c29","#e18727","#20854e","#0072b5")),
                show_row_dend = F,border = "black",
                column_names_gp = gpar(fontsize = 10),
                column_names_rot = 45,column_names_centered = F,
                rect_gp = gpar(col = "white", lwd = 1),
                col = c(colorRampPalette(colors = c("#0e7f84","white"))(length(bk1)),
                        colorRampPalette(colors = c("white","#F34E77"))(length(bk2))
                ),
                row_names_side = "left",
                row_split = term$g2,
                row_order = rownames(gs.exp.scale),
                row_title_rot = 0,
                row_title = NULL,
                right_annotation = row.anno,
                #breaks = c(bk1, bk2),
                width = unit(4, "cm"),
                height = unit(10, "cm"),
                heatmap_legend_param = list(
                  title = "gsva.score\n(row_scaled)"#,
                  #legend_direction = c("horizontal")
                )
  )
  pdf("./29.immune.heatmap.mean.pdf",family = "ArialMT")
  print(p2)
  dev.off()
  
  
}
  

  
  