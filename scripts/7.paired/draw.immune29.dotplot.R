
library(ggsci)
library(ggplot2)
library(openxlsx)
library(tidyverse)

library(openxlsx)
library(ggplot2)
library(tibble)
library(dplyr)
library(tidyr)
library(ggpubr)

setwd("/Users/stl/Documents/AmoyDX/6.work_dir/10.lung_immune_adu/01.update_2/RNA/11.APOBEC_paired/")

clinical_dna <- read.xlsx("/Users/stl/Documents/AmoyDX/6.work_dir/10.lung_immune_adu/01.update_2/DNA/01.baseline_DNA_oncoplot/clinical_dna.xlsx")
clinical_dna_pair <- clinical_dna[!is.na(clinical_dna$surgery),]
clinical_dna_pair_out <- clinical_dna_pair[,c("patient_id","sample_id","surgery","APOBEC_group")]
clinical_dna_pair_out <- clinical_dna_pair_out[order(clinical_dna_pair_out$APOBEC_group),]
colnames(clinical_dna_pair_out) <- c("patient","baseline","surgery","group")

data1 <- clinical_dna_pair_out[,c(1,2,4)]
data2 <- clinical_dna_pair_out[,c(1,3,4)]
colnames(data1) <- c("patient","sample","group")
colnames(data2) <- c("patient","sample","group")
data <- rbind(data1,data2)
data$patient <- paste0("pt",data$patient)

mpr <- read.table("/Users/stl/Documents/AmoyDX/6.work_dir/19.nsclc_paired_and_pbmc/01.analysis/02.paired_analysis/MPR/all_TME_score.xls",header = T,sep = "\t",check.names = F)
nonmpr <- read.table("/Users/stl/Documents/AmoyDX/6.work_dir/19.nsclc_paired_and_pbmc/01.analysis/02.paired_analysis/NonMPR/all_TME_score.xls",header = T,sep = "\t",check.names = F)
num <- as.character(as.numeric(sapply(strsplit(nonmpr$ID,"pt"),'[[',2)) +30)
nonmpr$ID <- paste0(rep("pt",length(num)),num)
all.tmp <- rbind(mpr,nonmpr)
all.tmp$Origin_cell <- paste0(all.tmp$Origin,"_",all.tmp$Cell_type)
all.tmp_select <- all.tmp[all.tmp$Sample %in% data$sample,]

all.tmp_select$group3 <- NA
for(i in 1:nrow(all.tmp_select)){
  sample <- all.tmp_select[i,"Sample"]
  all.tmp_select[i,"group3"] <- data[data$sample == sample,]$group
}

for(i in 1:nrow(all.tmp_select)){
  sample <- all.tmp_select[i,"Sample"]
  all.tmp_select[i,"ID"] <- data[data$sample == sample,]$patient
}

all.tmp_select <- dplyr::arrange(all.tmp_select,Origin_cell,group,group3,ID)
all <- all.tmp_select

all$new_group <- paste0(all$group,"_",all$group3)
all$new_group <- factor(all$new_group,levels = c("baseline_neg","surgery_neg","baseline_pos","surgery_pos"))


############################## significant
all$sig <- paste0(all$Origin,".",all$Cell_type)

## wilcox test significant 
tool_select <- data.frame()
tool_select <- all[NA,]
tool_select <- tool_select[1,]

data_for_dot_plot <- data.frame(celltype=unique(all$sig),neg_delta=NA,pos_delta=NA,
                                neg_pvalue=NA,pos_pvalue=NA) 
rownames(data_for_dot_plot) <- data_for_dot_plot$celltype
for(i in unique(all$sig)){
  tool <- all[all$sig == i,]
  if(sum(tool$Proportion) != 0){
    mpr <- tool[tool$group3 == "neg",]
    nonmpr <- tool[tool$group3 == "pos",]
    res1 <- wilcox.test(Proportion ~ new_group, data = mpr, paired = TRUE)
    res2 <- wilcox.test(Proportion ~ new_group, data = nonmpr, paired = TRUE)
    pvalue1 <- res1$p.value
    pvalue2 <- res2$p.value
    pvalue1 <- ifelse(is.na(pvalue1),9,pvalue1)
    pvalue2 <- ifelse(is.na(pvalue2),9,pvalue2)
    delta1 <- median(mpr[mpr$new_group=="surgery_neg",]$Proportion - mpr[mpr$new_group=="baseline_neg",]$Proportion)
    delta2 <- median(nonmpr[nonmpr$new_group=="surgery_pos",]$Proportion - nonmpr[nonmpr$new_group=="baseline_pos",]$Proportion)
    if(pvalue1<0.05 & is.na(pvalue2)){
      tool_select <- rbind(tool_select,tool)
      data_for_dot_plot[i,] <- c(i,delta1,delta2,pvalue1,pvalue2)
    }else if(is.na(pvalue1) & pvalue2 < 0.05){
      tool_select <- rbind(tool_select,tool)
      data_for_dot_plot[i,] <- c(i,delta1,delta2,pvalue1,pvalue2)
    }else if(pvalue1<0.05 | pvalue2 <0.05){
      tool_select <- rbind(tool_select,tool)
      data_for_dot_plot[i,] <- c(i,delta1,delta2,pvalue1,pvalue2)
    }
  }
}

tool_select <- tool_select[!is.na(tool_select$Sample),]
tool_select <- tool_select[!grepl("immuneAI",tool_select$Origin),]

data_for_dot_plot <- data_for_dot_plot[!is.na(data_for_dot_plot$neg_delta),]

## immune29
immune29_class <- read.xlsx("/Users/stl/Documents/AmoyDX/6.work_dir/10.lung_immune_adu/01.update_2/RNA/11.APOBEC_paired/immune29_class.xlsx")

immune29_for_dot_plot <- data_for_dot_plot[grepl("Immune29",data_for_dot_plot$celltype),]
immune29_for_dot_plot_select <- immune29_for_dot_plot[c(-2,-6,-7,-10,-11,-15,-16),]

immune29_for_dot_plot_select$celltype <- gsub("Immune29.","",immune29_for_dot_plot_select$celltype)
rownames(immune29_for_dot_plot_select) <- immune29_for_dot_plot_select$celltype 
immune29_for_dot_plot_select <- immune29_for_dot_plot_select %>% dplyr::left_join(immune29_class,
                                                 by=c("celltype"="to"))

colnames(immune29_for_dot_plot_select)[6] <- "group"
immune29_for_dot_plot_select$neg_pvalue <- as.numeric(immune29_for_dot_plot_select$neg_pvalue)
immune29_for_dot_plot_select$pos_pvalue <- as.numeric(immune29_for_dot_plot_select$pos_pvalue)
immune29_for_dot_plot_select$neg_pvalue_log10 <- -log10(as.numeric(immune29_for_dot_plot_select$neg_pvalue))
immune29_for_dot_plot_select$pos_pvalue_log10 <- -log10(as.numeric(immune29_for_dot_plot_select$pos_pvalue))

immune29_for_dot_plot_select$neg_pvalue_sig = ifelse(immune29_for_dot_plot_select$neg_pvalue < 0.05, "<0.05", ">=0.05")
immune29_for_dot_plot_select$pos_pvalue_sig = ifelse(immune29_for_dot_plot_select$pos_pvalue < 0.05, "<0.05", ">=0.05")
immune29_for_dot_plot_select$celltype = factor(immune29_for_dot_plot_select$celltype, levels = immune29_for_dot_plot_select$celltype)
p1<- ggplot(immune29_for_dot_plot_select, aes(x = neg_delta, y= celltype)) +
  geom_col(width = 0.1, #调整宽度，使柱子变成一根细细的'棍子'
           aes(fill= group)) + 
  geom_point(aes(size = abs(neg_pvalue_log10), fill = group,
                 color = neg_pvalue_sig),shape = 21,stroke = 1.5) +
  scale_color_manual(values = c("<0.05" = "black", ">=0.05" = "white"))+
  scale_size_continuous(range = c(2, 5)) +
  scale_fill_nejm()

neg_immune29 <- immune29_for_dot_plot_select[,c(1,2,4,6,7,9)]
colnames(neg_immune29) <- c("celltype","delta","pvalue","group","pvalue_log10","sig")

p2<- ggplot(immune29_for_dot_plot_select, aes(x = pos_delta, y= celltype)) +
  geom_col(width = 0.1, #调整宽度，使柱子变成一根细细的'棍子'
           aes(fill= group)) + 
  geom_point(aes(size = abs(pos_pvalue_log10), fill = group,
                 color = pos_pvalue_sig),shape = 21,stroke = 1.5) +
  scale_color_manual(values = c("<0.05" = "black", ">=0.05" = "white"))+
  scale_size_continuous(range = c(2, 5)) +
  scale_fill_nejm()

pos_immune29 <- immune29_for_dot_plot_select[,c(1,3,5,6,8,10)]
colnames(pos_immune29) <- c("celltype","delta","pvalue","group","pvalue_log10","sig")
df.merge = list(neg = neg_immune29, pos = pos_immune29)
df.merge = plyr::ldply(df.merge,data.frame)

colnames(df.merge)[1] <- "apobec"
df.merge$delta <- as.numeric(df.merge$delta)
p.merge = ggplot(df.merge, aes(x = delta, y= celltype)) +
  geom_col(width = 0.1, #调整宽度，使柱子变成一根细细的'棍子'
           aes(fill= group)) + 
  geom_point(aes(size = abs(pvalue_log10), fill = group,
                 color = sig),shape = 21,stroke = 1.5) +
  scale_color_manual(values = c("<0.05" = "black", ">=0.05" = "grey"))+
  scale_size_continuous(range = c(2, 5)) +
  scale_fill_nejm(breaks = c("Malignant cells",
                             #"Angiogenesis Fibrosis",
                             #"Pro-tumor Microenvirment",
                             "Anti-tumor immunity")) +
  labs(x = "GSVA delta") +
  geom_vline(xintercept = 0, lty = "dashed")+
  theme(#panel.grid = element_blank(), 
    panel.background = element_blank(), 
    strip.text =element_text(size = 12,
                             color = 'black'),
    panel.border = element_rect(color = "black",fill = NA), 
    axis.line = element_line(color = 'black'),
    axis.text.y = element_text(
      color = 'black'),
    axis.title.x = element_text(size = 15,
                                color = 'black')
  ) +
  facet_grid(.~apobec)

p0 = ggplot(pos_immune29, aes(y=celltype))+
  geom_text(aes(x = 0, label = celltype,color = group), 
            hjust = 1, fontface = "bold",size = 3)+
  scale_color_nejm()+
  theme_void()+
  theme(legend.position = "none")+
  coord_cartesian(xlim = c(-0.5, -0.02))

library(patchwork)
library(ggpubr)
p.final = p0 + p.merge + rremove("ylab")+ rremove("y.text") +rremove("grid")#+ 
# plot_layout(design = layout)
ggsave(p.final, file = "fge.plot.pdf", width = 8,height = 6)

### 2025.6.23 new Fge plot

p.merge = ggplot(df.merge, aes(x = delta, y= celltype)) +
  geom_col(width = 0.1 #调整宽度，使柱子变成一根细细的'棍子'
           ) + 
  geom_point(aes(size = abs(pvalue_log10),fill=ifelse(delta >0,"blue",
                                                      "red"),
                 color = sig),shape = 21,stroke = 1.5) +
  scale_color_manual(values = c("<0.05" = "black", ">=0.05" = "grey"))+
  scale_size_continuous(range = c(2, 5)) +
  #scale_fill_nejm(breaks = c("Malignant cells",
  #                           #"Angiogenesis Fibrosis",
  #                           #"Pro-tumor Microenvirment",
  #                           "Anti-tumor immunity")) +
  labs(x = "GSVA delta") +
  geom_vline(xintercept = 0, lty = "dashed")+
  theme(#panel.grid = element_blank(), 
    panel.background = element_blank(), 
    strip.text =element_text(size = 12,
                             color = 'black'),
    panel.border = element_rect(color = "black",fill = NA), 
    axis.line = element_line(color = 'black'),
    axis.text.y = element_text(
      color = 'black'),
    axis.title.x = element_text(size = 15,
                                color = 'black')
  ) +
  facet_grid(.~apobec)

p0 = ggplot(pos_immune29, aes(y=celltype))+
  geom_text(aes(x = 0, label = celltype), 
            hjust = 1, fontface = "bold",size = 3)+
  #scale_color_nejm()+
  theme_void()+
  theme(legend.position = "none")+
  coord_cartesian(xlim = c(-0.5, -0.02))

library(patchwork)
library(ggpubr)
p.final = p0 + p.merge + rremove("ylab")+ rremove("y.text") +rremove("grid")#+ 
# plot_layout(design = layout)
ggsave(p.final, file = "fge.plot_v2.pdf", width = 8,height = 6,family="ArialMT")



## immune28
immune28_class <- read.xlsx("/Users/stl/Documents/AmoyDX/6.work_dir/10.lung_immune_adu/01.update_2/RNA/11.APOBEC_paired/immune28_class.xlsx")

immune28_for_dot_plot <- data_for_dot_plot[grepl("Immune28",data_for_dot_plot$celltype),]


