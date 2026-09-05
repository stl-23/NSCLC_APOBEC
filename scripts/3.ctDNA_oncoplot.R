rm(list = ls())
library(ComplexHeatmap)
library(openxlsx)
library(ggsci)
library(dplyr)
library(tidyverse)
source("./createOncoMatrix1.R")

group <- read.xlsx("./lung_PL0_and_RNA(3).xlsx",sheet = 4)
group$response_new <- ifelse(group$response == "Non-MPR","Non.MPR","MPR")
group$response_new <- factor(group$response_new,levels = c("MPR","Non.MPR"))
group <- group[!is.na(group$PL0) & group$treatment != "Chemotherapy",]

tmb <- read.xlsx("./ctDNA/01.PL0_oncoplot/TMB.xlsx")

group <- group %>% dplyr::left_join(tmb,
                                    by=c("PL0"="sample")) 

library(stringr)
sig = read.xlsx("./ctDNA/05.sigma/lite.signature.xlsx")
head(sig)
#sig$tumor = str_extract(sig$tumor, "Q.vcf")
sig$tumor <- gsub(".vcf","",sapply(strsplit(sig$tumor,"/"),tail,1))
#View(sig)
sigma.dat = sig[,c("tumor","Signature_clock_ml", "Signature_APOBEC_ml",
                   "Signature_28_ml", "Signature_18_ml","Signature_4_ml",
                   "Signature_3_ml", "Signature_16_ml")]
#colnames(sigma.dat) <- gsub("_ml","",colnames(sigma.dat))
group <- group %>% dplyr::left_join(sigma.dat,
                                    by=c("PL0"="tumor"))

group$PDL1 <- ifelse(group$`PDL1(1=0,2=1-49%,3=>50%)` == 0 |group$`PDL1(1=0,2=1-49%,3=>50%)` == 1 ,"PD-L1-","PD-L1+")
group$age <- as.numeric(group$年龄)
group$smoking <- as.numeric(group$`吸烟史（0：无；1：有）`)
group$type <- ifelse(group$type == "SQC","Squamous","Non.Squamous")
group$stage <- group$治疗前临床分期

#group$APOBEC_group <- ifelse(group$Signature_APOBEC_ml > quantile(group$Signature_APOBEC_ml,0.33,na.rm = T),"+","-")
#group$APOBEC_group <- factor(group$APOBEC_group,levels = c("+","-"))

library(pROC)
roc2 <- roc(response=group$response_new,predictor=group$Signature_APOBEC_ml,
            #smoothed = TRUE,
            #smooth = TRUE,
            levels=c("Non.MPR","MPR"),
            direction = "<") ## 手动设定 mplc < im
roc2$auc
pdf("ctDNA_APOBEC_ROC.new.pdf",width = 4,height = 4,family = "ArialMT")
pROC::plot.roc(roc2,
               legacy.axes=TRUE,  ## 1-sepecificity if TRUE
               print.auc=TRUE,   ## add AUC
               print.thres.col="#002962", col="#002962",
               print.thres="best", print.thres.best.method="youden")  ## add best cutoff
dev.off()

bestcutoff <- coords(roc2, "best")  ## 默认为youden
bestcutoff <- bestcutoff$threshold
# 0.05931694

group$APOBEC_group.bestcut <- ifelse(group$Signature_APOBEC_ml > bestcutoff,"pos","neg")
table(group$response_new,group$APOBEC_group.bestcut)

group.pre <- group
group <- group[!is.na(group$APOBEC_group.bestcut),]

### ### 2024.10.18 去除掉1例样本p20，与45例ctDNA文章35例有bTMB的数据对应
### barplot
group$APOBEC_group.bestcut <- ifelse(group$APOBEC_group.bestcut == "+","APOBEC+","APOBEC-")
group$APOBEC_group.bestcut <- factor(group$APOBEC_group.bestcut,levels = c("APOBEC+","APOBEC-"))

library(ggtext) ## for annotate richtext
library(ggplot2)
pvalue0 <- round(fisher.test(table(group$APOBEC_group.bestcut,group$response_new))$p.value,3)
response_table <- as.data.frame(table(group$APOBEC_group.bestcut,group$response_new))
names(response_table)[1:2] = c("Group","Response")
#response_table$Response <- gsub("non.response","NonMPR",response_table$Response)
#response_table$Response <- gsub("response","MPR",response_table$Response)
response_table$count = ifelse(response_table$Freq == 0, "", paste0("n=", response_table$Freq))
response_table$Response <- gsub("Non.MPR","NonMPR",response_table$Response)
response_table$Response <- factor(response_table$Response,levels = c("NonMPR","MPR"))
response_table$Group <- gsub("neg","APOBEC-",response_table$Group)
response_table$Group <- gsub("pos","APOBEC+",response_table$Group)
response_table$Group <- factor(response_table$Group,levels = c("APOBEC+","APOBEC-"))

p0 = ggplot(response_table, aes(x = Group,y = Freq, fill = Response)) + 
  geom_bar(stat = "identity", position = "fill",width = 0.6,color = "black") + 
  #coord_flip()+
  scale_y_continuous(labels = scales::percent,breaks = seq(0,1, 0.2)) +
  geom_text(aes(label = count), 
            colour = "black",
            position = position_fill(vjust = 0.5),size = 6) +
  #scale_y_continuous(limits = c(0,1)) + 
  #scale_fill_manual(values = c( PD = "#bc3c28",PR = "#0072b5",SD= "#a27e7e")) +
  #scale_fill_manual(values = c( NonMPR = "#bc3c28",MPR = "#0072b5")) +
  #scale_fill_manual(values = c( NonMPR ="#DD9489",MPR = "#5FB155")) +
  #scale_fill_manual(values = c( NonMPR = "#B696B6",MPR = "#F6C8A8")) +
  scale_fill_manual(values = c( NonMPR = "#2485C8",MPR = "#A27E7E")) +
  #scale_fill_lancet() +
  geom_segment(x = 1, xend = 2, y = 1.05, yend = 1.05)+
  geom_segment(x = 1, xend = 1, y = 1.03, yend = 1.05)+
  geom_segment(x = 2, xend = 2, y = 1.03, yend = 1.05)+
  annotate("text", x=1.5, y=1.1, parse = TRUE,label= paste0("italic(pval)==", pvalue0)) +
  #annotate("richtext", x=1.5, y=1.1,label= paste0("<i>pval</i>=", pvalue0),angle=-90,
  #         label.color=NA,fill=NA) +
  #annotate("richtext", x=2.5, y=0.5,label= paste0("<i>pval</i> = ", pvalue0),angle=0,
  #         label.color=NA,fill=NA) +
  labs(x = element_blank(), y = "Percent",title = "") +
  #theme_bw()+
  theme(axis.text.x = element_text(colour = "black",size = 10),
        axis.text.y = element_text(colour = "black",size = 10),
        axis.title.y = element_text(colour = "black",size = 15),
        plot.title = element_text(hjust = 0.5),
        axis.line = element_line(colour = "grey"),
        panel.grid = element_blank(),panel.background = element_blank(),
        panel.border = element_rect(fill = NA,color="grey", linetype="solid")
  )
p0
#print(p)
ggsave(p0, file = paste0("ctDNA_APOBEC_bar_v7", ".pdf"), width = 4, height = 4,family="ArialMT")


### oncoplot
mat <- read.xlsx("./1.all_PL0/ctDNA_snv_filter_c.xlsx",sheet=1)
mat[is.na(mat$`1KGEAS`),"1KGEAS"] <- 0
mat[mat$`1KGEAS` == "-","1KGEAS"] <- 0
mat[is.na(mat$GnomAD),"GnomAD"] <- 0
mat[mat$GnomAD == "-","GnomAD"] <- 0
mat$`1KGEAS` <- as.numeric(mat$`1KGEAS`)
mat$GnomAD <- as.numeric(mat$GnomAD)

mat = mat %>% filter((`1KGEAS` <= 0.01) &
                               (GnomAD <= 0.01))


mat = mat %>% filter(!Type %in% c("-","3'UTR", "5'UTR",
                                                "CorruptedStart", "Extension", 
                                                "Synonymous_Substitution","Intronic",
                                                "FlankingRegion5"))

mat <- mat %>% dplyr::mutate(Type=case_when(Type=="FrameShift_Deletion" ~ "Frame_Shift_Del",
                                            Type=="FrameShift_Duplication" ~ "Frame_Shift_Ins",
                                            Type=="FrameShift_Insertion" ~ "Frame_Shift_Ins",
                                            Type=="FrameShift_Substitution" ~ "Frameshift_INDEL",
                                            Type=="nonFrameShift_Deletion" ~ "In_Frame_Del",
                                            Type=="nonFrameShift_Insertion" ~ "In_Frame_Ins",
                                            Type=="nonFrameShift_Substitution" ~ "In_Frame_INDEL",
                                            Type=="nonFrameShift_Substitution" ~ "In_Frame_INDEL",
                                            Type=="Nonsense_Mutation" ~ "Nonsense_Mutation",
                                            Type=="nonSynonymous_Substitution" ~ "Missense_Mutation",
                                            Type=="Splicing" ~ "Splice_Site",
                                            TRUE ~ NA_character_
                                            ))


mat <- mat[mat$Sample %in% group$PL0,]

## add group
mat$Stage <- apply(mat, 1, function(x) {
  group$stage[match(x[1], group$PL0)]
})
mat$Treatment <- apply(mat, 1, function(x) {
  group$treatment[match(x[1], group$PL0)]
})
mat$Response <- apply(mat, 1, function(x) {
  group$response[match(x[1], group$PL0)]
})
mat$Subtype <- apply(mat, 1, function(x) {
  group$type[match(x[1], group$PL0)]
})


mat <- data.table::as.data.table(mat)


mut_matrix <- createOncoMatrix(mat,g=unique(mat$Gene))

diffsample <- setdiff(unique(group$PL0),unique(mat$Sample))
add_blank <- data.frame(matrix(nrow = length(diffsample),ncol = length(colnames(mat))))
add_blank$X1 <- diffsample
colnames(add_blank) <- colnames(mat)
mat <- rbind(mat,add_blank)

tmp1 <- data.frame(matrix(nrow = nrow(mut_matrix$oncoMatrix),ncol = length(diffsample)))
rownames(tmp1) <- rownames(mut_matrix$oncoMatrix)
colnames(tmp1) <- diffsample
tmp2 <- tmp1
tmp1[is.na(tmp1)] <- ""
tmp2[is.na(tmp2)] <- 0
mut_matrix$oncoMatrix <- cbind(mut_matrix$oncoMatrix,tmp1)
mut_matrix$numericMatrix <- cbind(mut_matrix$numericMatrix,tmp2)


response <- group[group$response_new == "MPR",]
nonresponse <- group[group$response_new == "Non.MPR",]

snv_matrix <- mut_matrix$oncoMatrix
snv_matrix_response <- snv_matrix[,colnames(snv_matrix) %in% response$PL0]
snv_matrix_nonresponse <- snv_matrix[,colnames(snv_matrix) %in% nonresponse$PL0]

response <- response[match(colnames(snv_matrix_response),response$PL0),]
nonresponse <- nonresponse[match(colnames(snv_matrix_nonresponse),nonresponse$PL0),]
group <- rbind(response,nonresponse)


response_top30 <- head(snv_matrix_response,n=20)
nonresponse_top30 <- head(snv_matrix_nonresponse,n=20)

# draw
col = c(RColorBrewer::brewer.pal(n = 12, name = 'Paired'))
names(col) = c(
  'Nonsense_Mutation',
  'Missense_Mutation',
  'Frame_Shift_Del',
  'Frame_Shift_Ins',
  'Frameshift_INDEL',
  'In_Frame_Del',
  'In_Frame_INDEL',
  'In_Frame_Ins',
  'Multi_Hit',
  'Splice_Site'
)
col[9] <- "black"  ## Multi_Hit 


#load("./predata.rda")
#library(randomcoloR)

library(RColorBrewer)

alter_fun <- list(
  background = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.8, "mm"), h-unit(0.5, "mm"),
              #gp = gpar(fill = "#DDDDDD", col = NA))
              gp = gpar(fill = "white", col = NA))
  },
  Frame_Shift_Del = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.8, "mm"), h-unit(0.5, "mm"), #高度可调整h*0.33
              gp = gpar(fill = col["Frame_Shift_Del"], col = NA))
  },
  Frame_Shift_Ins = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.8, "mm"), h-unit(0.5, "mm"),
              gp = gpar(fill = col["Frame_Shift_Ins"], col = NA))
  },
  Frameshift_INDEL = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.8, "mm"), h-unit(0.5, "mm"),
              gp = gpar(fill = col["Frameshift_INDEL"], col = NA))
  },
  Multi_Hit = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.8, "mm"), h-unit(0.5, "mm"),
              gp = gpar(fill = col["Multi_Hit"], col = NA))
  },
  In_Frame_Del  = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.8, "mm"), h-unit(0.5, "mm"),
              gp = gpar(fill = col["In_Frame_Del"], col = NA))
  },
  In_Frame_INDEL  = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.8, "mm"), h-unit(0.5, "mm"),
              gp = gpar(fill = col["In_Frame_INDEL"], col = NA))
  },
  In_Frame_Ins  = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"),
              gp = gpar(fill = col["In_Frame_Ins"], col = NA))
  },
  Missense_Mutation  = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.8, "mm"), h-unit(0.5, "mm"),
              gp = gpar(fill = col["Missense_Mutation"], col = NA))
  },
  Nonsense_Mutation  = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.8, "mm"), h-unit(0.5, "mm"),
              gp = gpar(fill = col["Nonsense_Mutation"], col = NA))
  },
  Splice_Site  = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.8, "mm"), h-unit(0.5, "mm"),
              gp = gpar(fill = col["Splice_Site"], col = NA))
  }
  
)

response_tmb <- tmb[tmb$sample %in% response$PL0,]
response_tmb <- response_tmb[order(match(response_tmb$sample,colnames(response_top30))),]
response_tmb$sample <- as.character(response_tmb$sample)

nonresponse_tmb <- tmb[tmb$sample %in% nonresponse$PL0,]
nonresponse_tmb <- nonresponse_tmb[order(match(nonresponse_tmb$sample,colnames(nonresponse_top30))),]
nonresponse_tmb$sample <- as.character(nonresponse_tmb$sample)

###设置注释信息
# library(randomcoloR)
# tmpcolor <- distinctColorPalette(33) ## 挑选样本对应的颜色
library(randomcoloR)
stage_color <- distinctColorPalette(4)
names(stage_color) = names(table(group$stage))
#treatment_color <- c("#AEB38C","#E7E359","#E35786")
#response_color <- c("#5FB155","#7BBCE1","#DD9489")
response_color <- c("#5FB155","#DD9489")
names(response_color) = c("MPR","Non.MPR")

#type_color <- distinctColorPalette(20)[c(8,12)]
type_color <- c("#E1C0A6","#863CE2")
names(type_color) = names(table(group$type))

pdl1_color <- c("lightblue","#E35786")
names(pdl1_color) = names(table(group$PDL1))

top.anno1 <- columnAnnotation(
  TMB = anno_barplot(as.numeric(response_tmb$TMB),ylim = c(0,40),  axis_param = list(at = seq(0,40, by = 20)),
                     border = F),
  #Treatment = clinical.merge$Treatment,
  Group = response$response_new,
  Subtype = response$type,
  Stage = response$stage,
  PDL1 = response$PDL1,
  #institution = factor(clinical.merge$institution,levels = c("ZH","AD")),
  foo = anno_empty(border = F,height = unit(c(0.05), "cm")), ###加入空注释，以调整注释与热图之间的间距
  simple_anno_size = unit(c(0.5), "cm"), ###每个注释的宽度
  col = list(#Treatment = treatment_color,
             Group = response_color,
             Subtype = type_color,
             Stage = stage_color,
             PDL1 = pdl1_color
  ),
  #axis_param = list(at = c(0,5,10,15,20)),
  show_annotation_name = T,
  annotation_name_side = "left", ###注释的标签在左还是右
  annotation_name_gp = gpar(fontsize = 10) ####标签的参数修改
  #gap = unit(c(1, 1,1,0), "mm") ###各个注释之间的间距
)

top.anno2 <- columnAnnotation(
  TMB = anno_barplot(as.numeric(nonresponse_tmb$TMB),ylim = c(0,40),  axis_param = list(at = seq(0,40, by = 20)),
                     border = F),
  #Treatment = clinical.merge$Treatment,
  Group = nonresponse$response_new,
  Subtype = nonresponse$type,
  Stage = nonresponse$stage,
  PDL1 =  nonresponse$PDL1,
  #institution = factor(clinical.merge$institution,levels = c("ZH","AD")),
  foo = anno_empty(border = F,height = unit(c(0.05), "cm")), ###加入空注释，以调整注释与热图之间的间距
  simple_anno_size = unit(c(0.5), "cm"), ###每个注释的宽度
  col = list(#Treatment = treatment_color,
    Group = response_color,
    Subtype = type_color,
    Stage = stage_color,
    PDL1 = pdl1_color
  ),
  #axis_param = list(at = c(0,5,10,15,20)),
  show_annotation_name = T,
  annotation_name_side = "left", ###注释的标签在左还是右
  annotation_name_gp = gpar(fontsize = 10) ####标签的参数修改
  #gap = unit(c(1, 1,1,0), "mm") ###各个注释之间的间距
)

column_order1 = colnames(response_top30)
column_order2 = colnames(nonresponse_top30)


response_sigs <- response[,c(4,40:46)]
rownames(response_sigs) <- response_sigs$PL0
response_sigs <- response_sigs[,-1]
colnames(response_sigs) <- gsub("_ml","",colnames(response_sigs))
response_sigs <- apply(response_sigs,2,function(x){x*100})
response_sigs[is.na(response_sigs)] <- 0
response_sigs <- response_sigs[,c(2,1,3:7)]

nonresponse_sigs <- nonresponse[,c(4,40:46)]
rownames(nonresponse_sigs) <- nonresponse_sigs$PL0
nonresponse_sigs <- nonresponse_sigs[,-1]
colnames(nonresponse_sigs) <- gsub("_ml","",colnames(nonresponse_sigs))
nonresponse_sigs <- apply(nonresponse_sigs,2,function(x){x*100})
nonresponse_sigs[is.na(nonresponse_sigs)] <- 0
nonresponse_sigs <- nonresponse_sigs[,c(2,1,3:7)]

bottom.anno1 <- HeatmapAnnotation(
  empty = anno_empty(border = FALSE,height = unit(0, "cm")),
  Signature = anno_barplot(response_sigs,
                           gp = gpar(fill = c("#a27e7e","#1E4C9C","#345D82",
                                              "#3371B3","#5795C7","#81B5D5",
                                              "#AED4E5"),alpha = 0.7
                           ),
                           axis_param = list(
                             gp = gpar(fontsize = 10)),
                           bar_width = 1),
  height = unit(6, "cm"),
  annotation_legend_param = list(labels_gp = gpar(fontsize = 10),
                                 title_gp = gpar(fontsize = 10),
                                 legend_gp = gpar(fill = 1:6),
                                 grid_width = unit(0.5, "cm"),
                                 grid_height = unit(0.5, "cm")
  )
)

bottom.anno2 <- HeatmapAnnotation(
  empty = anno_empty(border = FALSE,height = unit(0, "cm")),
  Signature = anno_barplot(nonresponse_sigs,
                           gp = gpar(fill = c("#a27e7e","#1E4C9C","#345D82",
                                              "#3371B3","#5795C7","#81B5D5",
                                              "#AED4E5"),alpha = 0.7
                           ),
                           axis_param = list(
                             gp = gpar(fontsize = 10)),
                           bar_width = 1),
  height = unit(6, "cm"),
  annotation_legend_param = list(labels_gp = gpar(fontsize = 10),
                                 title_gp = gpar(fontsize = 10),
                                 legend_gp = gpar(fill = 1:6),
                                 grid_width = unit(0.5, "cm"),
                                 grid_height = unit(0.5, "cm")
  )
)

#group = rep(c("group1", "group2"), each = 7)
hlist1 = oncoPrint(response_top30,
                  alter_fun = alter_fun,
                  top_annotation = top.anno1,
                  #heatmap_height = unit(15, "cm"),
                  bottom_annotation = bottom.anno1,
                  #heatmap_width = unit(1000, "cm"),
                  right_annotation = NULL,
                  #left_annotation = left_anno,
                  col = col,
                  #gap = unit(0.1, "cm"),
                  get_type = function(x) strsplit(x, ";")[[1]],###按照分号区分重叠的类型
                  remove_empty_columns = F,
                  #column_order = col.order,
                  row_order = rownames(response_top30), ###行的顺序
                  remove_empty_rows = F, ###是否删除空的行
                  row_names_side = "left", #基因在左
                  pct_side = "right",###比例在右
                  show_pct = T, ##是否显示比例
                  pct_digits = 0, ##保留几位小数点
                  #row_labels = paste0(rownames(dat),"(10%)"),
                  #column_split = as.factor(clinical.merge$Response), ###列的分开
                  column_gap = unit(1, "mm"), ###分开的间隙
                  column_order = column_order1,
                  border=T,
                  show_column_names = F,
                  row_gap = unit(0, "mm"),
                  #row_title = NULL, ###行标签,设置为NULL即为空白
                  pct_gp = gpar(fontsize = 10,fontface="bold"), ##百分比字体
                  row_names_gp = gpar(fontsize = 10,fontface="bold"), ####基因名称字体
                  heatmap_legend_param = list(
                    #at = c(-2, 0, 2),
                    #labels = c("low", "zero", "high"),
                    title = "Mutations"####,图注标题
                    #legend_height = unit(4, "cm"),
                    #title_position = "lefttop-rot"
                  ),
                  #row_title_gp = gpar(col = c("red", "blue"), font = 1:2)
                  #left_annotation = left_anno,
                  column_title_gp = gpar(fontsize = 15), ###列标题字体修改
                  column_title = paste0("")##列标题
)

hlist2 = oncoPrint(nonresponse_top30,
                   alter_fun = alter_fun,
                   top_annotation = top.anno2,
                   #heatmap_height = unit(15, "cm"),
                   bottom_annotation = bottom.anno2,
                   #heatmap_width = unit(1000, "cm"),
                   right_annotation = NULL,
                   left_annotation = NULL,
                   col = col,
                   #gap = unit(0.1, "cm"),
                   get_type = function(x) strsplit(x, ";")[[1]],###按照分号区分重叠的类型
                   remove_empty_columns = F,
                   #column_order = col.order,
                   row_order = rownames(nonresponse_top30), ###行的顺序
                   remove_empty_rows = F, ###是否删除空的行
                   row_names_side = "left", #基因在左
                   pct_side = "right",###比例在右
                   show_pct = T, ##是否显示比例
                   pct_digits = 0, ##保留几位小数点
                   #row_labels = paste0(rownames(dat),"(10%)"),
                   #column_split = as.factor(clinical.merge$Response), ###列的分开
                   column_gap = unit(1, "mm"), ###分开的间隙
                   column_order = column_order2,
                   border=T,
                   show_column_names = F,
                   row_gap = unit(0, "mm"),
                   #row_title = NULL, ###行标签,设置为NULL即为空白
                   pct_gp = gpar(fontsize = 10,fontface="bold"), ##百分比字体
                   row_names_gp = gpar(fontsize = 10,fontface="bold"), ####基因名称字体
                   heatmap_legend_param = list(
                     #at = c(-2, 0, 2),
                     #labels = c("low", "zero", "high"),
                     title = "Mutations"####,图注标题
                     #legend_height = unit(4, "cm"),
                     #title_position = "lefttop-rot"
                   ),
                   #row_title_gp = gpar(col = c("red", "blue"), font = 1:2)
                   #left_annotation = left_anno,
                   column_title_gp = gpar(fontsize = 15), ###列标题字体修改
                   column_title = paste0("")##列标题
)

hlist_merge <- hlist1+hlist2
pdf("all_PL0_oncoplot.pdf",width = 12,height = 8)
draw(hlist_merge,
     heatmap_legend_list = Legend(labels = colnames(response_sigs), 
                                  title = "Signatures", 
                                  legend_gp = gpar(fill = c("#a27e7e","#1E4C9C","#345D82",
                                                                   "#3371B3","#5795C7","#81B5D5",
                                                                   "#AED4E5"),alpha = 0.7),
                                  title_gp = gpar(fontsize = 10),
                                  labels_gp = gpar(fontsize = 10),
                                  title_position = "topleft",ncol = 1)
)
dev.off()



