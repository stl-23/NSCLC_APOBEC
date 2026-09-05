library(openxlsx)
library(ComplexHeatmap)
library(tibble)
library(dplyr)
library(tidyverse)
source("createOncoMatrix1.R")
snv_filter <- read.xlsx("./dat.merge.filter.xlsx")

clinical <- read.xlsx("./lung_PL0_and_RNA(3).xlsx",sheet = 2)
clinical$DFS <- as.numeric(clinical$`DFS/月`)
clinical$OS <- as.numeric(clinical$`OS/月`)

clinical$DFS_status <- clinical$DFS状态
clinical$OS_status <- clinical$OS状态

clinical_dna <- clinical
clinical_dna$response_new <- ifelse(clinical_dna$response2 == "MPR","MPR","Non.MPR")
clinical_dna$type <- ifelse(clinical_dna$type == "SQC","Squamous","Non.Squamous")
clinical_dna$stage <- clinical_dna$治疗前临床分期
clinical_dna$PDL1 <- ifelse(clinical_dna$`PDL1(1=0,2=1-49%,3=>50%)` == 0 |clinical_dna$`PDL1(1=0,2=1-49%,3=>50%)` == 1 ,"PD-L1-","PD-L1+")
clinical_dna$PDL1 <- factor(clinical_dna$PDL1,levels = c("PD-L1-","PD-L1+"))

library(stringr)
sig = read.xlsx("./05.sigma/lite.signature.xlsx")
head(sig)
#sig$tumor = str_extract(sig$tumor, "Q.vcf")
sig$tumor <- gsub(".vcf","",sapply(strsplit(sig$tumor,"/"),tail,1))
#View(sig)
sigma.dat = sig[,c("tumor","Signature_clock_ml", "Signature_APOBEC_ml",
                   "Signature_28_ml", "Signature_18_ml","Signature_4_ml",
                   "Signature_3_ml", "Signature_16_ml")]
clinical_dna <- clinical_dna %>% dplyr::left_join(sigma.dat,
                                                  by=c("sample_id"="tumor"))
clinical_dna <- clinical_dna[!is.na(clinical_dna$Signature_APOBEC_ml),]


## cutoff
roc2 <- roc(response=clinical_dna$response_new,predictor=clinical_dna$Signature_APOBEC_ml,
            #smoothed = TRUE,
            #smooth = TRUE,
            levels=c("Non.MPR","MPR"),
            direction = "<") ## 手动设定 mplc < im
roc2$auc

pdf("Tissue_APOBEC_ROC.pdf",width = 4,height = 4,family = "ArialMT")
pROC::plot.roc(roc2,
               legacy.axes=TRUE,  ## 1-sepecificity if TRUE
               print.auc=TRUE,   ## add AUC
               print.thres.col="#002962", col="#002962",
               print.thres="best", print.thres.best.method="youden")  ## add best cutoff
dev.off()

bestcutoff <- coords(roc2, "best")  ## 默认为youden
bestcutoff <- bestcutoff$threshold
# 0.02603362

clinical_dna$APOBEC_group.bestcut <- ifelse(clinical_dna$Signature_APOBEC_ml > bestcutoff,"pos","neg")

## barplot
pvalue0 <- round(fisher.test(table(clinical_dna$APOBEC_group.bestcut,clinical_dna$response_new))$p.value,3)
response_table <- as.data.frame(table(clinical_dna$APOBEC_group.bestcut,clinical_dna$response_new))
names(response_table)[1:2] = c("Group","Response")
response_table$Response <- gsub("non.response","NonMPR",response_table$Response)
response_table$Response <- gsub("response","MPR",response_table$Response)
response_table$count = ifelse(response_table$Freq == 0, "", paste0("n=", response_table$Freq))
response_table$Response <- factor(response_table$Response,levels = c("NonMPR","MPR"))
response_table$Group <- factor(response_table$Group,levels = c("pos","neg"))

p0 = ggplot(response_table, aes(x = Group,y = Freq, fill = Response)) + 
  geom_bar(stat = "identity", position = "fill",width = 0.6,color = "black") +
  scale_y_continuous(labels = scales::percent,breaks = seq(0,1, 0.2)) +
  coord_flip()+
  geom_text(aes(label = count), 
            colour = "black",
            position = position_fill(vjust = 0.5),size = 6) +
  #scale_y_continuous(limits = c(0,1)) + 
  #scale_fill_manual(values = c( PD = "#bc3c28",PR = "#0072b5",SD= "#a27e7e")) +
  #scale_fill_manual(values = c( NonMPR = "#bc3c28",MPR = "#0072b5")) +
  #scale_fill_manual(values = c( NonMPR = "#B696B6",MPR = "#F6C8A8")) + 
  scale_fill_manual(values = c( NonMPR = "#2485C8",MPR = "#A27E7E")) +
  #annotate("text", x=1.5, y=1.1, parse = TRUE,label= paste0("italic(pval)==", pvalue0)) +
  annotate("richtext", x=2.5, y=0.5,label= paste0("<i>pval</i> = ", pvalue0),angle=0,
           label.color=NA,fill=NA) +
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
ggsave(p0, file = paste0("tissue_APOBEC_bar", ".pdf"), width = 5, height = 3)



## 
snv_filter <- snv_filter[snv_filter$Sample %in% clinical_dna$sample_id,]

snv_filter <- snv_filter %>% dplyr::mutate(Type=case_when(Type=="FrameShift_Deletion" ~ "Frame_Shift_Del",
                                            Type=="FrameShift_Duplication" ~ "Frame_Shift_Ins",
                                            Type=="FrameShift_Insertion" ~ "Frame_Shift_Ins",
                                            Type=="FrameShift_Substitution" ~ "Frameshift_INDEL",
                                            Type=="nonFrameShift_Deletion" ~ "In_Frame_Del",
                                            Type=="nonFrameShift_Insertion" ~ "In_Frame_Ins",
                                            Type=="nonFrameShift_Duplication" ~ "In_Frame_Ins",
                                            Type=="nonFrameShift_Substitution" ~ "In_Frame_INDEL",
                                            Type=="Nonsense_Mutation" ~ "Nonsense_Mutation",
                                            Type=="nonSynonymous_Substitution" ~ "Missense_Mutation",
                                            Type=="Splicing" ~ "Splice_Site",
                                            TRUE ~ NA_character_
))


snv_filter <- data.table::as.data.table(snv_filter)

maf.dat.mat <- createOncoMatrix(snv_filter,g=unique(snv_filter$Gene))
diffsample <- setdiff(unique(clinical_dna$sample_id),unique(snv_filter$Sample))
add_blank <- data.frame(matrix(nrow = length(diffsample),ncol = length(colnames(snv_filter))))
add_blank$X1 <- diffsample
colnames(add_blank) <- colnames(snv_filter)
snv_filter <- rbind(snv_filter,add_blank)

tmp1 <- data.frame(matrix(nrow = nrow(maf.dat.mat$oncoMatrix),ncol = length(diffsample)))
rownames(tmp1) <- rownames(maf.dat.mat$oncoMatrix)
colnames(tmp1) <- diffsample
tmp2 <- tmp1
tmp1[is.na(tmp1)] <- ""
tmp2[is.na(tmp2)] <- 0
maf.dat.mat$oncoMatrix <- cbind(maf.dat.mat$oncoMatrix,tmp1)
maf.dat.mat$numericMatrix <- cbind(maf.dat.mat$numericMatrix,tmp2)


# TMB
tmb <- as.data.frame(table(snv_filter$Sample)/1.8)
colnames(tmb) <- c("sample","TMB")

clinical_dna <-  clinical_dna %>% dplyr::left_join(tmb,by=c("sample_id"="sample"))

dat.list = maf.dat.mat$numericMatrix
dat.list = ifelse(dat.list == 0, 0,1)

## CNV
#cnv <- read.xlsx(".dat.cnv.xlsx",sheet = 1)
#cnv_filter <- cnv[cnv$Tag == "Pass",]

load("./cnv.filter.Rdata")
load("./cnv.filter.num.Rdata")

dat.list2 = cnv.dat.mat
dat.list2 = ifelse(dat.list2 == "", 0,1)
genelist2 = rowSums(dat.list2)
genelist2 = genelist2[order(genelist2, decreasing = T)]
genelist.all2 = genelist2
prop.genelist.all2 = genelist.all2/ncol(dat.list2)
row.order2 = names(genelist2)

cnv.dat.mat <-cnv.dat.mat[row.order2,]


#SNV
dat.list <- maf.dat.mat$numericMatrix
dat.list = ifelse(dat.list == 0, 0,1)
dat.list_cal <- as.data.frame(t(dat.list))
dat.list_cal <- dat.list_cal[clinical_dna$sample_id,]
dat.list_cal <- dat.list_cal[rownames(dat.list_cal) != "NA",]

## True
dat.list_cal <- as.data.frame(apply(dat.list_cal,2,function(x){return(as.numeric(x))}))
rownames(dat.list_cal) <- clinical_dna$sample_id
dat.list_cal$group <- clinical_dna$response_new

snv_pvalue <- data.frame(matrix(nrow = ncol(dat.list_cal),ncol = 2))
snv_pvalue$X1 <- colnames(dat.list_cal)
colnames(snv_pvalue) <- c("Gene","pvalue")
rownames(snv_pvalue) <- snv_pvalue$Gene
for(i in colnames(dat.list_cal)){
  my_table <- as.data.frame(table(dat.list_cal[,i],dat.list_cal$group))
  if(any(my_table$Freq) <5 | sum(my_table$Freq) < 40){
    pvalue <- fisher.test(table(dat.list_cal[,i],dat.list_cal$group))$p.value
  }else{
    pvalue <- chisq.test(table(dat.list_cal[,i],dat.list_cal$group),correct = F)$p.value
  }
  snv_pvalue[i,"pvalue"] <- pvalue
}
snv_pvalue <- arrange(snv_pvalue,pvalue)
snv_pvalue <- snv_pvalue[snv_pvalue$Gene != "group",]
snv_pvalue$FDR <- p.adjust(snv_pvalue$pvalue,method = "BH")
write.xlsx(snv_pvalue,file = "snv_pvalue_rm_NA.xlsx")

# CNV
dat.list_cal2 <- as.data.frame(t(dat.list2))
dat.list_cal2 <- dat.list_cal2[rownames(dat.list_cal2) %in% clinical_dna$sample_id,]
#group_cnv <- group[group$DNA %in% inter_sample,]
lacked2 <- setdiff(clinical_dna$sample_id,rownames(dat.list_cal2))

for(i in lacked2){
  rows <- rep(0,ncol(dat.list_cal2))
  names(rows) <- colnames(dat.list_cal2)
  dat.list_cal2 <- rbind(dat.list_cal2,rows)
}
rownames(dat.list_cal2)[c(21:33)] <- lacked2
dat.list_cal2 <- dat.list_cal2[clinical_dna$sample_id,]

identical(rownames(dat.list_cal2),clinical_dna$sample_id)
## True
dat.list_cal2 <- as.data.frame(apply(dat.list_cal2,2,function(x){return(as.numeric(x))}))
rownames(dat.list_cal2) <- clinical_dna$sample_id
dat.list_cal2$group <- clinical_dna$response_new

cnv_pvalue <- data.frame(matrix(nrow = ncol(dat.list_cal2),ncol = 2))
cnv_pvalue$X1 <- colnames(dat.list_cal2)
colnames(cnv_pvalue) <- c("Gene","pvalue")
rownames(cnv_pvalue) <- cnv_pvalue$Gene
for(i in colnames(dat.list_cal2)[-ncol(dat.list_cal2)]){
  my_table <- as.data.frame(table(dat.list_cal2[,i],dat.list_cal2$group))
  if(sum(dat.list_cal2[,i]) >0){
    if(any(my_table$Freq) <5 | sum(my_table$Freq) < 40){
      pvalue <- fisher.test(table(dat.list_cal2[,i],dat.list_cal2$group))$p.value
    }else{
      pvalue <- chisq.test(table(dat.list_cal2[,i],dat.list_cal2$group),correct = F)$p.value
    }
    cnv_pvalue[i,"pvalue"] <- pvalue
  }}
cnv_pvalue <- arrange(cnv_pvalue,pvalue)
cnv_pvalue <- cnv_pvalue[cnv_pvalue$Gene != "group",]
cnv_pvalue$FDR <- p.adjust(cnv_pvalue$pvalue,method = "BH")
write.xlsx(cnv_pvalue,file = "cnv_pvalue_rmNA.xlsx")
#write.table(cnv_pvalue,file = "cnv_pvalue.xls",quote = F,row.names = F,sep = "\t")


##### oncoplot
snv_matrix <- maf.dat.mat$oncoMatrix[,colnames(maf.dat.mat$oncoMatrix) %in% clinical_dna$sample_id]
snv_matrix <- as.data.frame(snv_matrix)


#snv_top20 <- head(snv_matrix,n=20)

cnv_matrix <- cnv.dat.mat[,colnames(cnv.dat.mat) %in% clinical_dna$sample_id]

for(i in lacked2){
  cnv_matrix[,i] <- ""
}

# group
response_group <- clinical_dna[clinical_dna$response_new == "MPR",]
nonresponse_group <- clinical_dna[clinical_dna$response_new == "Non.MPR",]

snv_matrix_response <- snv_matrix[,colnames(snv_matrix) %in% response_group$sample_id]
snv_matrix_nonresponse <- snv_matrix[,colnames(snv_matrix) %in% nonresponse_group$sample_id]

response_group <- response_group[match(colnames(snv_matrix_response),response_group$sample_id),]
nonresponse_group <- nonresponse_group[match(colnames(snv_matrix_nonresponse),nonresponse_group$sample_id),]
clinical_dna <- rbind(response_group,nonresponse_group)

cnv_matrx_response <- cnv_matrix[,colnames(snv_matrix_response)]
cnv_matrx_nonresponse <- cnv_matrix[,colnames(snv_matrix_nonresponse)]

snv_matrix_response_top20 <- head(snv_matrix_response,n=20)
cnv_matrx_response_top15 <- head(cnv_matrx_response,n=15)
response_top20 <- rbind(snv_matrix_response_top20,cnv_matrx_response_top15)

snv_matrix_nonresponse_top20 <- head(snv_matrix_nonresponse,n=20)
cnv_matrx_nonresponse_top15 <- head(cnv_matrx_nonresponse,n=15)
nonresponse_top20 <- rbind(snv_matrix_nonresponse_top20,cnv_matrx_nonresponse_top15 )



response_group <- clinical_dna[clinical_dna$response_new == "MPR",]
nonresponse_group <- clinical_dna[clinical_dna$response_new == "Non.MPR",]

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
  'Splice_Site',
  "amplification"
)
col[9] <- "black"  ## Multi_Hit 
col[11] <- "#976666" ## amplication

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
  },
  amplification  = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.8, "mm"), h-unit(0.5, "mm"),
              gp = gpar(fill = col["amplification"], col = NA))
  }
  
)

cut_group <- c(rep("SNV",20),rep("CNV",15))
cut_group <- factor(cut_group,levels = c("SNV","CNV"))

response_tmb <- tmb[tmb$sample %in% response_group$sample_id,]
response_tmb <- response_tmb[order(match(response_tmb$sample,colnames(response_top20))),]
response_tmb$sample <- as.character(response_tmb$sample)

nonresponse_tmb <- tmb[tmb$sample %in% nonresponse_group$sample_id,]
nonresponse_tmb <- nonresponse_tmb[order(match(nonresponse_tmb$sample,colnames(nonresponse_top20))),]
nonresponse_tmb$sample <- as.character(nonresponse_tmb$sample)


library(randomcoloR)

stage_color <- distinctColorPalette(3)
#treatment_color <- c("#AEB38C","#E7E359","#E35786")
#response_color <- c("#5FB155","#7BBCE1","#DD9489")
response_color <- c("#5FB155","#DD9489")
#type_color <- distinctColorPalette(20)[c(8,10,12)]
#type_color <- c("#E1C0A6","#863CE2","#D185DB")
type_color <- c("#E1C0A6","#863CE2")
names(stage_color) = names(table(clinical_dna$stage))
#names(treatment_color) = names(table(clinical.merge$Treatment))
names(response_color) = c("MPR","Non.MPR")
names(type_color) = names(table(clinical_dna$type))

pdl1_color <- c("lightblue","#E35786")
names(pdl1_color) = names(table(clinical_dna$PDL1))

top.anno1 <- columnAnnotation(
  TMB = anno_barplot(as.numeric(response_tmb$TMB),ylim = c(0,20),  axis_param = list(at = seq(0,20, by = 10)),
                     border = F),
  Group = as.character(response_group$response_new),
  Subtype = as.character(response_group$type),
  Stage = as.character(response_group$stage),
  PDL1 = as.character(response_group$PDL1),
  #institution = factor(clinical.merge$institution,levels = c("ZH","AD")),
  foo = anno_empty(border = F,height = unit(c(0.03), "cm")), ###加入空注释，以调整注释与热图之间的间距
  simple_anno_size = unit(c(0.3), "cm"), ###每个注释的宽度
  col = list(Group = response_color,
             Subtype = type_color,
             Stage = stage_color,
             PDL1 = pdl1_color
  ),
  #axis_param = list(at = c(0,5,10,15,20)),
  show_annotation_name = T,
  annotation_name_side = "left", ###注释的标签在左还是右
  annotation_name_gp = gpar(fontsize = 10) ####标签的参数修改
  #gap = unit(c(0.5,0.5), "mm") ###各个注释之间的间距
)

top.anno2 <- columnAnnotation(
  TMB = anno_barplot(as.numeric(nonresponse_tmb$TMB),ylim = c(0,20),  axis_param = list(at = seq(0,20, by = 10)),
                     border = F),
  Group = as.character(nonresponse_group$response_new),
  Subtype = as.character(nonresponse_group$type),
  Stage = as.character(nonresponse_group$stage),
  PDL1 = as.character(nonresponse_group$PDL1),
  #institution = factor(clinical.merge$institution,levels = c("ZH","AD")),
  foo = anno_empty(border = F,height = unit(c(0.03), "cm")), ###加入空注释，以调整注释与热图之间的间距
  simple_anno_size = unit(c(0.3), "cm"), ###每个注释的宽度
  col = list(Group = response_color,
             Subtype = type_color,
             Stage = stage_color,
             PDL1 = pdl1_color
  ),
  #axis_param = list(at = c(0,5,10,15,20)),
  show_annotation_name = T,
  annotation_name_side = "left", ###注释的标签在左还是右
  annotation_name_gp = gpar(fontsize = 10) ####标签的参数修改
  #gap = unit(c(0.5,0.5), "mm") ###各个注释之间的间距
)



left_anno = rowAnnotation(foo = anno_block(gp = gpar(fill = c(4,2)),
                                           #  labels = c("Mutations"), 
                                           labels_gp = gpar(col = "white", fontsize = 15)))


response_sigs <- response_group[,c(2,44:50)]
rownames(response_sigs) <- response_sigs$sample_id
response_sigs <- response_sigs[,-1]
colnames(response_sigs) <- gsub("_ml","",colnames(response_sigs))
response_sigs <- apply(response_sigs,2,function(x){x*100})
response_sigs <- response_sigs[,c(2,1,3:7)]

nonresponse_sigs <- nonresponse_group[,c(2,44:50)]
rownames(nonresponse_sigs) <- nonresponse_sigs$sample_id
nonresponse_sigs <- nonresponse_sigs[,-1]
colnames(nonresponse_sigs) <- gsub("_ml","",colnames(nonresponse_sigs))
nonresponse_sigs <- apply(nonresponse_sigs,2,function(x){x*100})
nonresponse_sigs[is.na(nonresponse_sigs)] <- 0
nonresponse_sigs <- nonresponse_sigs[,c(2,1,3:7)]

library(viridis)
# https://waldyrious.net/viridis-palette-generator/
bottom.anno1 <- HeatmapAnnotation(
  empty = anno_empty(border = FALSE,height = unit(0, "cm")),
  Signature = anno_barplot(response_sigs,
                           #gp = gpar(fill = c("#fde725","#90d743",
                          #                     "#35b779","#21918c", 
                           #                   "#31688e","#443983","#440154")
                          #                    ),
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
                           #gp = gpar(fill = c("#fde725","#90d743",
                          #                    "#35b779","#21918c", 
                           #                   "#31688e","#443983","#440154")
                          # ),
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

hlist1 = oncoPrint(response_top20,
                   alter_fun = alter_fun,
                   top_annotation = top.anno1,
                   #heatmap_height = unit(15, "cm"),
                   #bottom_annotation = bottom.anno1,
                   #heatmap_width = unit(1000, "cm"),
                   left_annotation = left_anno,
                   right_annotation = NULL,
                   #left_annotation = left_anno,
                   bottom_annotation = bottom.anno1,
                   col = col,
                   #gap = unit(0.1, "cm"),
                   get_type = function(x) strsplit(x, ";")[[1]],###按照分号区分重叠的类型
                   remove_empty_columns = F,
                   #column_order = rownames(mut_top20),
                   #column_split = col_group,
                   row_order =  rownames(response_top20), ###行的顺序
                   row_split = cut_group,
                   remove_empty_rows = F, ###是否删除空的行
                   row_names_side = "left", #基因在左
                   pct_side = "right",###比例在右
                   show_pct = T, ##是否显示比例
                   pct_digits = 0, ##保留几位小数点
                   #row_labels = paste0(rownames(dat),"(10%)"),
                   #column_split = as.factor(clinical.merge$group), ###列的分开
                   #column_gap = unit(2, "mm"), ###分开的间隙
                   #column_order = colnames(dat2_top10),
                   border=T,
                   show_column_names = F,
                   row_gap = unit(2, "mm"),
                   #row_title = NULL, ###行标签,设置为NULL即为空白
                   pct_gp = gpar(fontsize = 10,fontface="bold"), ##百分比字体
                   row_names_gp = gpar(fontsize = 10,fontface="bold"), ####基因名称字体
                   #column_names_gp = gpar(fontsize=8,fontface="bold"),
                   heatmap_legend_param = list(
                     #at = c(-2, 0, 2),
                     #labels = c("low", "zero", "high"),
                     title = "Alterations"####,图注标题
                     #legend_height = unit(4, "cm"),
                     #title_position = "lefttop-rot"
                   ),
                   #row_title_gp = gpar(col = c("red", "blue"), font = 1:2)
                   #left_annotation = left_anno,
                   column_title_gp = gpar(fontsize = 15), ###列标题字体修改
                   column_title = paste0("")##列标题
)



hlist2 = oncoPrint(nonresponse_top20,
                   alter_fun = alter_fun,
                   top_annotation = top.anno2,
                   #heatmap_height = unit(15, "cm"),
                   #bottom_annotation = bottom.anno,
                   #heatmap_width = unit(1000, "cm"),
                   left_annotation = NULL,
                   right_annotation = NULL,
                   #left_annotation = left_anno,
                   bottom_annotation = bottom.anno2,
                   col = col,
                   #gap = unit(0.1, "cm"),
                   get_type = function(x) strsplit(x, ";")[[1]],###按照分号区分重叠的类型
                   remove_empty_columns = F,
                   #column_order = rownames(mut_top20),
                   #column_split = col_group,
                   row_order =  rownames(response_top20), ###行的顺序
                   row_split = cut_group,
                   remove_empty_rows = F, ###是否删除空的行
                   row_names_side = "left", #基因在左
                   pct_side = "right",###比例在右
                   show_pct = T, ##是否显示比例
                   pct_digits = 0, ##保留几位小数点
                   #row_labels = paste0(rownames(dat),"(10%)"),
                   #column_split = as.factor(clinical.merge$group), ###列的分开
                   #column_gap = unit(2, "mm"), ###分开的间隙
                   #column_order = colnames(dat2_top10),
                   border=T,
                   show_column_names = F,
                   row_gap = unit(2, "mm"),
                   #row_title = NULL, ###行标签,设置为NULL即为空白
                   pct_gp = gpar(fontsize = 10,fontface="bold"), ##百分比字体
                   row_names_gp = gpar(fontsize = 10,fontface="bold"), ####基因名称字体
                   #column_names_gp = gpar(fontsize=8,fontface="bold"),
                   heatmap_legend_param = list(
                     #at = c(-2, 0, 2),
                     #labels = c("low", "zero", "high"),
                     title = "Alterations"####,图注标题
                     #legend_height = unit(4, "cm"),
                     #title_position = "lefttop-rot"
                   ),
                   #row_title_gp = gpar(col = c("red", "blue"), font = 1:2)
                   #left_annotation = left_anno,
                   column_title_gp = gpar(fontsize = 15), ###列标题字体修改
                   column_title = paste0("")##列标题
)

hlist_merge <- hlist1+hlist2
pdf("oncoplot_snv_cnv_v6.pdf",width = 12,height = 8)
draw(hlist_merge,
      heatmap_legend_list = Legend(labels = colnames(response_sigs), 
                                   title = "Signatures", 
                                   legend_gp = gpar(fill = c("#a27e7e","#1E4C9C","#345D82",
                                                             "#3371B3","#5795C7","#81B5D5",
                                                             "#AED4E5"),alpha=0.7
                                                      ),
                                   title_gp = gpar(fontsize = 10),
                                   labels_gp = gpar(fontsize = 10),
                                   title_position = "topleft",ncol = 1)
)
dev.off()




