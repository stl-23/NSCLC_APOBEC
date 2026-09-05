library(openxlsx)
library(tibble)
library(dplyr)
library(tidyverse)
library(finalfit)

load("../baseline_DNA_mut.Rdata")
#clinical <- read.xlsx("/Users/stl/Documents/AmoyDX/6.work_dir/10.lung_immune_adu/00.data/lung_38_PL0_and_34_RNA_v2.xlsx")
clinical <- read.xlsx("./lung_PL0_and_RNA(3).xlsx",sheet = 1)
clinical_immuno <- clinical  ## clincial factors use all samples


clinical_immuno$response_new <- ifelse(clinical_immuno$response == "Non-MPR","non.response","response")
clinical_immuno$response_new <- factor(clinical_immuno$response_new,levels = c("response","non.response"))
clinical_immuno$stage <- clinical_immuno$治疗前临床分期

tmb <- read.xlsx("./TMB.xlsx")
clinical_immuno <-  clinical_immuno %>% dplyr::left_join(tmb,by=c("sample_id"="sample"))

library(stringr)
sig = read.xlsx("./05.sigma/lite.signature.xlsx")
head(sig)
#sig$tumor = str_extract(sig$tumor, "Q.vcf")
sig$tumor <- gsub(".vcf","",sapply(strsplit(sig$tumor,"/"),tail,1))
#View(sig)
sigma.dat = sig[,c("tumor","Signature_clock_ml", "Signature_APOBEC_ml",
                   "Signature_28_ml", "Signature_18_ml","Signature_4_ml",
                   "Signature_3_ml", "Signature_16_ml")]
#colnames(sigma.dat) <- gsub("_ml","",colnames(sigma.dat))
clinical_immuno <- clinical_immuno %>% dplyr::left_join(sigma.dat,
                                                        by=c("sample_id"="tumor"))

clinical_immuno$Response_new <- ifelse(clinical_immuno$response_new == "non.response",1,0)
clinical_immuno$性别.1男 <- trimws(clinical_immuno$性别.1男)
clinical_immuno$gender <- ifelse(clinical_immuno$性别 == "1",1,0)
clinical_immuno$age <- as.numeric(clinical_immuno$年龄)
clinical_immuno$smoking <- as.numeric(clinical_immuno$`吸烟史（0：无；1：有）`)
clinical_immuno$Cancer_types <- clinical_immuno$type

clinical_immuno$smoking <- factor(clinical_immuno$smoking,levels = c(0,1))

clinical_immuno$smoking.factor <- ifelse(clinical_immuno$smoking == 1,"Yes","No")
clinical_immuno$smoking.factor <- factor(clinical_immuno$smoking.factor,levels = c("No","Yes"))

clinical_immuno$clinical.N.new.factor <- clinical_immuno$clinical.N.new
clinical_immuno$clinical.N.new.factor <- ifelse(clinical_immuno$clinical.N.new.factor ==0 | clinical_immuno$clinical.N.new.factor ==1 ,"N0/1",
                                                ifelse(clinical_immuno$clinical.N.new.factor ==2,"N2","N3"))
clinical_immuno$clinical.N.new.factor <- factor(clinical_immuno$clinical.N.new.factor,levels = c("N0/1","N2","N3"))


clinical_immuno$clinical.T.new.factor <- clinical_immuno$clinical.T.new
clinical_immuno$clinical.T.new.factor <- ifelse(clinical_immuno$clinical.T.new.factor ==1,"T1",
                                                ifelse(clinical_immuno$clinical.T.new.factor ==2,"T2",
                                                       ifelse(clinical_immuno$clinical.T.new.factor ==3,"T3","T4")))



clinical_immuno$clinical.T.new.factor <- factor(clinical_immuno$clinical.T.new.factor,levels = c("T1","T2","T3","T4"))

clinical_immuno$Cancer_types.factor <- ifelse(clinical_immuno$Cancer_types == "SQC","Squamous","Non.squamous")
clinical_immuno$Cancer_types.factor <- factor(clinical_immuno$Cancer_types.factor,levels = c("Non.squamous","Squamous"))

clinical_immuno$PDL1.factor <- ifelse(clinical_immuno$`PDL1(1=0,2=1-49%,3=>50%)` == 0 |clinical_immuno$`PDL1(1=0,2=1-49%,3=>50%)` == 1 ,"-","+")
clinical_immuno$PDL1.factor <- factor(clinical_immuno$PDL1.factor,levels = c("-","+"))

clinical_immuno$PDL1.factor2 <- ifelse(clinical_immuno$`PDL1(1=0,2=1-49%,3=>50%)` == 0 |clinical_immuno$`PDL1(1=0,2=1-49%,3=>50%)` == 1 ,"1",
                                       ifelse(clinical_immuno$`PDL1(1=0,2=1-49%,3=>50%)` == 2,"2","3"))
clinical_immuno$PDL1.factor2 <- factor(clinical_immuno$PDL1.factor2,levels = c("1","2","3"))

clinical_immuno$TMB_group <- ifelse(clinical_immuno$TMB >= 10,"High(>=10)","Low(<10)")
clinical_immuno$TMB_group <- factor(clinical_immuno$TMB_group,levels = c("Low(<10)","High(>=10)"))

save(clinical_immuno,file = "clinical_immuno.rdata")

## 
explanatory1 <- colnames(clinical_immuno)[c(59,49,52,54,53,55,56,58)]
dependent1 <- "response_new"
Uni_glm1 <- clinical_immuno %>% 
  finalfit(dependent1,explanatory1) 


explanatory2 <- colnames(clinical_dna)[c(59,49,52,54,53,55,56,58)] ## APOBEC and clinical
Uni_glm2 <- clinical_dna %>% 
  finalfit(dependent2,explanatory1)


## single and multiple together

rs_forest1 <- Uni_glm1
library(stringr)
rs_forest1$OR = str_extract(rs_forest1$`OR (univariable)`, "^.+?(?=\\()") %>% as.numeric()
rs_forest1$OR.95L = str_extract(rs_forest1$`OR (univariable)`, "(?<=\\().+?(?=-)")%>% as.numeric()
rs_forest1$OR.95H = str_extract(rs_forest1$`OR (univariable)`, "(?<=-).+?(?=,)")%>% as.numeric()
rs_forest1$P = str_extract(rs_forest1$`OR (univariable)`, "(?<=\\, ).+?(?=\\))")

rs_forest1 <- rs_forest1[-c(17:18,21:22),]
rs_forest1 <- rs_forest1[c(17:18,1:16),]

rs_forest2 <- Uni_glm2
library(stringr)
rs_forest2$OR = str_extract(rs_forest2$`OR (multivariable)`, "^.+?(?=\\()") %>% as.numeric()
rs_forest2$OR.95L = str_extract(rs_forest2$`OR (multivariable)`, "(?<=\\().+?(?=-)")%>% as.numeric()
rs_forest2$OR.95H = str_extract(rs_forest2$`OR (multivariable)`, "(?<=-).+?(?=,)")%>% as.numeric()
rs_forest2$P = str_extract(rs_forest2$`OR (multivariable)`, "(?<=\\, ).+?(?=\\))")

identical(rs_forest1$`Dependent: response_new`,rs_forest2$`Dependent: response_new`)

rs_forest <- cbind(rs_forest1[,c(1,2,5)],rs_forest2[,c(6),drop=FALSE])


library(stringr)
rs_forest$OR.uni = str_extract(rs_forest$`OR (univariable)`, "^.+?(?=\\()") %>% as.numeric()
rs_forest$OR.95L.uni = str_extract(rs_forest$`OR (univariable)`, "(?<=\\().+?(?=-)")%>% as.numeric()
rs_forest$OR.95H.uni = str_extract(rs_forest$`OR (univariable)`, "(?<=-).+?(?=,)")%>% as.numeric()
rs_forest$P.uni = str_extract(rs_forest$`OR (univariable)`, "(?<=\\, ).+?(?=\\))")

rs_forest$OR.mul = str_extract(rs_forest$`OR (multivariable)`, "^.+?(?=\\()") %>% as.numeric()
rs_forest$OR.95L.mul = str_extract(rs_forest$`OR (multivariable)`, "(?<=\\().+?(?=-)")%>% as.numeric()
rs_forest$OR.95H.mul = str_extract(rs_forest$`OR (multivariable)`, "(?<=-).+?(?=,)")%>% as.numeric()
rs_forest$P.mul = str_extract(rs_forest$`OR (multivariable)`, "(?<=\\, ).+?(?=\\))")

dataset <- data.frame(
  Characteristic = rs_forest$`Dependent: response_new`,
  Subgroup = rs_forest[,2],
  OR.uni = rs_forest$OR.uni,
  Lower.uni= rs_forest$OR.95L.uni,
  Upper.uni=rs_forest$OR.95H.uni,
  p.Value.uni = rs_forest$P.uni,
  estimate_lab1 = paste0(rs_forest$OR.uni, "(", rs_forest$OR.95L.uni, 
                         "-", rs_forest$OR.95H.uni,")"),
  OR.mul = rs_forest$OR.mul,
  Lower.mul= rs_forest$OR.95L.mul,
  Upper.mul=rs_forest$OR.95H.mul,
  p.Value.mul = rs_forest$P.mul,
  estimate_lab2 = paste0(rs_forest$OR.mul, "(", rs_forest$OR.95L.mul, 
                         "-", rs_forest$OR.95H.mul,")")
)

#dataset$Characteristic[c(1,3,5,7,9,11,13)] <- c("Risk","Age","Stage","Type","PARPi",
#                                                "BRCA.status","ORD.status")

dataset[,13] <- paste(rep(" ", 25), collapse = " ")
colnames(dataset)[c(7,12,13)] <- c("OR(95%CI).uni","OR(95%CI).mul","") 

dataset$`OR(95%CI).uni`[dataset$`OR(95%CI).uni` == "NA(NA-NA)"] = "Ref"
dataset$p.Value.uni[is.na(dataset$p.Value.uni)] <- "-"
dataset$`OR(95%CI).mul`[dataset$`OR(95%CI).mul` == "NA(NA-NA)"] = "Ref"
dataset$p.Value.mul[is.na(dataset$p.Value.mul)] <- "-"
dataset[c(2,11,15),9] <- c(0.001,0.001,0.001)
dataset <- dataset[-c(10:12),]

dataset$Characteristic[c(1,3,4,6,10,13,15)] <- c("NSCLC APOBEC","Age","Smoking","Clinical T",
                                              "Cancer type","PD-L1","TMB")

library(forestploter)
## https://blog.csdn.net/dege857/article/details/127859291
f_theme <- forest_theme(
  core = list(bg_params=list(fill = c("white"))), ## backgroud color
  base_size = 10,  #文本的大小
  # Confidence interval point shape, line type/color/width
  ci_pch = 15,   #可信区间点的形状
  ci_col = "#762a83",    #CI的颜色
  ci_fill = "#0072B5FF",     #ci颜色填充
  ci_alpha = 0.8,        #ci透明度
  ci_lty = 1,            #CI的线型
  ci_lwd = 1.5,          #CI的线宽
  ci_Theight = 0.2, # Set an T end at the end of CI  ci的高度，默认是NULL
  # Reference line width/type/color   参考线默认的参数，中间的竖的虚线
  refline_lwd = 1,       #中间的竖的虚线
  refline_lty = "dashed",
  refline_col = "grey20",
  # Vertical line width/type/color  垂直线宽/类型/颜色   可以添加一条额外的垂直线，如果没有就不显示
  vertline_lwd = 1,              #可以添加一条额外的垂直线，如果没有就不显示
  vertline_lty = "dashed",
  vertline_col = "grey20",
  # Change summary color for filling and borders   更改填充和边框的摘要颜色
  summary_fill = "yellow",       #汇总部分大菱形的颜色
  summary_col = "#4575b4",
  # Footnote font size/face/color  脚注字体大小/字体/颜色
  footnote_cex = 0.6,
  footnote_fontface = "italic",
  footnote_col = "red",
  title_just = "center"  ## 标题居中
)


p3 <-forest(dataset[,c(1:2,13,7,6,13,12,11)],
            est=list(dataset$OR.uni,
                     dataset$OR.mul),
            lower=list(as.numeric(dataset$Lower.uni),
                       as.numeric(dataset$Lower.mul)),
            upper=list(as.numeric(dataset$Upper.uni),
                       as.numeric(dataset$Upper.mul)),
            sizes = 1.5,
            ticks_at = c(0.1,0.5,1,8,64),
            ci_column = c(3,6),
            x_trans = "log2",
            #xlab = "Odds ratio (log2 scale)",
            arrow_lab = c("MPR","NonMPR"),
            title = "",
            theme = f_theme
)

pdf("clinical_apobec.Logistic_v8.pdf",width = 12,height = 5,family = "ArialMT")
print(p3)
dev.off()
