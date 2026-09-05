library(openxlsx)
library(ggsci)
library(dplyr)
library(tidyverse)


load("./01.PL0_oncoplot/1.ctDNA_oncoplot_rm_NSCLC/add_2_apobec_neg_samples/group_2025.5.19.rdata")


group$APOBEC_group.bestcut <- ifelse(group$APOBEC_group.bestcut == "APOBEC-","Low","High")

group$APOBEC_group.bestcut <- factor(group$APOBEC_group.bestcut,levels = c("Low","High"))
group$smoking.factor <- ifelse(group$smoking == 1,"Yes","No")
group$smoking.factor <- factor(group$smoking.factor,levels = c("No","Yes"))

group$clinical.N.new.factor <- group$clinical.N.new
group$clinical.N.new.factor <- factor(group$clinical.N.new.factor)

group$clinical.T.new.factor <- group$clinical.T.new

group$clinical.T.new.factor <- factor(group$clinical.T.new.factor)

group$Cancer_types.factor <- group$type
group$Cancer_types.factor <- factor(group$Cancer_types.factor,levels = c("Non.Squamous","Squamous"))

group$PDL1.factor <- ifelse(group$PDL1 == "PD-L1-","-","+")
group$PDL1.factor <- factor(group$PDL1.factor,levels = c("-","+"))
group$TMB_group <- ifelse(group$TMB >=16,"High(>=16)","Low(<16)")
group$TMB_group <- factor(group$TMB_group,levels =c("Low(<16)","High(>=16)"))



source("./logistic.R")

group_new <- group[,c(37,53,49,54,56,57,58,59)]
colnames(group_new) <- c("Response","Blood.APOBEC","Age","Smoking","Clinical.T","Cancer.type","PDL1","bTMB")

explanatory1 <- colnames(group_new)[c(2:8)]
dependent1 <- "Response"


All_logistic2(data = group_new,response_var = dependent1,predictor_vars = explanatory1,
              neg_prefer = "MPR",pos_prefer = "Non.MPR",width=6,height=5,cellwidth=2,cellheight=3,
              filename = "ctDNA_clinical_logistic.pdf")


