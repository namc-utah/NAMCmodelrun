library(randomForest)
library(labdsv)
library(vegan)


#AREMP
setwd("Z:\\buglab\\OE_Modeling\\NAMC_Supported_OEmodels\\AREMP2014\\ModelFiles")
load("AREMP_OE.rdata")
load("AREMP_MMI.rdata")


test_preds=read.csv(file="Z:\\buglab\\OE_Modeling\\NAMC_Supported_OEmodels\\AREMP2014\\ModelFiles\\InputsAndResults\\2021\\Test_preds.csv",row.names="SampleID")
bugsOTU_matrix_list=read.csv(file="Z:\\buglab\\OE_Modeling\\NAMC_Supported_OEmodels\\AREMP2014\\ModelFiles\\InputsAndResults\\2021\\bugsOTU_matrix.csv")
#Matching predictors and bugs
bugsOTU_matrix_matrix=matrify(bugsOTU_matrix_list)
bugsOTU_matrix_matrix=bugsOTU_matrix_matrix[row.names(test_preds),]
if(any((row.names(bugsOTU_matrix_matrix)==row.names(test_preds))=="FALSE")=="TRUE"){print("Warning - predictor/bug mismatch")}

#Checking distributions of new predictors against reference predictors
boxplot(refpreds$WSA_SQKM,test_preds$WSA_SQKM, names=c("Ref","New"),ylab="Watershed Area (km2)")
boxplot(refpreds$TMAX_WS,test_preds$TMAX_WS, names=c("Ref","New"),ylab="TMAX_WS (C)")
boxplot(refpreds$ELVmean_WS,test_preds$ELVmean_WS, names=c("Ref","New"),ylab="ELVmean_WS (m)")
boxplot(refpreds$ELVmax_WS,test_preds$ELVmax_WS, names=c("Ref","New"),ylab="ELVmax_WS (m)")
boxplot(refpreds$PMIN_WS,test_preds$PMIN_WS, names=c("Ref","New"),ylab="PMIN_WS (mm)")

boxplot(OE.assess.test$OE.scores$OoverE,ylab="New sample O/E scores")
plot(OE.assess.test$OE.scores$OoverE~resamp_vec,ylab="New sample O/E scores",xlab="Resample count")

setwd("Z:\\buglab\\OE_Modeling\\NAMC_Supported_OEmodels\\AREMP2014\\ModelFiles")

#Matching predictors and bugs
test_metrics_raw=test_metrics_raw[,colnames(ref_metrics_adj)]
test_preds=test_preds[row.names(test_metrics_raw),]
if(any((row.names(test_metrics_raw)==row.names(test_preds))=="FALSE")=="TRUE"){print("Warning - predictor/bug mismatch")}

#Checking distributions of new predictors against reference predictors
boxplot(ref_preds$ELVmean_WS,test_preds$ELVmean_WS, names=c("Ref","New"),ylab="ELVmean_WS (m)")
boxplot(ref_preds$ELVmin_WS,test_preds$ELVmin_WS, names=c("Ref","New"),ylab="ELVmin_WS (m)")
boxplot(ref_preds$KFCT_WS,test_preds$KFCT_WS, names=c("Ref","New"),ylab="KFCT_WS")
boxplot(ref_preds$PMIN_PT,test_preds$PMIN_PT, names=c("Ref","New"),ylab="PMIN_PT")
boxplot(ref_preds$PMIN_WS,test_preds$PMIN_WS, names=c("Ref","New"),ylab="PMIN_WS")
boxplot(ref_preds$rh_WS,test_preds$rh_WS, names=c("Ref","New"),ylab="rh_WS")
boxplot(ref_preds$TMAX_WS,test_preds$TMAX_WS, names=c("Ref","New"),ylab="TMAX_WS")
boxplot(ref_preds$TMEAN_PT,test_preds$TMEAN_PT, names=c("Ref","New"),ylab="TMEAN_PT")
boxplot(ref_preds$TMEAN_WS,test_preds$TMEAN_WS, names=c("Ref","New"),ylab="TMEAN_WS")

boxplot(ref_scores,test_MMI_scores,mdeg_cal_scores,mdeg_val_scores,ylab="MMI score",names=c("Ref","Test","Mdeg cal","Mdeg val"))
abline(h=thresholds[1],lwd=3,col="red")
abline(h=thresholds[2],lwd=3,col="blue")


#NV

load("OE_MMI_models.rdata")
#modeling, adjusting, rescaling
predictors=test_pred
bugsOTU_matrix=bugsOTU_matrix[row.names(predictors),]
namecheck=row.names(predictors)==row.names(bugsOTU_matrix)
if (any(namecheck=="FALSE")) stop("Predictors do not match bug metrics")

# OR

load("MWCF.rdata")
load('Nov05model_MWCF_16jan13.Rdata')
source("model.predict.v4.1.R")
library(plyr)

load("WCCP.RData")
source("model.predict.v4.1.R")
load('Nov05model_WCCP_16jan13.RData')

#PIBO
library(cluster); library(Hmisc);
library(randomForest);
library(plyr);
library(dplyr);
source("\\\\share1.bluezone.usu.edu\\miller\\\\buglab\\OE_Modeling\\NAMC_Supported_OEmodels\\PIBO\\PIBO_RF2020\\model.predict.RanFor.r")
load("\\\\share1.bluezone.usu.edu\\miller\\\\buglab\\OE_Modeling\\NAMC_Supported_OEmodels\\PIBO\\PIBO_RF2020\\Benkendorf.RF.Model.Version1.Rdata")

#UT

library(randomForest)
library(labdsv)
library(vegan)
load('Model files\\UTDEQ_15_OE_model.rdata')


#WY
library(gtools); library(MASS); library(cluster); library(Hmisc);
source("model.predict.v4.1.r");
load('WYRIVPACS2012.Rdata')


#WW
load(paste0("My.RF.Model.Version1.Rdata"))
source(paste0("model.predict.RanFor.4.2.r"));
