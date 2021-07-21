###Apply AREMP model to test sites:
setwd("Z:\\buglab\\OE_Modeling\\NAMC_Supported_OEmodels\\AREMP2014\\ModelFiles")
load("AREMP_OE.rdata")

library(randomForest)
library(labdsv)
library(vegan)

#Read in test data:
#test_preds=read.csv(file="test_preds.csv",row.names="SampleID")
#test_bugs_list=read.csv(file="test_bugs_list.csv")

##test_preds=read.csv(file="Z:\\buglab\\OE_Modeling\\NAMC_Supported_OEmodels\\AREMP2014\\InputsAndResults\\2015\\habitat_input.csv",row.names="SampleID")
##test_bugs_list=read.csv(file="Z:\\buglab\\OE_Modeling\\NAMC_Supported_OEmodels\\AREMP2014\\InputsAndResults\\2015\\OE_bugs_list_input.csv")

test_preds=read.csv(file="Z:\\buglab\\OE_Modeling\\NAMC_Supported_OEmodels\\AREMP2014\\ModelFiles\\InputsAndResults\\2021\\Test_preds.csv",row.names="SampleID")
test_bugs_list=read.csv(file="Z:\\buglab\\OE_Modeling\\NAMC_Supported_OEmodels\\AREMP2014\\ModelFiles\\InputsAndResults\\2021\\Test_bugs.csv")

#Matching predictors and bugs
test_bugs_matrix=matrify(test_bugs_list)
test_bugs_matrix=test_bugs_matrix[row.names(test_preds),]
if(any((row.names(test_bugs_matrix)==row.names(test_preds))=="FALSE")=="TRUE"){print("Warning - predictor/bug mismatch")}

#Checking distributions of new predictors against reference predictors
boxplot(refpreds$WSA_SQKM,test_preds$WSA_SQKM, names=c("Ref","New"),ylab="Watershed Area (km2)")
boxplot(refpreds$TMAX_WS,test_preds$TMAX_WS, names=c("Ref","New"),ylab="TMAX_WS (C)")
boxplot(refpreds$ELVmean_WS,test_preds$ELVmean_WS, names=c("Ref","New"),ylab="ELVmean_WS (m)")
boxplot(refpreds$ELVmax_WS,test_preds$ELVmax_WS, names=c("Ref","New"),ylab="ELVmax_WS (m)")
boxplot(refpreds$PMIN_WS,test_preds$PMIN_WS, names=c("Ref","New"),ylab="PMIN_WS (mm)")


#Resample bugs to 300 counts:
resamp_vec=rowSums(test_bugs_matrix)
resamp_vec[resamp_vec>300]=300
test_bugs_300cnt=rrarefy(test_bugs_matrix,resamp_vec)

OE.assess.test<-model.predict.RanFor.4.2(bugcal.pa,grps.final,preds.final, ranfor.mod=cluspred.rf,prednew=test_preds,bugnew=test_bugs_300cnt,Pc=0.5,Cal.OOB=FALSE)

boxplot(OE.assess.test$OE.scores$OoverE,ylab="New sample O/E scores")
plot(OE.assess.test$OE.scores$OoverE~resamp_vec,ylab="New sample O/E scores",xlab="Resample count")

write.csv(file="Z:\\buglab\\OE_Modeling\\NAMC_Supported_OEmodels\\AREMP2014\\ModelFiles\\InputsAndResults\\2021\\2142_OE_results.csv",cbind(OE.assess.test$OE.scores,resamp_vec))





##################################
########MMI

setwd("Z:\\buglab\\OE_Modeling\\NAMC_Supported_OEmodels\\AREMP2014\\ModelFiles")
load("AREMP_MMI.rdata")
library(randomForest)

##test_metrics_raw=read.csv(file="Z:\\buglab\\OE_Modeling\\NAMC_Supported_OEmodels\\AREMP2014\\InputsAndResults\\2015\\MMI_bugs_input.csv",row.names="SampleID")
##test_preds=read.csv(file="Z:\\buglab\\OE_Modeling\\NAMC_Supported_OEmodels\\AREMP2014\\InputsAndResults\\2015\\habitat_input.csv",row.names="SampleID")

test_metrics_raw=read.csv(file="Z:\\buglab\\OE_Modeling\\NAMC_Supported_OEmodels\\AREMP2014\\ModelFiles\\InputsAndResults\\2021\\Test_metrics_raw.csv",row.names="SampleID")
test_preds=read.csv(file="Z:\\buglab\\OE_Modeling\\NAMC_Supported_OEmodels\\AREMP2014\\ModelFiles\\InputsAndResults\\2021\\Test_preds.csv",row.names="SampleID")

#test_preds=read.csv(file="test_preds.csv",row.names="SampleID")
#test_metrics=read.csv(file="test_metrics_raw.csv",row.names="SampleID")

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


#Making metric predictions for new sites
test_metrics_prd=matrix(ncol=0,nrow=dim(test_metrics_raw)[1])
for(n in 1:length(rf_models)){
  model=get(rf_models[n])
  if(model$rsq[1500]>=0.1){
    metric_prd=predict(model,newdata=test_preds)}
  if(model$rsq[1500]<0.1){
    metric_prd=rep(0,times=dim(test_metrics_raw)[1])}
  test_metrics_prd=cbind(test_metrics_prd,metric_prd)
}
colnames(test_metrics_prd)=colnames(test_metrics_raw)

test_metrics_adj=test_metrics_raw-test_metrics_prd
test_metrics_rs=metricMatrixRescale(test_metrics_adj,ref_metrics_adj,mdeg_metrics_adj_cal)

test_MMI_scores=rowSums(test_metrics_rs)/6
boxplot(ref_scores,test_MMI_scores,mdeg_cal_scores,mdeg_val_scores,ylab="MMI score",names=c("Ref","Test","Mdeg cal","Mdeg val"))
abline(h=thresholds[1],lwd=3,col="red")
abline(h=thresholds[2],lwd=3,col="blue")


write.csv(file="Z:\\buglab\\OE_Modeling\\NAMC_Supported_OEmodels\\AREMP2014\\ModelFiles\\InputsAndResults\\2021\\2142_MMI_scores.csv",test_MMI_scores)
