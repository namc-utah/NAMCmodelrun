# R code to run WYRIVPACS2012; 
# Provided by Eric Hargett, WYDEQ; 
# Edited by Christian Perry 20180312
# Edited by TWA February 2019



setwd("\\\\share1.bluezone.usu.edu\\miller\\\\buglab\\OE_Modeling\\NAMC_Supported_OEmodels\\Wyoming\\Model Files\\")
library(gtools); library(MASS); library(cluster); library(Hmisc);

# Load necessary model and .Rdata files:
source("model.predict.v4.1.r");
load('WYRIVPACS2012.Rdata');

# Load reference bug and pred data:
predall<-read.table("WY_REFPRED.txt",row.names="SAMPLE",header=T,sep="\t");
bugall<-read.table("WY_REFBUG.txt",row.names="SAMPLE",header=T,sep="\t");
bugall<-bugall[,-1]; # delete first column of bugall;

#check sample(row) alignment of bug and predictor data;
row.names(bugall)==row.names(predall);
#If samples are not aligned, fix by aligning bugs data to predictor data, since latter is sorted by sample type;
bugall<-bugall[row.names(predall),];
row.names(bugall)==row.names(predall);

#Create a Presence/absence (1/0) matrix (site by taxa) for the bugs;
bugall.pa<-bugall;
bugall.pa[bugall.pa>0]<-1;

# Separate calibration ("C") and validation ("V") sites;
predcal<-predall[predall[,'TYPE']=='C',];  #pred - calibration sites;
pred.vld<-predall[substr(as.character(predall[,'TYPE']),1,1)=='V',];  #pred - validation sites;

bugcal<-bugall[predall[,'TYPE']=='C',]; #Bug Abundance matrix, calibration sites;
bug.vld<-bugall[substr(as.character(predall[,'TYPE']),1,1)=='V',];  #Bug Abundance matrix, validation sites;

bugcal.pa<-bugall.pa[predall[,'TYPE']=='C',]; #Bug presence/absence matrix, calibration sites;
bug.vld.pa<-bugall.pa[substr(as.character(predall[,'TYPE']),1,1)=='V',];  #Bug presence/absence matrix, validation sites;

#Data sets complete and aligned;
########################################;

# Load test data:

Boxdata <- ("\\\\share1.bluezone.usu.edu\\miller\\\\buglab\\OE_Modeling\\NAMC_Supported_OEmodels\\Wyoming\\InputsAndResults\\CurrentRun")
setwd(Boxdata)
getwd()
pred.test <- read.csv("WY_Habitat.csv",row.names="SampleID",header=T)
bug.test <- read.csv("WY_Bugs.csv",row.names="SampleID",header=T)
bug.test[bug.test>0]<-1; #convert bug test matrix to presence/absence;


#Drop all samples/sites that do not not have complete data for the model predictors;
pred.test<-pred.test[complete.cases(pred.test[,preds.final]),];
bug.test <-bug.test[row.names(pred.test),];

#makes predictions for test data;
OE.assess.test<-model.predict.v4.1(bugcal.pa,grps.final,preds.final, grpmns,covpinv,prednew=pred.test,bugnew=bug.test,Pc=0.5);
OE.assess.test$OE.scores;
write.table(OE.assess.test$OE.scores,'WY_OE_BLM2110.csv',sep=",",col.names=NA)


######### This section of code oututs the siteXtaxon probability of capture matrix.
# Example predictions: For nonreference sites in the Wyoming DEQ data set that are labeled "T" (see Step 1);

pred.test<-read.table("WY_REFPRED.txt",row.names="SAMPLE",header=T,sep="\t");  #predictor data-test sites;

#pred.test<-predall[as.character(predall[,'Mdl704'])=='N_lc',];  #predictor data - test sites;

bug.test.pa<-read.table("WY_REFBUG.txt",row.names="SAMPLE",header=T,sep="\t"); #load bug test matrix;
bug.test.pa[bug.test.pa>0]<-1; #convert bug test matrix to presence/absence;
#write.table(bug.test.pa,'WYbugtestpa.csv',sep=",",col.names=NA);

#bug.test.pa<-bugall.pa[as.character(predall[,'Mdl704'])=='N_lc',]; #Bug presence/absence matrix, test sites;

#Drop all samples/sites that do not not have complete data for the model predictors;
pred.test<-pred.test[complete.cases(pred.test[,preds.final]),];
bug.test.pa<-bug.test.pa[row.names(pred.test),];

#makes predictions for test data;
OE.assess.test<-model.predict.v4.1(bugcal.pa,grps.final,preds.final, grpmns,covpinv,prednew=pred.test,bugnew=bug.test.pa,Pc=0.5);

# look at O/E and BC scores of test-data samples;
OE.assess.test$OE.scores;
write.table(OE.assess.test,'testjul9refzero.csv',sep=",",col.names=NA)???




