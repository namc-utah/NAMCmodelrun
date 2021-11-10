# R code to run WYRIVPACS2012; 
# Provided by Eric Hargett, WYDEQ; 
# Edited by Christian Perry 20180312
# Edited by TWA February 2019

WY_model<-function(test_bugs,test_preds,rf.model){
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


#makes predictions for test data;
OE<-model.predict.v4.1(bugcal.pa,grps.final,preds.final, grpmns,covpinv,prednew=test_preds,bugnew=test_bugs,Pc=0.5)
return(OE)
}