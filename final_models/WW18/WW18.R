# Welcome to the westwide non-midge 2017 OE model script.
WW_model<-function(bugsOTU_matrix,test_preds,rf_model){
OE<-model.predict.RanFor.4.2(bugcal.pa[,c(2:222)],grps.final,preds.final,ranfor.mod=rf.mod.best.from.rFVR, prednew=test_preds, bugnew=bugsOTU_matrix, Pc=0.5, Cal.OOB=FALSE);
return(OE)
}
