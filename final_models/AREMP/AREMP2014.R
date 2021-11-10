###Apply AREMP model to test sites:
AREMP_OE_model<-function(test_bugs,test_preds,rf_model){
  OE<-model.predict.RanFor.4.2(bugcal.pa,grps.final,preds.final, ranfor.mod=cluspred.rf,prednew=test_preds,bugnew=test_bugs,Pc=0.5,Cal.OOB=FALSE)
 return(OE)
}


##################################
########MMI
#Making metric predictions for new sites
AREMP_MMI_model<-function(test_bugs,test_preds,rf_model){
  test_bugs_prd=matrix(ncol=0,nrow=dim(test_bugs_raw)[1])
  for(n in 1:length(rf_models)){
    model=get(rf_models[n])
    if(model$rsq[1500]>=0.1){
      metric_prd=predict(model,newdata=test_preds)}
    if(model$rsq[1500]<0.1){
      metric_prd=rep(0,times=dim(test_bugs_raw)[1])}
    test_bugs_prd=cbind(test_bugs_prd,metric_prd)
  }
  colnames(test_bugs_prd)=colnames(test_bugs_raw)
  
  test_bugs_adj=test_bugs_raw-test_bugs_prd
  test_bugs_rs=metricMatrixRescale(test_bugs_adj,ref_metrics_adj,mdeg_metrics_adj_cal)
  
  MMI=rowSums(test_bugs_rs)/6
  return(MMI)
}
