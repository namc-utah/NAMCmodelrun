###Apply AREMP model to test sites:

#### version 4.2 function
AREMP_OE_model<-function(bugnew,prednew,ranfor.mod){
  OE<-model.predict.RanFor.4.2(bugcal.pa,grps.final,preds.final, ranfor.mod,prednew,bugnew,Pc=0.5,Cal.OOB=FALSE)
  return(OE)
}
UTDEQ_2015_model<-function(bugnew,prednew,ranfor.mod){
  OE<-model.predict.RanFor.4.2(bugcal.pa,grps.final, preds.final,ranfor.mod, prednew,bugnew,Pc=0.5,Cal.OOB=FALSE)
  return(OE)
}
WW_model<-function(bugnew,prednew,ranfor.mod){
  OE<-model.predict.RanFor.4.2(bugcal.pa[,c(2:222)],grps.final,preds.final,ranfor.mod, prednew, bugnew, Pc=0.5, Cal.OOB=FALSE);
  return(OE) # bugcal.pa data set needs reviewed and standardized. need to set sampleiDs and the row names but I dont know why we are excluding taxa
}


#### version 4.1 function
OR_MWCF_model<- function(bugnew,prednew,ranfor.mod){
  #make predictions for test data;
  OE<-model.predict.v4.1(bugcal.pa,grps.final,preds.final, grpmns,covpinv,prednew,bugnew,Pc=0.5);
  return(OE)
}
OR_WCCP_model<- function (bugnew,prednew,ranfor.mod){
  #make predictions for test data;
  OE<-model.predict.v4.1(bugcal.pa,grps.final,preds.final, grpmns,covpinv,prednew,bugnew,Pc=0.5);
  return(OE)
}
WY_model<-function(bugnew,prednew,ranfor.mod){
  OE<-model.predict.v4.1(bugcal.pa,grps.final,preds.final, grpmns,covpinv,prednew,bugnew,Pc=0.5)
  return(OE)
}


##### version 0 function
PIBO_model<-function (bugnew,prednew,ranfor.mod){
  # merges test bug data with a null set of all OTU bug names in order for VanSickle code to work properly.
  bugnew.pa <- merge(x=bug.otu, y=bugnew, all=TRUE)         
  # converts na data produced from columns above into 0s;
  bugnew.pa[is.na(bugnew.pa)]<- 0 
  
  # makes predictions for test data;
  OE<-model.predict.RanFor(bugall, grps.final, preds.final, ranfor.mod, prednew, bugnew.pa, Pc=0.5);
  
  ###is this needed????
  # Append Sampleid to results
  OE$OE.scores = data.frame(bugnew.pa$sampleId, OE$OE.scores)
  names(OE$OE.scores) = c('sampleId',names(OE$OE.scores)[2:ncol(OE$OE.scores)] )
  OE$Capture.Probs = data.frame(bugnew.pa$sampleId, OE$Capture.Probs)
  names(OE$Capture.Probs) = c('sampleId',names(OE$Capture.Probs)[2:ncol(OE$Capture.Probs)] )
  
  return(OE)
}