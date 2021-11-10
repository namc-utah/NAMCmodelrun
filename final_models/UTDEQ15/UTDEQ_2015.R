
UTDEQ_2015_model<-function(test_bugs,test_preds,rf.model){

OE<-model.predict.RanFor.4.2(bugcal.pa,grps.final, preds.final,ranfor.mod=rfmod.final, prednew=test_preds,bugnew=test_bugs,Pc=0.5,Cal.OOB=FALSE)

return(OE)
}