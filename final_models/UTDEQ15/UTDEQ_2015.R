
UTDEQ_2015_model<-function(bugsOTU_matrix,test_preds,rf.model){

OE<-model.predict.RanFor.4.2(bugcal.pa,grps.final, preds.final,ranfor.mod=rfmod.final, prednew=test_preds,bugnew=bugsOTU_matrix,Pc=0.5,Cal.OOB=FALSE)

return(OE)
}
