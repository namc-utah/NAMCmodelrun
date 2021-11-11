NV_MMI_model<-function(test_bugs_metrics,test_preds,rf_model){
  ####adjust metrics for natural variability and rescale
  INSET.raw=test_bugs_metrics$INSET
  INSET.pred=predict(INSET.rf, newdata=test_preds, type="response")
  INSET.adj=INSET.raw-INSET.pred
  INSET.rs=100*((INSET.adj--20.23451)/(9.53285--20.23451))
  INSET.rs[INSET.rs>100]=100
  INSET.rs[INSET.rs<0]=0
  
  PER_CFA.raw=test_bugs_metrics$PER_CFA
  PER_CFA.pred=predict(PER_CFA.rf, newdata=test_preds, type="response")
  PER_CFA.adj=PER_CFA.raw-PER_CFA.pred
  PER_CFA.rs=100*((PER_CFA.adj--26.84895)/(20.15404--26.84895))
  PER_CFA.rs[PER_CFA.rs>100]=100
  PER_CFA.rs[PER_CFA.rs<0]=0
  
  PER_EPHEA.raw=test_bugs_metrics$PER_EPHEA
  PER_EPHEA.pred=predict(PER_EPHEA.rf, newdata=test_preds, type="response")
  PER_EPHEA.adj=PER_EPHEA.raw-PER_EPHEA.pred
  PER_EPHEA.rs=100*((PER_EPHEA.adj--37.05147)/(29.21069--37.05147))
  PER_EPHEA.rs[PER_EPHEA.rs>100]=100
  PER_EPHEA.rs[PER_EPHEA.rs<0]=0
  
  NONSET.raw=test_bugs_metrics$NONSET
  NONSET.pred=predict(NONSET.rf, newdata=test_preds, type="response")
  NONSET.adj=NONSET.raw-NONSET.pred
  NONSET.rs=100*((NONSET.adj--3.733934)/(3.773724--3.733934))
  NONSET.rs[NONSET.rs>100]=100
  NONSET.rs[NONSET.rs<0]=0
  
  CLINGER.raw=test_bugs_metrics$CLINGER
  CLINGER.pred=predict(CLINGER.rf, newdata=test_preds, type="response")
  CLINGER.adj=CLINGER.raw-CLINGER.pred
  CLINGER.rs=100*((CLINGER.adj--8.310543)/(5.719022--8.310543))
  CLINGER.rs[CLINGER.rs>100]=100
  CLINGER.rs[CLINGER.rs<0]=0
  
  PER_PLECA.raw=test_bugs_metrics$PER_PLECA
  PER_PLECA.pred=predict(PER_PLECA.rf, newdata=test_preds, type="response")
  PER_PLECA.adj=PER_PLECA.raw-PER_PLECA.pred
  PER_PLECA.rs=100*((PER_PLECA.adj--6.212933)/(15.33518--6.212933))
  PER_PLECA.rs[PER_PLECA.rs>100]=100
  PER_PLECA.rs[PER_PLECA.rs<0]=0
  
  SHDIVER.raw=test_bugs_metrics$SHDIVER #SHDIVER not adjusted by natural gradients
  SHDIVER.rs=100*((SHDIVER.raw-0.8980395)/(3.18087-0.8980395))
  SHDIVER.rs[SHDIVER.rs>100]=100
  SHDIVER.rs[SHDIVER.rs<0]=0
  
  #Calculate MMI scores for new samples
  test_bugs_metrics.rs=data.frame(INSET.rs, PER_EPHEA.rs, SHDIVER.rs, PER_CFA.rs, PER_PLECA.rs, NONSET.rs, CLINGER.rs, row.names=row.names(test_bugs_metrics))
  MMI=rowSums(test_bugs_metrics.rs)/7
  return(MMI)
}

