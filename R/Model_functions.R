# vansickles 3 functions are in seperate files


#' Oregon NBR null model
#'
#' @param bugnew
#'
#' @return OoverE
#' @export
#'
#' @examples
OR_NBR_model<-function(bugnew){
  ref_taxa <- c(
  'Baetis',
  'Brachycentrus',
  'Chironominae',
  'Diphetor_hageni',
  'Epeorus',
  'Optioservus',
  'Orthocladiinae',
  'Rhyacophila',
  'Trombidiformes',
  'Zaitzevia'
) #list of 10 "reference" taxa

bugnew$TaxaToCount <- ifelse(bugnew$otuName %in% ref_taxa, 1, 0)

OE <-bugnew %>% dplyr::group_by(sampleId) %>% dplyr::summarize(O = sum(TaxaToCount))
OE$E <-
  7.56 #null E is always the same. Does not account for environmental factors
OE$OoverE <- OE$O / OE$E

return(OE)
}





# MMIs
# each MMI will need its own function because of the variable name and number of randomforest model inputs....

###this function requires yet another function metricMatrixRescale which pulls in additional things...
#need to think about how to turn this into a function or if it should be a package and if this is a more sustainable way of writing MMIs.
#it is more efficient but not as transparent as NV code
#' AREMP MMI
#'
#' @param bugnew
#' @param prednew
#' @param CLING_rich.rf
#' @param DIPT_rich.rf
#' @param LLT_rich.rf
#' @param NON_INSECT_rich.rf
#' @param PER_EPT.rf
#' @param PER_INTOL.rf
#' @param rf_models
#' @param mdeg_metrics_adj_cal
#' @param ref_metrics_adj
#'
#' @return MMI value
#' @export
#'
#' @examples
AREMP_MMI_model<-function(bugnew,prednew,CLING_rich.rf,DIPT_rich.rf,LLT_rich.rf,NON_INSECT_rich.rf,PER_EPT.rf,PER_INTOL.rf,rf_models,mdeg_metrics_adj_cal,ref_metrics_adj){
  bugnew_prd=matrix(ncol=0,nrow=dim(bugnew)[1])
  for(n in 1:length(rf_models)){
    model=get(rf_models[n])
    if(model$rsq[1500]>=0.1){
      metric_prd=predict(model,prednew)}
    if(model$rsq[1500]<0.1){
      metric_prd=rep(0,times=dim(bugnew)[1])}
    bugnew_prd=cbind(bugnew_prd,metric_prd)
  }
  colnames(bugnew_prd)=colnames(bugnew)

  bugnew_adj=bugnew-bugnew_prd
  bugnew_rs=metricMatrixRescale(metrics=bugnew_adj,ref_metrics=ref_metrics_adj,mostdeg_metrics=mdeg_metrics_adj_cal)
  bugnew_rs=as.data.frame(bugnew_rs)
  bugnew_rs$MMI=rowSums(bugnew_rs)/6
  return(bugnew_rs)
}


#' MMI metric rescaling
#'
#' @param metrics
#' @param ref_metrics
#' @param mostdeg_metrics
#'
#' @return rescaled standardized metrics to min and max of reference and most degraded data
#' @export
#'
#' @examples
metricMatrixRescale<-function(metrics,ref_metrics,mostdeg_metrics){
  if(any(colnames(metrics)==colnames(ref_metrics))==FALSE){stop("Columns in new metrics must match columns in ref metrics")}
  if(any(colnames(metrics)==colnames(mostdeg_metrics))==FALSE){stop("Columns in new metrics must match columns in most deg metrics")}
  if(any(colnames(mostdeg_metrics)==colnames(ref_metrics))==FALSE){stop("Columns in ref metrics must match columns in most deg metrics")}
  metrics_rs=matrix(nrow=dim(metrics)[1],ncol=0)
  for(n in 1:dim(metrics)[2]){
    metric=metrics[,n]
    ref_metric=ref_metrics[,n]
    mostdeg_metric=mostdeg_metrics[,n]
    if(mean(ref_metric)>mean(mostdeg_metric)){
      min=quantile(mostdeg_metric,0.05)
      max=quantile(ref_metric,0.95)
      metric_rs=(metric-min)/(max-min)}
    if(mean(ref_metric)<mean(mostdeg_metric)){
      min=quantile(ref_metric,0.05)
      max=quantile(mostdeg_metric,0.95)
      metric_rs=1-((metric-min)/(max-min))}
    metric_rs[metric_rs>1]=1
    metric_rs[metric_rs<0]=0
    metrics_rs=cbind(metrics_rs,metric_rs)}
  colnames(metrics_rs)=colnames(metrics)
  row.names(metrics_rs)=rownames(metrics)
  return(metrics_rs)}



#' NV MMI
#'
#' @param bugnew
#' @param prednew
#' @param CLINGER.rf
#' @param INSET.rf
#' @param NONSET.rf
#' @param PER_CFA.rf
#' @param PER_EPHEA.rf
#' @param PER_PLECA.rf
#'
#' @return MMI score
#' @export
#'
#' @examples
NV_MMI_model<-function(bugnew,prednew,CLINGER.rf,INSET.rf,NONSET.rf,PER_CFA.rf,PER_EPHEA.rf,PER_PLECA.rf){
  ####adjust metrics for natural variability and rescale
  INSET.raw=bugnew$INSET
  INSET.pred=predict(INSET.rf, prednew, type="response")
  INSET.adj=INSET.raw-INSET.pred
  INSET.rs=100*((INSET.adj--20.23451)/(9.53285--20.23451))
  INSET.rs[INSET.rs>100]=100
  INSET.rs[INSET.rs<0]=0

  PER_CFA.raw=bugnew$PER_CFA
  PER_CFA.pred=predict(PER_CFA.rf, prednew, type="response")
  PER_CFA.adj=PER_CFA.raw-PER_CFA.pred
  PER_CFA.rs=100*((PER_CFA.adj--26.84895)/(20.15404--26.84895))
  PER_CFA.rs[PER_CFA.rs>100]=100
  PER_CFA.rs[PER_CFA.rs<0]=0

  PER_EPHEA.raw=bugnew$PER_EPHEA
  PER_EPHEA.pred=predict(PER_EPHEA.rf, prednew, type="response")
  PER_EPHEA.adj=PER_EPHEA.raw-PER_EPHEA.pred
  PER_EPHEA.rs=100*((PER_EPHEA.adj--37.05147)/(29.21069--37.05147))
  PER_EPHEA.rs[PER_EPHEA.rs>100]=100
  PER_EPHEA.rs[PER_EPHEA.rs<0]=0

  NONSET.raw=bugnew$NONSET
  NONSET.pred=predict(NONSET.rf, prednew, type="response")
  NONSET.adj=NONSET.raw-NONSET.pred
  NONSET.rs=100*((NONSET.adj--3.733934)/(3.773724--3.733934))
  NONSET.rs[NONSET.rs>100]=100
  NONSET.rs[NONSET.rs<0]=0

  CLINGER.raw=bugnew$CLINGER
  CLINGER.pred=predict(CLINGER.rf, prednew, type="response")
  CLINGER.adj=CLINGER.raw-CLINGER.pred
  CLINGER.rs=100*((CLINGER.adj--8.310543)/(5.719022--8.310543))
  CLINGER.rs[CLINGER.rs>100]=100
  CLINGER.rs[CLINGER.rs<0]=0

  PER_PLECA.raw=bugnew$PER_PLECA
  PER_PLECA.pred=predict(PER_PLECA.rf, prednew, type="response")
  PER_PLECA.adj=PER_PLECA.raw-PER_PLECA.pred
  PER_PLECA.rs=100*((PER_PLECA.adj--6.212933)/(15.33518--6.212933))
  PER_PLECA.rs[PER_PLECA.rs>100]=100
  PER_PLECA.rs[PER_PLECA.rs<0]=0

  SHDIVER.raw=bugnew$SHDIVER #SHDIVER not adjusted by natural gradients
  SHDIVER.rs=100*((SHDIVER.raw-0.8980395)/(3.18087-0.8980395))
  SHDIVER.rs[SHDIVER.rs>100]=100
  SHDIVER.rs[SHDIVER.rs<0]=0

  #Calculate MMI scores for new samples
  bugnew.rs=data.frame(INSET.rs, PER_EPHEA.rs, SHDIVER.rs, PER_CFA.rs, PER_PLECA.rs, NONSET.rs, CLINGER.rs, row.names=row.names(bugnew))
  bugnew.rs$MMI=rowSums(bugnew.rs)/7
  return(bugnew.rs)
}

AZ_MMI_model<-function(bugnew,prednew){

  thresholds<-textConnection("metricId	metric	warm	cold
351	Totaltaxa	37	38
354	  Trichopterataxa	9	 NA
352	  Ephemeropterataxa	9 NA
356	  Dipterataxa	10	11
375	  Scrapertaxa	7	11
199	  Percentscraper	23.7	45.1
176	  PercentEphemeroptera	70  NA
448	  PercentDominantTaxon	19.1	NA
444	  HilsenhoffBioticIndex	4.89	4.21
177	percentplecoptera		19.1  NA
381	intoleranttaxa		6 NA
")
  metricThresholds <- read.table(thresholds, header = TRUE, stringsAsFactors = FALSE)

 j=dplyr::left_join(bugnew,metricThresholds, by='metricId')

 decreasers=subset(j,metricId %in% c(351,354,352,356,375,199,176,177,381))
 increasers=subset(j,metricId %in% c(444,448))
 decreasers$coldmetric=(decreasers$metricValue/decreasers$cold)*100
 decreasers$warmmetric=(decreasers$metricValue/decreasers$warm)*100
 increasers$coldmetric=ifelse(increasers$metricId==444,(10-increasers$metricValue)/(10-increasers$cold)*100,NA)
 increasers$warmmetric=ifelse(increasers$metricId==444,(10-increasers$metricValue)/(10-increasers$warm)*100,NA)
 increasers$warmmetric=ifelse(increasers$metricId==448,(100-increasers$metricValue)/(100-increasers$warm)*100,increasers$warmmetric)
 metrics=rbind(decreasers,increasers)
 #join in prednew to modelResults
 metrics=dplyr::left_join(metrics,prednew, by="sampleId")
 colddf=subset(metrics,metrics$ElevCat>1524)
 warmdf=subset(metrics,metrics$ElevCat<1524)
 coldm=setNames(stats::aggregate(coldmetric~sampleId,data=colddf,FUN=mean,na.rm=TRUE),c("sampleId","MMI"))
 warmm=setNames(stats::aggregate(warmmetric~sampleId,data=warmdf,FUN=mean,na.rm=TRUE),c("sampleId","MMI"))
 modelResults=rbind(coldm,warmm)
}
