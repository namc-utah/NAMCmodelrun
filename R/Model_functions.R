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



#' AZ perennial MMI
#'
#' @param bugnew
#' @param prednew
#'
#' @return MMI score
#' @export
#'
#' @examples
AZ_perennial_MMI_model<-function(bugnew,prednew){

  # MATS IBI ----

  ## Warm Water ----
  # Store reference numbers.  From the 2006 QAPP.  Appendix A. (also in 2021 QAPP)
  R.W.NumTaxa <- 37
  R.W.NumEphTaxa <- 9
  R.W.NumTriTaxa <- 9
  R.W.NumScraperTaxa <- 7
  R.W.PctScraper <- 23.7
  R.W.HBI <- 4.89
  R.W.NumDipTaxa <- 10
  R.W.PctEphm <- 70
  R.W.PctDom <- 19.1

  ### A. Metrics ----

  #### 1 - Total Taxa ----
  m.w.taxa <- bugnew %>%
    filter(is.na(exclude)) %>% # Exclude is for where OTU_ADEQ is NA due to taxonomist not getting low enough AND if that bug also identified to lower level
    group_by(SampleID) %>%
    summarize(NumTaxa = n()) %>%
    mutate(M.NumTaxa = ifelse(NumTaxa >= R.W.NumTaxa, 100, (NumTaxa/R.W.NumTaxa)* 100))

  #### 2 - Number of Ephemeroptera Taxa ----
  m.w.ephtaxa <- bugnew %>%
    filter(is.na(exclude)) %>%
    filter(Order == "Ephemeroptera") %>%
    group_by(SampleID) %>%
    summarize(NumEphTaxa = n()) %>%
    select(SampleID, NumEphTaxa) %>%
    mutate(M.NumEphTaxa = ifelse(NumEphTaxa >= R.W.NumEphTaxa, 100, (NumEphTaxa/R.W.NumEphTaxa)* 100))

  #### 3 - Number of Trichoptera Taxa ----
  m.w.tritaxa <- bugnew %>%
    filter(is.na(exclude)) %>%
    filter(Order == "Trichoptera") %>%
    group_by(SampleID) %>%
    summarize(NumTriTaxa = n()) %>%
    select(SampleID, NumTriTaxa) %>%
    mutate(M.NumTriTaxa = ifelse(NumTriTaxa >= R.W.NumTriTaxa, 100, (NumTriTaxa/R.W.NumTriTaxa)* 100))

  #### 4 - Number of Scraper Taxa ----
  m.w.scrapertaxa <- bugnew %>%
    filter(is.na(exclude)) %>%
    filter(Family != "Chironomidae") %>% # Family FFG for all Chironomidae is "UN" or null in master taxa
    filter(FFG == "Scraper") %>%
    group_by(SampleID) %>%
    summarize(NumScraperTaxa = n()) %>%
    select(SampleID, NumScraperTaxa) %>%
    mutate(M.NumScraperTaxa = ifelse(NumScraperTaxa >= R.W.NumScraperTaxa, 100, (NumScraperTaxa/R.W.NumScraperTaxa)* 100))

  #### 5 - Percent Scrapers ----
  m.w.pctscraper <- bugnew %>%
    mutate(FFG = ifelse(Family == "Chironomidae", NA, FFG)) %>% # numbers are not filtered out as in scraper taxa
    group_by(SampleID, FFG) %>%
    summarise(ffg_num = sum(Individuals)) %>%
    group_by(SampleID) %>%
    mutate(total_ind = sum(ffg_num)) %>%
    filter(FFG == "Scraper") %>%
    mutate(PctScraper = (ffg_num/total_ind) * 100) %>%
    select(SampleID, PctScraper) %>%
    mutate(M.PctScraper = ifelse(PctScraper >= R.W.PctScraper, 100, (PctScraper/R.W.PctScraper)* 100))

  #### 6 - Hilsenhoff Biotic Index ----
  m.w.hbi <- bugnew %>%
    filter(!is.na(TolVal)) %>% # no na's
    group_by(SampleID) %>%
    summarise(HBI = sum(TolVal * Individuals)/sum(Individuals)) %>% # HBI equation
    mutate(M.HBI = ifelse(HBI <= R.W.HBI, 100, (10 - HBI)/(10 - R.W.HBI)*100))

  #### 7 - Number of Diptera Taxa ----
  m.w.diptaxa <- bugnew %>%
    filter(is.na(exclude)) %>%
    filter(Order == "Diptera") %>%
    group_by(SampleID) %>%
    summarize(NumDipTaxa = n()) %>%
    select(SampleID, NumDipTaxa) %>%
    mutate(M.NumDipTaxa = ifelse(NumDipTaxa >= R.W.NumDipTaxa, 100, (NumDipTaxa/R.W.NumDipTaxa)* 100))

  #### 8 - Percent Ephemeroptera ----
  m.w.pctephm <- bugnew %>%
    group_by(SampleID, Order) %>%
    summarise(ephmnum = sum(Individuals)) %>%
    group_by(SampleID) %>%
    mutate(PctEphm = (ephmnum/sum(ephmnum))*100) %>%
    filter(Order == "Ephemeroptera") %>%
    select(SampleID, PctEphm) %>%
    mutate(M.PctEphm = ifelse(PctEphm >= R.W.PctEphm, 100, (PctEphm/R.W.PctEphm)* 100))

  #### 9 - Percent individuals in dominant taxon ----

  # Identify Dominant Taxon with the most individuals.
  m.w.domtaxa <- bugnew %>%
    arrange(desc(Individuals)) %>%
    group_by(SampleID) %>%
    mutate(PctDom = (Individuals/sum(Individuals))*100) %>%
    slice(1) %>% # slice picks x number of rows in the group to return. Here 1 represents the dominant taxa
    select(SampleID, PctDom) %>%
    mutate(M.PctDom = ifelse(PctDom <= R.W.PctDom, 100, (100 - PctDom)/(100 - R.W.PctDom)* 100))

  ### B. Index Normalization ----

  # Join metrics together
  m.w.all <- m.w.hbi %>% #6
    full_join(m.w.diptaxa, by = c("SampleID")) %>% #7
    full_join(m.w.domtaxa, by = c("SampleID")) %>%
    full_join(m.w.taxa, by = c("SampleID")) %>%  #1
    full_join(m.w.ephtaxa, by = c("SampleID")) %>% #2
    full_join(m.w.tritaxa, by = c("SampleID")) %>% #3
    full_join(m.w.scrapertaxa, by = c("SampleID")) %>% #4
    full_join(m.w.pctscraper, by = c("SampleID")) %>% #5
    full_join(m.w.pctephm, by = c("SampleID")) #8

  ### C. Index Calculation ----
  # Straight up average of normalized scores

  #NA in some final scores due to NAs in metrics?
  #find out why
  ibi.w <- m.w.all %>%
    left_join(prednew, by = "SampleID") %>%
    filter(ELEV_SITE <= 1524) %>%
    #Num Dip Taxa not being coerced to 0, even though it is metric used in both C and W models...?
    replace_na(list(M.HBI = 0, M.NumChironTaxa = 0, M.PctDom = 0, M.NumTaxa = 0, M.NumEphTaxa = 0, M.NumTriTaxa = 0, M.NumScraperTaxa = 0, M.PctScraper = 0, M.PctEphm = 0)) %>%
    mutate(modelId=236,IBI = (M.HBI + M.NumDipTaxa + M.PctDom + M.NumTaxa + M.NumEphTaxa + M.NumTriTaxa + M.NumScraperTaxa + M.PctScraper + M.PctEphm)/9) %>%
    select(SampleID, NumTaxa, NumTriTaxa, NumEphTaxa, NumDipTaxa, NumScraperTaxa, PctScraper, PctEphm, PctDom, HBI, everything())



  ##### seperate MMI ######

  ## Cold Water ----
  # Store reference numbers.  From the 2006 QAPP.  Appendix A. (also in 2021 QAPP)
  R.C.NumTaxa <- 38
  R.C.NumDipTaxa <- 11
  R.C.NumIntolTaxa <- 6
  R.C.NumScraperTaxa <- 11
  R.C.PctScraper <- 45.1
  R.C.PctPlec <- 19.1
  R.C.HBI <- 4.23

  ### A. Metrics ----

  #### 1 - Total Taxa ----
  m.c.taxa <- bugnew %>%
    filter(is.na(exclude)) %>%
    group_by(SampleID) %>%
    summarize(NumTaxa = n()) %>%
    mutate(M.NumTaxa = ifelse(NumTaxa >= R.C.NumTaxa, 100, (NumTaxa/R.C.NumTaxa)* 100))

  #### 2 - Number of Diptera Taxa ----
  # ADEQ method from QAPP.  Uses genus for all insects but midges which go to family (so family just counted as one taxa)
  m.c.diptaxa <- bugnew %>%
    filter(is.na(exclude)) %>%
    filter(Order == "Diptera") %>%
    group_by(SampleID) %>%
    summarize(NumDipTaxa = n()) %>%
    select(SampleID, NumDipTaxa) %>%
    mutate(M.NumDipTaxa = ifelse(NumDipTaxa >= R.C.NumDipTaxa, 100, (NumDipTaxa/R.C.NumDipTaxa)* 100))

  #### 3 - Number of Intolerant Taxa ----
  # Number of taxa having a tolerance value of <= 3
  m.c.intoltaxa <- bugnew %>%
    filter(is.na(exclude)) %>%
    filter(TolVal <= 3) %>% # no na's
    group_by(SampleID) %>%
    summarise(NumIntolTaxa = n()) %>%
    mutate(M.NumIntolTaxa = ifelse(NumIntolTaxa >= R.C.NumIntolTaxa, 100, (NumIntolTaxa/R.C.NumIntolTaxa)* 100))

  #### 4 - Number of Scraper Taxa ----
  m.c.scrapertaxa <- bugnew %>%
    filter(is.na(exclude)) %>%
    filter(Family != "Chironomidae") %>% # Family FFG for all Chironomidae is "UN" or null in master taxa
    filter(FFG == "Scraper") %>%
    group_by(SampleID) %>%
    summarize(NumScraperTaxa = n()) %>%
    select(SampleID, NumScraperTaxa) %>%
    mutate(M.NumScraperTaxa = ifelse(NumScraperTaxa >= R.C.NumScraperTaxa, 100, (NumScraperTaxa/R.C.NumScraperTaxa)* 100))

  #### 5 - Percent Scrapers ----
  m.c.pctscraper <- bugnew %>%
    mutate(FFG = ifelse(Family == "Chironomidae", NA, FFG)) %>% # Abundances numbers are not filtered out like scraper taxa
    group_by(SampleID, FFG) %>%
    summarise(ffg_num = sum(Individuals)) %>%
    group_by(SampleID) %>%
    mutate(total_ind = sum(ffg_num)) %>%
    filter(FFG == "Scraper") %>%
    mutate(PctScraper = (ffg_num/total_ind) * 100) %>%
    select(SampleID, PctScraper) %>%
    mutate(M.PctScraper = ifelse(PctScraper >= R.C.PctScraper, 100, (PctScraper/R.C.PctScraper)* 100))

  #### 6 - Percent Plecoptera ----
  m.c.pctplec <- bugnew %>%
    group_by(SampleID, Order) %>%
    summarise(plecnum = sum(Individuals)) %>%
    group_by(SampleID) %>%
    mutate(PctPlec = (plecnum/sum(plecnum))*100) %>%
    filter(Order == "Plecoptera") %>%
    select(SampleID, PctPlec) %>%
    mutate(M.PctPlec = ifelse(PctPlec >= R.C.PctPlec, 100, (PctPlec/R.C.PctPlec)* 100))

  #### 7 - Hilsenhoff Biotic Index ----
  m.c.hbi <- bugnew %>%
    filter(!is.na(TolVal)) %>% # no na's
    group_by(SampleID) %>%
    summarise(HBI = sum(TolVal * Individuals)/sum(Individuals)) %>%
    mutate(M.HBI = ifelse(HBI <= R.C.HBI, 100, (10 - HBI)/(10 - R.C.HBI)*100))

  ### B. Index Normalization ----

  # Join metrics together
  m.c.all <- m.c.hbi %>% #7
    full_join(m.c.taxa, by = c("SampleID")) %>%  #1
    full_join(m.c.diptaxa, by = c("SampleID")) %>% #2
    full_join(m.c.intoltaxa, by = c("SampleID")) %>% #3
    full_join(m.c.scrapertaxa, by = c("SampleID")) %>% #4
    full_join(m.c.pctscraper, by = c("SampleID")) %>% #5
    full_join(m.c.pctplec, by = c("SampleID")) #6

  ### C. Index Calculation ----
  # Straight up average of normalized scores
  ibi.c <- m.c.all %>%
    left_join(prednew, by = "SampleID") %>%
    filter(ELEV_SITE>1524) %>%
    replace_na(list(M.HBI = 0, M.NumTaxa = 0, M.DipTaxa = 0, M.NumIntolTaxa = 0, M.NumScraperTaxa = 0, M.PctScraper = 0, M.PctPlec = 0)) %>%
    mutate(modelId=169,IBI = (M.HBI + M.NumTaxa + M.NumDipTaxa + M.NumIntolTaxa + M.NumScraperTaxa + M.PctScraper + M.PctPlec)/7) %>%
    select(StationID, CollDate, NumTaxa, NumDipTaxa, NumIntolTaxa, NumScraperTaxa, PctScraper, PctPlec, HBI, everything())

  ## IBI ----
  ibi <- ibi.w %>%
    bind_rows(ibi.c)
  #add sampleId to the final table

  return(ibi)
}


# AZ_MMI_model<-function(bugnew,prednew){
#
#   thresholds<-textConnection("metricId	metric	warm	cold
# 351	Totaltaxa	37	38
# 354	  Trichopterataxa	9	 NA
# 352	  Ephemeropterataxa	9 NA
# 356	  Dipterataxa	10	11
# 375	  Scrapertaxa	7	11
# 199	  Percentscraper	23.7	45.1
# 176	  PercentEphemeroptera	70  NA
# 448	  PercentDominantTaxon	19.1	NA
# 444	  HilsenhoffBioticIndex	4.89	4.21
# 177	percentplecoptera		19.1  NA
# 381	intoleranttaxa		6 NA
# ")
#   metricThresholds <- read.table(thresholds, header = TRUE, stringsAsFactors = FALSE)
#
#  j=dplyr::left_join(bugnew,metricThresholds, by='metricId')
#
#  decreasers=subset(j,metricId %in% c(351,354,352,356,375,199,176,177,381))
#  increasers=subset(j,metricId %in% c(444,448))
#  decreasers$coldmetric=(decreasers$metricValue/decreasers$cold)*100
#  decreasers$warmmetric=(decreasers$metricValue/decreasers$warm)*100
#  increasers$coldmetric=ifelse(increasers$metricId==444,(10-increasers$metricValue)/(10-increasers$cold)*100,NA)
#  increasers$warmmetric=ifelse(increasers$metricId==444,(10-increasers$metricValue)/(10-increasers$warm)*100,NA)
#  increasers$warmmetric=ifelse(increasers$metricId==448,(100-increasers$metricValue)/(100-increasers$warm)*100,increasers$warmmetric)
#  metrics=rbind(decreasers,increasers)
#  #join in prednew to modelResults
#  metrics=dplyr::left_join(metrics,prednew, by="sampleId")
#  colddf=subset(metrics,metrics$ElevCat>1524)
#  warmdf=subset(metrics,metrics$ElevCat<1524)
#  coldm=setNames(stats::aggregate(coldmetric~sampleId,data=colddf,FUN=mean,na.rm=TRUE),c("sampleId","MMI"))
#  warmm=setNames(stats::aggregate(warmmetric~sampleId,data=warmdf,FUN=mean,na.rm=TRUE),c("sampleId","MMI"))
#  modelResults=rbind(coldm,warmm)
# }

rarify<-function(inbug, sample.ID, abund, subsiz){
  start.time=proc.time();
  outbug<-inbug;
  sampid<-unique(inbug[,sample.ID]);
  nsamp<-length(sampid);
  #parameters are set up;
  #zero out all abundances in output data set;
  outbug[,abund]<-0;
  #loop over samples, rarify each one in turn;

  for(i in 1:nsamp) { ;
    #extract current sample;
    isamp<-sampid[i];
    flush.console();
    print(as.character(isamp));
    onesamp<-inbug[inbug[,sample.ID]==isamp,];
    onesamp<-data.frame(onesamp,row.id=seq(1,dim(onesamp)[[1]])); #add sequence numbers as a new column;
    #expand the sample into a vector of individuals;
    samp.expand<-rep(x=onesamp$row.id,times=onesamp[,abund]);
    nbug<-length(samp.expand); #number of bugs in sample;
    #vector of uniform random numbers;
    ranvec<-runif(n=nbug);
    #sort the expanded sample randomly;
    samp.ex2<-samp.expand[order(ranvec)];
    #keep only the first piece of ranvec, of the desired fised count size;
    #if there are fewer bugs than the fixed count size, keep them all;
    if(nbug>subsiz){subsamp<-samp.ex2[1:subsiz]} else{subsamp<-samp.ex2};
    #tabulate bugs in subsample;
    subcnt<-table(subsamp);
    #define new subsample frame and fill it with new reduced counts;
    newsamp<-onesamp;
    newsamp[,abund]<-0;
    newsamp[match(newsamp$row.id,names(subcnt),nomatch=0)>0,abund]<-as.vector(subcnt);
    outbug[outbug[,sample.ID]==isamp,abund]<-newsamp[,abund];
  }; #end of sample loop;

  elaps<-proc.time()-start.time;
  print(c("Rarefaction complete. Number of samples = ",nsamp),quote=F);
  print(c("Execution time (sec)= ", elaps[1]),quote=F);
  outbug; #return subsampled data set as function value;
} #end of function;


#' OR MMI
#'
#' @param bugnew
#' @param prednew
#' @param rfmod_nt_habitat_rheo
#' @param rfmod_pi_EPTNoHydro
#' @param rfmod_pi_ti_stenocold_cold_cool
#' @param rfmod_pt_tv_intol
#'
#' @return MMI score
#' @export
#'
#' @examples
OR_MMI_model<-function(bugnew,prednew,rfmod_nt_habitat_rheo,rfmod_pi_EPTNoHydro,rfmod_pi_ti_stenocold_cold_cool,rfmod_pt_tv_intol){

  #join bug data to predictors
  Drfdat=dplyr::left_join(bugnew,prednew,by='SampleID')

  # which rf models to use
  rfmodels <- c('rfmod_pt_tv_intol', 'rfmod_nt_habitat_rheo', 'rfmod_pt_ti_stenocold_cold_cool',
                'rfmod_pi_EPTNoHydro')
  ## test site predictions -------------------------------------------------------------------------------------------
  Dpredictions=list()
  for (i in 1:length(rfmodels)){
    tryCatch({Dpredictions[[paste0("E.",rfmodels[i])]]<- round(predict(eval(parse(text =paste0(rfmodels[i]))), Drfdat, type = "response"),digits=4)
    }, error =function (e){
      cat(paste0("/n/tERROR calculating: ",paste0(names(Drfdat)[i],"_pred"),"/n"))
      str(e,indent.str = "   "); cat("/n")
    })
  }
  predictionsdf=as.data.frame(do.call(cbind,Dpredictions))
  #join predictions into master dataframe
  Drfdat2=cbind(Drfdat,predictionsdf)

  ## CALCULATE RESIDUALS ---------------------------------------------------------------------------------------------
  resid=list()
  for (i in 2:5){
    tryCatch({resid[[paste0(colnames(Drfdat2)[i],"_resid")]]=Drfdat2[,i]- Drfdat2[,paste0("E.rfmod_",colnames(Drfdat2)[i])]

    }, error =function (e){
      cat(paste0("/n/tERROR calculating: ",paste0(Drfdat2[i],"_resid"),"/n"))
      str(e,indent.str = "   "); cat("/n")
    })

  }
  residualsdf=as.data.frame(do.call(cbind,resid))
  rfdat_all_final4=cbind(Drfdat2,residualsdf)


  # Select SAMPLEID + 4 metrics.residuals
  candmetrics <- rfdat_all_final4 %>%
    select(SAMPLEID, pt_tv_intol_resid, nt_habitat_rheo_resid,
           pt_ti_stenocold_cold_cool_resid, pi_EPTNoHydro_resid)

  metric_rs <- candmetrics |>
    mutate(pt_tv_intol_resid = (pt_tv_intol_resid - -50.54222) /  (14.13465 -  -50.54222) ,
           nt_habitat_rheo_resid = (nt_habitat_rheo_resid - -24.334100) /  (8.299599 -  -24.334100) ,
           pt_ti_stenocold_cold_cool_resid = (pt_ti_stenocold_cold_cool_resid - -38.35066) /  (15.15093 -  -38.35066),
           pi_EPTNoHydro_resid = (pi_EPTNoHydro_resid - -54.76710) /  (24.54683 -  -54.76710)) %>%
    mutate(pt_tv_intol_resid = case_when(pt_tv_intol_resid <0 ~ 0,
                                         pt_tv_intol_resid >1 ~ 1,
                                         TRUE ~ pt_tv_intol_resid)) %>%
    mutate(nt_habitat_rheo_resid = case_when(nt_habitat_rheo_resid <0 ~ 0,
                                             nt_habitat_rheo_resid >1 ~ 1,
                                             TRUE ~ nt_habitat_rheo_resid)) %>%
    mutate(pt_ti_stenocold_cold_cool_resid = case_when(pt_ti_stenocold_cold_cool_resid <0 ~ 0,
                                                       pt_ti_stenocold_cold_cool_resid >1 ~ 1,
                                                       TRUE ~ pt_ti_stenocold_cold_cool_resid)) %>%
    mutate( pi_EPTNoHydro_resid = case_when( pi_EPTNoHydro_resid <0 ~ 0,
                                             pi_EPTNoHydro_resid >1 ~ 1,
                                             TRUE ~  pi_EPTNoHydro_resid)) %>%
    mutate(MMI = (pt_tv_intol_resid + nt_habitat_rheo_resid +
                    pt_ti_stenocold_cold_cool_resid+pi_EPTNoHydro_resid)  /4)

  #Throw MMI result and MMI metrics into a list for getting out of the function
  metric_rs=dplyr::left_join(bugnew, metric_rs, by="SampleID")
  return(metric_rs)
}
