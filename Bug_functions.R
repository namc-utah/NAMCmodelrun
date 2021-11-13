
OE_bug_matrix<-function(sampleId,translationId,fixedCount){
  bugsOTU = query(
    "sampleMetrics",
    translationId = def_models$translationId,
    fixedCount = def_models$fixedCount,
    sampleIds = def_model_results$sampleId
  )
  rarefiedOTUTaxa = subset(bugsOTU, metricName == "Rarefied Taxa")
  rarefiedOTUTaxa = NAMCr::json.expand(rarefiedOTUTaxa, "metricValue")
  sumrarefiedOTUTaxa = rarefiedOTUTaxa  %>%
    dplyr::group_by(sampleId, OTUName) %>%
    dplyr::summarize(sumSplitCount = sum(splitCount)) # why are multiple records exported here per OTU???
  
  sumrarefiedOTUTaxa$presence = ifelse(sumrarefiedOTUTaxa$sumSplitCount >=1, 1, 0)
  bugnew = tidyr::pivot_wider(sumrarefiedOTUTaxa,id_cols = "sampleId", names_from = "OTUName",values_from = "presence")
}


MMI_metrics<-function(sampleId,fixedCount){
  # need to add list of required metrics to model table as well as random forest input file
  #get needed list of metrics from new end point
  # NV and AREMP relative abundances were OTU standardized too!
  bugsMetrics = query("sampleMetrics",
                      def_models$fixedCount,
                      sampleIds = def_model_results$SampleId)
  modelInfo = query("modelInfo", modelId = def_model_results$modelId)
  #MMI_metrics = subset(bugsMetrics, metricId %in% ()) # store translations to metric names in database
  # need to replace metric name with a model specific metric abbreviation
  bugnew = tidyr::pivot_wider(MMI_metrics,id_cols = "sampleId",names_from = "metricName",values_from = "metricValue")
}




CSCI_bug <- function(sampleId, translationId){
  
  bugsTranslation = query(
    "sampleTaxaTranslation",
    translationId = def_models$translationId,
    sampleIds = def_model_results$SampleId
  )
  sites = query(
    "samples",
    include = c("sampleId", 'siteId'),
    sampleIds = def_model_results$SampleId
  )
  bugRaw = query(
    "sampleRawTaxa",
    include = c(),
    translationId = def_models$translationId,
    sampleIds = def_model_results$SampleId
  )# raw query with pivoted taxonomy, and join translation name but not roll it up.... then summ in here
  taxonomy = query("taxonomy hearchy", sampleId)
  # can be
  CSCIbugs = bugsRaw
  CSCIbugs$SampleID =
    CSCIbugs$StationCode = sites$siteId#### need to get site info!!!
  CSCIbugs$FinalID = CSCIbugs$OTUName
  CSCIbugs$LifeStageCode = ifelse(CSCIbugs$lifeStage)
  # #LifeStageCode = case
  # when class <> 'insecta' then 'X'
  # when class is null then 'X'
  # else Lifestage
  # end,
  CSCIbugs$Distinct == 0
  CSCIbugs$BAResult = sum(CSCIbugs$splitCount, CSCIbugs$bigRare)# have query add this as a column
}


CO_bug_export<-function(boxId, translationId){
  CObugs=NAMCr::query("sampleTaxaTranslation",
                      translationId=translationId,
                      boxId=boxId)
  sites =query("samples", 
               include= c("sampleId",'siteId','sampleDate',"waterbodyName"),
               boxId = boxId
  )
  
  CObugs$Project="NAMC report"
  CObugs$Station=sites$siteId
  CObugs$Name=sites$siteName
  CObugs$Location=sites$waterbodyName
  CObugs$CollDate=sites$sampleDate
  CObugs$Organism=CObugs$OTUName
  CObugs$Individuals=CObugs$splitCount
  CObugs$Stage=CObugs$lifestage
  CObugs$CommentsTaxa
  CObugs$RepNum
  CObugs$Grids
  CObugs$CommentsSample
  CObugs$CommentsRep
}


#PREDATOR NBR pseudo-O/E
#(The PREDitctive Assessment Tool for ORegon)
#based on FORTRAN Subsample.exe
#code by Andrew Caudillo, Lab Manager of NAMC
#modified from code ideas proposed by David Fowler
#Essentially a null O/E model (poor model; ref O/E SD is 0.29)

OR_NBR_bug <- function(sampleId, translationId, fixedCount) {
  bugsOTU = query(
    "sampleMetrics",
    translationId = def_models$translationId,
    fixedCount = def_models$fixedCount,
    sampleIds = def_model_results$sampleId
  )
  rarefiedOTUTaxa = subset(bugsOTU, metricName == "Rarefied Taxa")
  rarefiedOTUTaxa = NAMCr::json.expand(rarefiedOTUTaxa, "metricValue")
  sumrarefiedOTUTaxa = rarefiedOTUTaxa  %>%
    dplyr::group_by(sampleId, OTUName) %>%
    dplyr::summarize(sumSplitCount = sum(splitCount)) # why are multiple records exported here per OTU???
  bugnew=sumrarefiedOTUTaxa
return(bugnew)
  }
  
