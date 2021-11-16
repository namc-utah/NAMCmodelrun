#' OE bug data matrix export
#'
#' @param sampleId 
#' @param translationId 
#' @param fixedCount 
#'
#' @return bug data sample by taxa presence/absence matrix
#' @export
#'
#' @examples
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
return(bugnew)
  }


#' MMI metrics
#'
#' @param sampleId 
#' @param fixedCount 
#'
#' @return sample by metrics dataframe, only metrics needed for a given model
#' @export 
#'
#' @examples
MMI_metrics<-function(sampleId,fixedCount, modelId){
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
return(bugnew)
  }




#' CA CSCI bug data export
#'
#' @param sampleId 
#' @param translationId 
#'
#' @return raw bug data translated by CSCI translation and with column names formatted for CSCI model function
#' @export
#'
#' @examples
CSCI_bug <- function(sampleId, translationId){
# get needed data from the APIs  
  bugRaw = query(
    "sampleTaxaInfo",
    translationId = def_models$translationId,
    sampleIds = def_model_results$SampleId
  )# raw query with pivoted taxonomy, and join translation name but not roll it up.... then summ in here
  
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
  
 # join that data together into a single dataframe   
  CSCIbugs=dplyr::left_join(bugRaw,bugsTranslation, by=c("taxonomyId", "SampleId"))
  CSCIbugs=dplyr::left_join(CSCIbugs,sites, by='SampleId')
  
  # rename columns
  CSCIbugs$SampleID = CSCIbugs$sampleId
  CSCIbugs$StationCode = CSCIbugs$siteId#### need to get site info!!!
  CSCIbugs$FinalID = CSCIbugs$OTUName
  CSCIbugs$LifeStageCode=ifelse(CSCIbugs$Class=="Insecta","X",ifelse(is.na(CSCIbugs$lifeStageAbbreviation))==TRUE,X,CSCIbugs$lifeStageAbbreviation)
  CSCIbugs$Distinct == 0
  CSCIbugs$BAResult = sum(CSCIbugs$rawCount, CSCIbugs$rawBigRareCount)# have query add this as a column... use the mutate function?

  # subset columns
  CSCIbugs=CSCIbugs[,c('SampleID','StationCode','FinalID','LifeStageCode','Distinct','BAResult')]
  return(CSCIbugs)
  }


#' CO EDAS 2017 bug data export
#'
#' @param boxId 
#' @param translationId 
#'
#' @return raw bug data translated by CO EDAS translation and with column names formatted for CO EDAS access database
#' @export
#'
#' @examples
CO_bug_export<-function(boxId, translationId){
  bugRaw = query(
    "sampleTaxaInfo",
    translationId = def_models$translationId,
    sampleIds = def_model_results$SampleId
  )# raw query with pivoted taxonomy, and join translation name but not roll it up.... then summ in here
  
  bugsTranslation = query(
    "sampleTaxaTranslation",
    translationId = def_models$translationId,
    sampleIds = def_model_results$SampleId
  )
  samples = query(
    "samples",
    include = c("sampleId",'siteId','sampleDate',"siteName", "sampleMethod","habitat","area"),# possibly add waterbody name to the samples query
    sampleIds = def_model_results$SampleId
  )
  
  # join that data together into a single dataframe   
  CObugs=dplyr::left_join(bugRaw,bugsTranslation, by=c("taxonomyId", "SampleId"))
  CObugs=dplyr::left_join(CObugs,sites, by='SampleId')
  
  CObugs$Project="NAMC report"
  CObugs$Station=CObugs$sampleId
  CObugs$Name=CObugs$siteId
  CObugs$Location=CObugs$siteName # previously used waterbody name.. use that if we export this data for use by CO state
  CObugs$CollDate=CObugs$sampleDate
  CObugs$Organism=CObugs$OTUName
  CObugs$Individuals=CObugs$rawCount
  CObugs$Stage=CObugs$lifestageAbbreviation
  CObugs$CommentsTaxa
  CObugs$RepNum==1
  CObugs$Grids=NA
  CObugs$CommentsSample==paste0("sampleMethod: ",CObugs$sampleMethod,"habitat: ",CObugs$habitat,"area: ",CObugs$area)
  CObugs$CommentsRep==""
return(CObugs)
  }


#Essentially a null O/E model (poor model; ref O/E SD is 0.29)

#' OR NBR eastern
#'
#' @param sampleId 
#' @param translationId 
#' @param fixedCount 
#'
#' @return list of rarefied and translated taxa per sample 
#' @export
#'
#' @examples
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
  
