#' OE bug data matrix export by box
#'
#' @param boxId
#' @param translationId
#' @param fixedCount
#'
#' @return bug data sample by taxa presence/absence matrix
#' @export
#'
#' @examples
OE_bug_matrix_box<-function(boxId,translationId,fixedCount){
  bugsOTU = NAMCr::query("sampleTaxaTranslationRarefied",
                         translationId = def_models$translationId,
                         fixedCount = def_models$fixedCount,
                         boxId = boxId
  )
  sumrarefiedOTUTaxa = bugsOTU  %>%
    dplyr::group_by(sampleId, otuName) %>%
    dplyr::summarize(sumSplitCount = sum(splitCount)) # why are multiple records exported here per OTU???

  sumrarefiedOTUTaxa$presence = ifelse(sumrarefiedOTUTaxa$sumSplitCount >=1, 1, 0)
  bugnew = tidyr::pivot_wider(sumrarefiedOTUTaxa,id_cols = "sampleId", names_from = "otuName",values_from = "presence")
  bugnew[is.na(bugnew)]<-0
  return(bugnew)
}


#' MMI metrics by box
#'
#' @param boxId
#' @param fixedCount
#'
#' @return sample by metrics dataframe, only metrics needed for a given model
#' @export
#'
#' @examples
MMI_metrics<-function(boxId,translationId,fixedCount){
  # need to add list of required metrics to model table as well as random forest input file
  #get needed list of metrics from new end point
  # NV and AREMP relative abundances were OTU standardized too!
  bugsMetrics = NAMCr::query("sampleMetrics",
                             fixedCount=def_models$fixedCount,
                             translationId=def_models$translationId,
                             boxId = boxId)
  modelInfo = NAMCr::query("modelInfo", modelId = def_model_results$modelId)
  #MMI_metrics = subset(bugsMetrics, metricId %in% ()) # store translations to metric names in database
  # need to replace metric name with a model specific metric abbreviation
  bugnew = tidyr::pivot_wider(MMI_metrics,id_cols = "sampleId",names_from = "metricName",values_from = "metricValue")
  return(bugnew)
}




#' CA CSCI bug data export by box
#'
#' @param boxId
#' @param translationId
#'
#' @return raw bug data translated by CSCI translation and with column names formatted for CSCI model function
#' @export
#'
#' @examples
CSCI_bug_box <- function(boxId){
  # get needed data from the APIs
  bugRaw = NAMCr::query(
    "sampleTaxa",
    boxId = boxId
  )# raw NAMCr::query with pivoted taxonomy, and join translation name but not roll it up.... then summ in here

  bugsTranslation = NAMCr::query(
    "sampleTaxaTranslation",
    translationId = 9,
    boxId = boxId
  )
  sites = NAMCr::query(
    "samples",
    include = c("sampleId", 'siteId'),
    boxId = boxId
  )

  # join that data together into a single dataframe
  CSCIbugs=dplyr::left_join(bugRaw,bugsTranslation, by=c("taxonomyId", "sampleId"))
  CSCIbugs=dplyr::left_join(CSCIbugs,sites, by='sampleId')

  # rename columns
  CSCIbugs$boxId = CSCIbugs$boxId
  CSCIbugs$StationCode = CSCIbugs$siteId#### need to get site info!!!
  CSCIbugs$FinalID = CSCIbugs$otuName
  CSCIbugs$class<-ifelse(is.na(CSCIbugs$class)==TRUE,"X",CSCIbugs$class)
  CSCIbugs$LifeStageCode=ifelse(is.na(CSCIbugs$lifeStageAbbreviation)==TRUE,"X",ifelse(CSCIbugs$class=="Insecta","X",CSCIbugs$lifeStageAbbreviation))
  CSCIbugs$Distinct = 0
  CSCIbugs=CSCIbugs %>% mutate(BAResult=splitCount.y+bigRareCount)
  CSCIbugs=subset(CSCIbugs,is.na(otuName)==FALSE)
  # subset columns
  CSCIbugs=CSCIbugs[,c('sampleId','StationCode','FinalID','LifeStageCode','Distinct','BAResult')]
  return(CSCIbugs)
}


#' CO EDAS 2017 bug data export by box
#'
#' @param boxId
#'
#' @return raw bug data translated by CO EDAS translation and with column names formatted for CO EDAS access database
#' @export
#'
#' @examples
CO_bug_export_box<-function(boxId){
  bugRaw = NAMCr::query(
    "sampleTaxa",
    translationId = 8,
    boxId = boxId
  )# raw NAMCr::query with pivoted taxonomy, and join translation name but not roll it up.... then summ in here

  bugsTranslation = NAMCr::query(
    "sampleTaxaTranslation",
    translationId = 8,
    boxId=BoxId
  )
  samples = NAMCr::query(
    "samples",
    include = c("boxId",'siteId','sampleDate',"siteName", "sampleMethod","habitat","area"),# possibly add waterbody name to the samples NAMCr::query
    boxId=boxId
  )

  # join that data together into a single dataframe
  CObugs=dplyr::left_join(bugRaw,bugsTranslation, by=c("taxonomyId", "sampleId"))
  CObugs=dplyr::left_join(CObugs,sites, by='sampleId')

  CObugs$Project="NAMC report"
  CObugs$Station=CObugs$boxId
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
  #write excel file to workspace
  return(CObugs)
}


#Essentially a null O/E model (poor model; ref O/E SD is 0.29)

#' OR NBR eastern by box
#'
#' @param boxId
#' @param translationId
#' @param fixedCount
#'
#' @return list of rarefied and translated taxa per sample
#' @export
#'
#' @examples
OR_NBR_bug_box <- function(boxId, translationId, fixedCount) {
  bugsOTU = NAMCr::query(
    "sampleTaxaTranslationRarefied",
    translationId = def_models$translationId,
    fixedCount = def_models$fixedCount,
    boxId = boxId
  )
  sumrarefiedOTUTaxa = bugsOTU  %>%
    dplyr::group_by(sampleId, otuName) %>%
    dplyr::summarize(sumSplitCount = sum(splitCount)) # why are multiple records exported here per OTU???
  bugnew=sumrarefiedOTUTaxa
  return(bugnew)
}

