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
  bugsOTU = NAMCr::query("sampleTaxaTranslationRarefied",
                        translationId = def_models$translationId,
                        fixedCount = def_models$fixedCount,
                        sampleIds = def_model_results$sampleId
                        )
    sumrarefiedOTUTaxa = bugsOTU  %>%
    dplyr::group_by(sampleId, otuName) %>%
    dplyr::summarize(sumSplitCount = sum(splitCount)) # why are multiple records exported here per OTU???

  sumrarefiedOTUTaxa$presence = ifelse(sumrarefiedOTUTaxa$sumSplitCount >=1, 1, 0)
  bugnew = tidyr::pivot_wider(sumrarefiedOTUTaxa,id_cols = "sampleId", names_from = "otuName",values_from = "presence")
  bugnew[is.na(bugnew)]<-0
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
MMI_metrics<-function(sampleId,translationId,fixedCount){
  # need to add list of required metrics to model table as well as random forest input file
  #get needed list of metrics from new end point
  # NV and AREMP relative abundances were OTU standardized too!
  bugsMetrics = NAMCr::query("sampleMetrics",
                             fixedCount=def_models$fixedCount,
                             translationId=def_models$translationId,
                             sampleIds = def_model_results$SampleId)
  #modelInfo = NAMCr::query("modelInfo", modelId = modelId)
  #MMI_metrics = subset(bugsMetrics, metricId %in% ()) # store translations to metric names in database
  # need to replace metric name with a model specific metric abbreviation
  bugnew = tidyr::pivot_wider(bugsMetrics,id_cols = "sampleId",names_from = "metricName",values_from = "metricValue")
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
CSCI_bug <- function(sampleId){
# get needed data from the APIs
  bugRaw = NAMCr::query(
    "sampleTaxa",
    sampleIds = def_model_results$SampleId
  )# raw NAMCr::query with pivoted taxonomy, and join translation name but not roll it up.... then summ in here

  bugsTranslation = NAMCr::query(
    "sampleTaxaTranslation",
    translationId = 9,
    sampleIds = def_model_results$SampleId
  )
  sites = NAMCr::query(
    "samples",
    include = c("sampleId", 'siteId'),
    sampleIds = def_model_results$SampleId
  )

 # join that data together into a single dataframe
  CSCIbugs=dplyr::left_join(bugRaw,bugsTranslation, by=c("taxonomyId", "SampleId"))
  CSCIbugs=dplyr::left_join(CSCIbugs,sites, by='SampleId')

  # rename columns
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
  bugsOTU = NAMCr::query(
    "sampleTaxaTranslationRarefied",
    translationId = def_models$translationId,
    fixedCount = def_models$fixedCount,
    sampleIds = def_model_results$sampleId
  )
  sumrarefiedOTUTaxa = bugsOTU  %>%
    dplyr::group_by(sampleId, otuName) %>%
    dplyr::summarize(sumSplitCount = sum(splitCount)) # why are multiple records exported here per OTU???
  bugnew=sumrarefiedOTUTaxa
return(bugnew)
  }

