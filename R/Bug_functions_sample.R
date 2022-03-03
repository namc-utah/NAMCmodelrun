#' OE bug data matrix export
#'
#' @param sampleIds
#' @param translationId
#' @param fixedCount
#'
#' @return bug data sample by taxa presence/absence matrix
#' @export
#'
#' @examples
OE_bug_matrix<-function(sampleIds,translationId,fixedCount){
  bugsOTU = NAMCr::query("sampleTaxaTranslationRarefied",
                        sampleIds = sampleIds,
                        translationId = translationId,
                        fixedCount = fixedCount
                        )
    sumrarefiedOTUTaxa = bugsOTU  %>%
    dplyr::group_by(sampleId, otuName) %>%
    dplyr::summarize(sumSplitCount = sum(splitCount)) # why are multiple records exported here per OTU???

  sumrarefiedOTUTaxa$presence = ifelse(sumrarefiedOTUTaxa$sumSplitCount >=1, 1, 0)
  bugnew = tidyr::pivot_wider(sumrarefiedOTUTaxa,id_cols = "sampleId", names_from = "otuName",values_from = "presence")
  bugnew[is.na(bugnew)]<-0
  bugnew=as.data.frame(bugnew)
  rownames(bugnew)<-bugnew$sampleId
  bugnew<-bugnew[,-1]
return(bugnew)
  }


#' MMI metrics
#'
#' @param sampleIds
#' @param translationId
#' @param fixedCount
#'
#' @return sample by metrics dataframe, only metrics needed for a given model
#' @export
#'
#' @examples
MMI_metrics<-function(sampleIds,translationId,fixedCount){
  # need to add list of required metrics to model table as well as random forest input file
  #get needed list of metrics from new end point
  # NV and AREMP relative abundances were OTU standardized too!
  bugsMetrics = NAMCr::query("sampleMetrics",
                             sampleIds=sampleIds,
                             translationId=translationId,
                             fixedCount=fixedCount
                             )
  #modelInfo = NAMCr::query("modelInfo", modelId = def_model$modelId)
  #MMI_metrics = subset(bugsMetrics, metricId %in% ()) # store translations to metric names in database
  # need to replace metric name with a model specific metric abbreviation
  if (def_models$translationId==13){

    # NV MMI metrics
    MMI_metrics = subset(bugsMetrics, metricId %in% c(279,#insect richness
                                                      283,# noninsect richness
                                                      295,#clinger richness
                                                      444, #shannons diversity
                                                      203, # collector filterer density, 438 relative abundance
                                                      178,#percent  Ephemeoptera
                                                      179))# percent plecoptera
    MMI_metrics$metricValue=as.numeric(MMI_metrics$metricValue)
    MMI_metrics$metricModelName=ifelse(MMI_metrics$metricId==71,"INSET",
                                       ifelse(MMI_metrics$metricId==73,"NONSET",
                                              ifelse(MMI_metrics$metricId==47,"CLINGER",
                                                     ifelse(MMI_metrics$metricId==22,"SHDIVER",
                                                            ifelse(MMI_metrics$metricId==42,"CFA",
                                                                          ifelse(MMI_metrics$metricId==362,"EPHEA",
                                                                                 ifelse(MMI_metrics$metricId==363,"PLECA",NA)
                                                                          ))))))

    bugnew = tidyr::pivot_wider(MMI_metrics,id_cols = "sampleId",names_from = "metricModelName",values_from = "metricValue")
    bugnew$PER_CFA=bugnew$CFA*100
    bugnew$PER_EPHEA=bugnew$EPHEA*100
    bugnew$PER_PLECA=bugnew$PLECA*100
    bugnew=bugnew[,c("sampleId","INSET","NONSET","CLINGER","SHDIVER","PER_CFA","PER_EPHEA","PER_PLECA")]
  }else if (def_models$translationId==23){

    #AREMP MMI metrics
    MMI_metrics = subset(bugsMetrics, metricId %in% c(295,#clinger richness
                                                      186,#EPT
                                                      271,#diptera richness
                                                      208,#intolerant
                                                      283,# noninsect richness
                                                      294#long lived taxa richness
    ))
    MMI_metrics$metricValue=as.numeric(MMI_metrics$metricValue)
    MMI_metrics$metricModelName=ifelse(MMI_metrics$metricId==47,"CLING_rich",
                                       ifelse(MMI_metrics$metricId==26,"EPT",
                                              ifelse(MMI_metrics$metricId==61,"DIPT_rich",
                                                     ifelse(MMI_metrics$metricId==33,"INTOL",
                                                                   ifelse(MMI_metrics$metricId==73,"NON_INSECT_rich",
                                                                          ifelse(MMI_metrics$metricId==48,"LLT_rich",NA)
                                                                   )))))
    bugnew = tidyr::pivot_wider(MMI_metrics,id_cols = "sampleId",names_from = "metricModelName",values_from = "metricValue")
    bugnew$PER_EPT=bugnew$EPT*100
    bugnew$PER_INTOL=bugnew$INTOL*100
    bugnew=bugnew[,c("sampleId","CLING_rich","PER_EPT","DIPT_rich","PER_INTOL","NON_INSECT_rich","LLT_rich")]
  }else{

  }


  bugnew=as.data.frame(bugnew)
  rownames(bugnew)<-bugnew$sampleId
  bugnew<-bugnew[,-1]
  return(bugnew)
}



#' CA CSCI bug data export
#'
#' @param sampleIds
#'
#' @return raw bug data translated by CSCI translation and with column names formatted for CSCI model function
#' @export
#'
#' @examples
CSCI_bug <- function(sampleIds){
# get needed data from the APIs
  bugRaw = NAMCr::query(
    "sampleTaxa",
    sampleIds = sampleIds
  )# raw NAMCr::query with pivoted taxonomy, and join translation name but not roll it up.... then summ in here

  bugsTranslation = NAMCr::query(
    "sampleTaxaTranslation",
    translationId = 9,
    sampleIds = sampleIds
  )
  sites = NAMCr::query(
    "samples",
    include = c("sampleId", 'siteName'),
    sampleIds = sampleIds
  )

 # join that data together into a single dataframe
  CSCIbugs=dplyr::left_join(bugRaw,bugsTranslation, by=c("taxonomyId", "SampleId"))
  CSCIbugs=dplyr::left_join(CSCIbugs,sites, by='SampleId')

  # rename columns
  # rename columns
  CSCIbugs$boxId = CSCIbugs$boxId
  CSCIbugs$StationCode = CSCIbugs$siteName#### need to get site info!!!
  CSCIbugs$FinalID = CSCIbugs$otuName
  CSCIbugs$class<-ifelse(is.na(CSCIbugs$class)==TRUE,"X",CSCIbugs$class)
  CSCIbugs$LifeStageCode=ifelse(is.na(CSCIbugs$lifeStageAbbreviation)==TRUE,"X",ifelse(CSCIbugs$class!="Insecta","X",CSCIbugs$lifeStageAbbreviation))
  CSCIbugs$Distinct = 0
  CSCIbugs=CSCIbugs %>% mutate(BAResult=splitCount.x+bigRareCount)
  CSCIbugs=subset(CSCIbugs,is.na(otuName)==FALSE)
  names(CSCIbugs)[names(CSCIbugs)=="sampleId"]<-"SampleID"
  # subset columns
  CSCIbugs=CSCIbugs[,c('SampleID','StationCode','FinalID','LifeStageCode','Distinct','BAResult')]
  bugnew= CSCIbugs %>% dplyr::group_by(SampleID,StationCode,FinalID,LifeStageCode,Distinct) %>% dplyr::summarize(BAResult=sum(BAResult))
  return(bugnew)
  }

#' CO EDAS 2017 bug data export by box
#'
#' @param sampleIds
#'
#' @return raw bug data translated by CO EDAS translation and with column names formatted for CO EDAS access database
#' @export
#'
#' @examples
CO_bug_export<-function(sampleIds){
  bugRaw = NAMCr::query(
    "sampleTaxa",
    sampleIds=sampleIds
  )# raw NAMCr::query with pivoted taxonomy, and join translation name but not roll it up.... then summ in here

  bugsTranslation = NAMCr::query(
    "sampleTaxaTranslation",
    translationId = 8,
    sampleIds=sampleIds
  )
  samples = NAMCr::query(
    "samples",
    include = c("boxId","sampleId",'siteId','sampleDate',"siteName", "sampleMethod","habitatName","area"),# possibly add waterbody name to the samples NAMCr::query
    sampleIds=sampleIds
  )

  # join that data together into a single dataframe
  CObugs=dplyr::left_join(bugRaw,bugsTranslation, by=c("taxonomyId", "sampleId"))
  CObugs=dplyr::left_join(CObugs,samples, by='sampleId')

  CObugs$Project=paste0(boxId)
  CObugs$Station=CObugs$sampleId
  CObugs$Name=CObugs$siteId
  CObugs$Location=CObugs$siteName # previously used waterbody name.. use that if we export this data for use by CO state
  CObugs$CollDate=CObugs$sampleDate
  CObugs$Organism=CObugs$otuName
  CObugs$Individuals=CObugs$splitCount.y
  CObugs$Stage=CObugs$lifeStageAbbreviation
  CObugs$CommentsTaxa=paste0("taxonomyId: ",CObugs$taxonomyId)
  CObugs$RepNum=1
  CObugs$Grids=NA
  CObugs$CommentsSample=paste0("sampleMethod: ",CObugs$sampleMethod, "area: ",CObugs$area)
  CObugs$CommentsRep=paste0("habitat: ",CObugs$habitatName)
  CObugs=CObugs[,c("Project","Station","Name","Location","CollDate","Organism","Individuals","Stage","CommentsTaxa","RepNum","Grids","CommentsSample","CommentsRep")]
  #write excel file to workspace
  write.csv(CObugs,file = paste0("CObugs","boxId_",CObugs$Project[1],"_",Sys.Date(),".csv"),row.names=FALSE)
  cat(paste("csv with CObugs has been written out to your current working directory.",
            "Convert this csv to excel 2003 and import into CO EDAS access database to compute the CSCI score.",
            "Follow instructions in this pdf Box\\NAMC\\OE_Modeling\\NAMC_Supported_OEmodels\\CO\\Documentation\\EDAS2017\\Tutorial Guide to EDAS_Version 1.7.pdf",
            "to import bug and habitat data, harmonize taxa list, rarefy and compute MMI",
            "then read resulting excel file back into R to save results in the database.", sep="\n"))
  return(CObugs)
}


#Essentially a null O/E model (poor model; ref O/E SD is 0.29)

#' OR NBR eastern
#'
#' @param sampleIds
#' @param translationId
#' @param fixedCount
#'
#' @return list of rarefied and translated taxa per sample
#' @export
#'
#' @examples
OR_NBR_bug <- function(sampleIds, translationId, fixedCount) {
  bugsOTU = NAMCr::query(
    "sampleTaxaTranslationRarefied",
    translationId = translationId,
    fixedCount = fixedCount,
    sampleIds = sampleIds
  )
  sumrarefiedOTUTaxa = bugsOTU  %>%
    dplyr::group_by(sampleId, otuName) %>%
    dplyr::summarize(sumSplitCount = sum(splitCount)) # why are multiple records exported here per OTU???
  bugnew=sumrarefiedOTUTaxa
return(bugnew)
  }

