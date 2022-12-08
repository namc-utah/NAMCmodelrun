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
MMI_metrics<-function(sampleIds,translationId,fixedCount,modelId){
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
  if (def_models$modelId==3){

    # NV MMI metrics
    MMI_metrics = subset(bugsMetrics, metricId %in% c(276,#insect richness
                                                      280,# noninsect richness
                                                      291,#clinger richness
                                                      439, #shannons diversity
                                                      113, # collector filter relative abundance
                                                      89,#percent  Ephemeoptera
                                                      90))# percent plecoptera
    MMI_metrics$metricValue=as.numeric(MMI_metrics$metricValue)
    MMI_metrics$metricModelName=ifelse(MMI_metrics$metricId==276,"INSET",
                                       ifelse(MMI_metrics$metricId==280,"NONSET",
                                              ifelse(MMI_metrics$metricId==291,"CLINGER",
                                                     ifelse(MMI_metrics$metricId==439,"SHDIVER",
                                                            ifelse(MMI_metrics$metricId==113,"CFA",
                                                                          ifelse(MMI_metrics$metricId==89,"EPHEA",
                                                                                 ifelse(MMI_metrics$metricId==90,"PLECA",NA)
                                                                          ))))))

    bugnew = tidyr::pivot_wider(MMI_metrics,id_cols = "sampleId",names_from = "metricModelName",values_from = "metricValue")
    bugnew$PER_CFA=bugnew$CFA*100
    bugnew$PER_EPHEA=bugnew$EPHEA*100
    bugnew$PER_PLECA=bugnew$PLECA*100
    bugnew=bugnew[,c("sampleId","INSET","NONSET","CLINGER","SHDIVER","PER_CFA","PER_EPHEA","PER_PLECA")]
  }else if (def_models$modelId==8){

    #AREMP MMI metrics
    MMI_metrics = subset(bugsMetrics, metricId %in% c(291,#clinger richness
                                                      97,#EPT
                                                      268,#diptera richness
                                                      118,#intolerant
                                                      280,# noninsect richness
                                                      290#long lived taxa richness
    ))
    MMI_metrics$metricValue=as.numeric(MMI_metrics$metricValue)
    MMI_metrics$metricModelName=ifelse(MMI_metrics$metricId==291,"CLING_rich",
                                       ifelse(MMI_metrics$metricId==97,"EPT",
                                              ifelse(MMI_metrics$metricId==268,"DIPT_rich",
                                                     ifelse(MMI_metrics$metricId==118,"INTOL",
                                                                   ifelse(MMI_metrics$metricId==280,"NON_INSECT_rich",
                                                                          ifelse(MMI_metrics$metricId==290,"LLT_rich",NA)
                                                                   )))))
    bugnew = tidyr::pivot_wider(MMI_metrics,id_cols = "sampleId",names_from = "metricModelName",values_from = "metricValue")
    bugnew$PER_EPT=bugnew$EPT*100
    bugnew$PER_INTOL=bugnew$INTOL*100
    bugnew=bugnew[,c("sampleId","CLING_rich","PER_EPT","DIPT_rich","PER_INTOL","NON_INSECT_rich","LLT_rich")]

  }else if (def_models$modelId==136){

  # arid west modeled insect richness
    MMI_metrics = subset(bugsMetrics, metricId %in% c(364,#unique insect richness
                                                      362)) #unique richness midges
    MMI_metrics$metricValue=as.numeric(MMI_metrics$metricValue)
    MMI_metrics$metricAbbreviation=gsub("-","_",MMI_metrics$metricAbbreviation)
    bugnew = tidyr::pivot_wider(MMI_metrics,id_cols = "sampleId",names_from = "metricAbbreviation",values_from = "metricValue")
    bugnew$UniqueRichness_Insecta=bugnew$UniqueRichness_Insecta-bugnew$UniqueRichness_Chironomidae
  }else if (def_models$modelId==169){
    MMI_metrics = subset(bugsMetrics, metricId %in% c(351,#total richness
                                                      352,#ephem richness
                                                      354,#tricoptera richness
                                                      356,# diptera richness
                                                      375,#scraper richness
                                                      381,#intolerant taxa
                                                      444,#hilsenhoff
                                                      176,#ephem relative abundance
                                                      177,#plecoptera relative abundance
                                                      199,#scraper relative abundance
                                                      448#dominant taxa abundance
                                                      ))
    MMI_metrics$metricValue=as.numeric(MMI_metrics$metricValue)
    MMI_metrics$metricAbbreviation=gsub("-","_",MMI_metrics$metricAbbreviation)
    MMI_metrics$metricValue=ifelse(MMI_metrics$metricId %in% c(176,177,199,448),MMI_metrics$metricValue*100,MMI_metrics$metricValue)# converting relative abundance to %
      bugnew=MMI_metrics
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
  CSCIbugs=dplyr::left_join(bugRaw,bugsTranslation, by=c("taxonomyId", "sampleId"))
  CSCIbugs=dplyr::left_join(CSCIbugs,sites, by='sampleId')

  # rename columns
  # rename columns
  CSCIbugs$boxId = CSCIbugs$boxId
  CSCIbugs$StationCode = CSCIbugs$sampleId#### need to get site info!!!
  CSCIbugs$FinalID = CSCIbugs$otuName
  CSCIbugs$class<-ifelse(is.na(CSCIbugs$class)==TRUE,"X",CSCIbugs$class)
  CSCIbugs$LifeStageCode=ifelse(CSCIbugs$lifeStageAbbreviation=='U',"L",CSCIbugs$lifeStageAbbreviation)
  CSCIbugs$LifeStageCode=ifelse(is.na(CSCIbugs$LifeStageCode)==TRUE,"X",ifelse(CSCIbugs$class!="Insecta","X",CSCIbugs$LifeStageCode))
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

  CObugs$Project=CObugs$boxId
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


#' AZ EDAS bug data export
#'
#' @param sampleIds
#'
#' @return raw bug data translated by AZ EDAS translation and with column names formatted for AZ EDAS access database
#' @export
#'
#' @examples
AZ_bug_export<-function(sampleIds){
  bugRaw = NAMCr::query(
    "sampleTaxaUnambiguous",
    sampleIds=sampleIds
  )# unique unrarefied taxa NAMCr::query with pivoted taxonomy, and join translation name but not roll it up.... then summ in here

  bugRare =
  bugsTranslation = NAMCr::query(
  "sampleTaxaTranslation",
  translationId = 26,
  sampleIds=sampleIds
  )
  samples = NAMCr::query(
    "samples",
    include = c("boxId","sampleId",'siteId','sampleDate',"siteName", "waterbodyName","sampleMethod","habitatName","area"),# possibly add waterbody name to the samples NAMCr::query
    sampleIds=sampleIds
  )

  # join that data together into a single dataframe
  AZbugs=dplyr::left_join(bugRaw,bugsTranslation, by=c("taxonomyId", "sampleId"))
  AZbugs=dplyr::left_join(bugRaw,samples, by='sampleId')

  AZbugs$StationID=AZbugs$siteName
  AZbugs$WaterbodyName=AZbugs$waterbodyName
  AZbugs$ActivityID=AZbugs$sampleId
  AZbugs$RepNum=1
  AZbugs$CommentsSample=NA
  AZbugs$CollDate=AZbugs$sampleDate
  AZbugs$CollMeth=AZbugs$sampleMethod
  AZbugs$CorrectionFactor=AZbugs$labSplit
  AZbugs$FinalID=AZbugs$scientificName
  AZbugs$Individuals=AZbugs$splitCount
  AZbugs$Stage="L"
  AZbugs$LargeRare='No'
  AZbugs$Habitat=ifelse(AZbugs$habitatName=='Targeted Riffle','Riffle','Multi-Habitat')
  AZbugs$Lab='NAMC'
  AZbugs$LabID=NA

  AZbugs2=AZbugs[,c("StationID","WaterbodyName","ActivityID","RepNum","CollDate","CommentsSample","CollMeth","CorrectionFactor","FinalID","Individuals","Stage","LargeRare","Habitat","Lab","LabID")]
  #write excel file to workspace
  write.csv(AZbugs2,file = paste0("AZbugs","boxId_",boxId,"_",Sys.Date(),".csv"),row.names=FALSE)
  cat(paste("csv with AZbugs has been written out to your current working directory.",
            "convert this csv to excel 2003 and import into AZ EDAS access database to compute the IBI score.",
            "/OE_Modeling/NAMC_Supported_OEmodels/AZ/Benthic Data Bulk Upload_SOP.doc",
             "to import bug and habitat data, harmonize taxa list, rarefy and compute MMI",
            "then read resulting excel file back into R to save results in the database.", sep="\n"))
  return(AZbugs)
}

