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
      #NRSA xeric MMI
  }else if (def_models$modelId==367){
   taxa=query("sampleTaxa",sampleIds=sampleIds)
   count=taxa %>%
     dplyr::group_by(sampleId) %>%
   dplyr::summarize(count = sum(splitCount))

   taxatop=taxa %>%
     group_by(sampleId) %>%
     top_n(5,splitCount)
   taxatopcount=taxatop %>%
     dplyr::summarize(topcount = sum(splitCount))


   taxatopcountfinal=left_join(count,taxatopcount,by="sampleId")
   taxatopcountfinal$pcttop5taxa=taxatopcountfinal$topcount/taxatopcountfinal$count

     MMI_metrics = subset(bugsMetrics, metricId %in% c(448,351,448,360,380,384,390,628))
                         MMI_metrics$metricValue=as.numeric(MMI_metrics$metricValue)
                         MMI_metrics$metricAbbreviation=gsub("-","_",MMI_metrics$metricAbbreviation)
                         bugnew = tidyr::pivot_wider(MMI_metrics,id_cols = "sampleId",names_from = "metricAbbreviation",values_from = "metricValue")
                         bugnew=left_join(bugnew,taxatopcountfinal, by="sampleId")
                         bugnew$pctTolerantTaxa=bugnew$UniqueRichness_TolerantTaxa/bugnew$UniqueRichness*100
                         bugnew$pctclinger=bugnew$UniqueRichness_Habit_prim_Clinger/bugnew$UniqueRichness*100
                         bugnew$noninsect_stand=ifelse(1-(bugnew$RelativeAbundance_NonInsects-0.0333)/(0.36-0.033)>1,1,
                                                       ifelse(1-(bugnew$RelativeAbundance_NonInsects-0.0333)/(0.36-0.033)<0,0,
                                                              1-(bugnew$RelativeAbundance_NonInsects-0.0333)/(0.36-0.033)))
                         bugnew$top5_stand=ifelse(1-(bugnew$pcttop5taxa-0.406)/(0.823-0.406)>1,1,
                                                  ifelse(1-(bugnew$pcttop5taxa-0.406)/(0.823-0.406)<0,0,
                                                         1-(bugnew$pcttop5taxa-0.406)/(0.823-0.406)))
                         bugnew$scraper_stand=ifelse(bugnew$UniqueRichness_Feed_prim_abbrev_SC/7>1,1,bugnew$UniqueRichness_Feed_prim_abbrev_SC/7)

                         bugnew$clinger_stand=ifelse((bugnew$pctclinger-15.8)/(65.8-15.8)>1,1,
                                                     ifelse((bugnew$pctclinger-15.8)/(65.8-15.8)<0,0,
                                                            (bugnew$pctclinger-15.8)/(65.8-15.8)))
                         bugnew$EPT_stand=ifelse((bugnew$pctclinger-1)/(18-1)>1,1,
                                                 ifelse((bugnew$pctclinger-1)/(18-1)<0,0,
                                                        (bugnew$pctclinger-1)/(18-1)))
                         bugnew$tolerant_stand=ifelse(1-(bugnew$pctTolerantTaxa-3.57)/(36.4-3.57)>1,1,
                                                      ifelse(1-(bugnew$pctTolerantTaxa-3.57)/(36.4-3.57)<0,0,
                                                             1-(bugnew$pctTolerantTaxa-3.57)/(36.4-3.57)))

                           t=mutate(bugnew,MMI = rowMeans(select(bugnew,c(noninsect_stand:tolerant_stand))))
                           t$MMI=t$MMI*100
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
  bugRaw =
    NAMCr::query(
    "sampleTaxaTranslationRarefied",
    sampleIds=sampleIds,
    translationId=8,
    fixedCount=300
  )# raw NAMCr::query with pivoted taxonomy, and join translation name but not roll it up.... then summ in here

  bugLives = NAMCr::query(
    "sampleTaxa",
   sampleIds=sampleIds
  )
  samples = NAMCr::query(
    "samples",
    include = c("boxId","sampleId",'siteId','sampleDate',"siteName", "sampleMethod","habitatName","area"),# possibly add waterbody name to the samples NAMCr::query
    sampleIds=sampleIds
  )
bugRaw<-plyr::join(bugRaw,bugLives[,c('taxonomyId','lifeStageAbbreviation')],by='taxonomyId')
bugRaw$burner<-paste(bugRaw$sampleId,bugRaw$otuName,bugRaw$splitCount)
bugRaw<-bugRaw[!duplicated(bugRaw$burner),]
  # join that data together into a single dataframe
  #CObugs=dplyr::left_join(bugRaw,bugsTranslation, by=c("taxonomyId", "sampleId"))
  CObugs=dplyr::left_join(bugRaw,samples, by='sampleId')
  CObugs2<-as.data.frame(matrix(nrow=nrow(CObugs),ncol=13))
  names(CObugs2)<-c("Project","Station","Name","Location","CollDate","Organism","Individuals","Stage","CommentsTaxa","RepNum","Grids","CommentsSample","CommentsRep")#13 is number for bug import
  CObugs2$Project=CObugs$boxId
  CObugs2$Station=CObugs$sampleId
  CObugs2$Name=CObugs$siteId
  CObugs2$Location=CObugs$siteName # previously used waterbody name.. use that if we export this data for use by CO state
  CObugs2$CollDate=CObugs$sampleDate
  CObugs2$Organism=CObugs$otuName
  CObugs2$Individuals=CObugs$splitCount
  CObugs2$Stage=CObugs$lifeStageAbbreviation
  CObugs2$CommentsTaxa=paste0("taxonomyId: ",CObugs$taxonomyId)
  CObugs2$RepNum=1
  CObugs2$Grids=NA
  CObugs2$CommentsSample=paste0("sampleMethod: ",CObugs$sampleMethod, "area: ",CObugs$area)
  CObugs2$CommentsRep=paste0("habitat: ",CObugs$habitatName)
  #CObugs=CObugs[,c("Project","Station","Name","Location","CollDate","Organism","Individuals","Stage","CommentsTaxa","RepNum","Grids","CommentsSample","CommentsRep")]
  #write excel file to workspace
  write.csv(CObugs2,file = paste0(CO_path,"CObugs","boxId_",CObugs2$Project[1],"_",Sys.Date(),".csv"),row.names=FALSE)
  cat(paste("csv with CO bugs has been written out to the CO EDAS imports folder.",
            "Convert this csv to excel 2003 and import into CO EDAS access database to compute the CSCI score.",
            "Follow instructions in this pdf Box\\NAMC\\OE_Modeling\\NAMC_Supported_OEmodels\\CO\\Documentation\\EDAS2017\\Tutorial Guide to EDAS_Version 1.7.pdf",
            "to import bug and habitat data, harmonize taxa list, rarefy and compute MMI",
            "then read resulting excel file back into R to save results in the database.", sep="\n"))
  return(CObugs)
}
CO_v5_bug_export<-function(sampleIds){
  bugRaw =
    NAMCr::query(
      "sampleTaxaTranslationRarefied",
      sampleIds=sampleIds,
      translationId=92,
      fixedCount=300
    )# raw NAMCr::query with pivoted taxonomy, and join translation name but not roll it up.... then summ in here

  bugLives = NAMCr::query(
    "sampleTaxa",
    sampleIds=sampleIds
  )
  samples = NAMCr::query(
    "samples",
    include = c("boxId","sampleId",'siteId','sampleDate',"siteName", "sampleMethod","habitatName","area"),# possibly add waterbody name to the samples NAMCr::query
    sampleIds=sampleIds
  )
  bugRaw<-plyr::join(bugRaw,bugLives[,c('taxonomyId','lifeStageAbbreviation')],by='taxonomyId')
  bugRaw$burner<-paste(bugRaw$sampleId,bugRaw$otuName,bugRaw$splitCount)
  bugRaw<-bugRaw[!duplicated(bugRaw$burner),]
  # join that data together into a single dataframe
  #CObugs=dplyr::left_join(bugRaw,bugsTranslation, by=c("taxonomyId", "sampleId"))
  CObugs=dplyr::left_join(bugRaw,samples, by='sampleId')
  CObugs2<-as.data.frame(matrix(nrow=nrow(CObugs),ncol=13))
  #names(CObugs2)<-c("Project","Station","Name","Location","CollDate","Organism","Individuals","Stage","CommentsTaxa","RepNum","Grids","CommentsSample","CommentsRep")#13 is number for bug import
  CObugs2$Project=CObugs$boxId
  CObugs2$Station=CObugs$sampleId
  CObugs2$Name=CObugs$siteId
  CObugs2$Location=CObugs$siteName # previously used waterbody name.. use that if we export this data for use by CO state
  CObugs2$CollDate=CObugs$sampleDate
  CObugs2$Organism=CObugs$otuName
  CObugs2$Individuals=CObugs$splitCount
  CObugs2$Stage=CObugs$lifeStageAbbreviation
  CObugs2$CommentsTaxa=paste0("taxonomyId: ",CObugs$taxonomyId)
  CObugs2$RepNum=1
  CObugs2$Grids=NA
  CObugs2$CommentsSample=paste0("sampleMethod: ",CObugs$sampleMethod, "area: ",CObugs$area)
  CObugs2$CommentsRep=paste0("habitat: ",CObugs$habitatName)
  #CObugs=CObugs[,c("Project","Station","Name","Location","CollDate","Organism","Individuals","Stage","CommentsTaxa","RepNum","Grids","CommentsSample","CommentsRep")]
  #write excel file to workspace
  write.csv(CObugs2,file = paste0(CO_path,"CObugs","boxId_",CObugs2$Project[1],"_",Sys.Date(),".csv"),row.names=FALSE)
  cat(paste("csv with CO bugs has been written out to the CO EDAS imports folder.",
            "Convert this csv to excel 2003 and import into CO EDAS access database to compute the CSCI score.",
            "Follow instructions in this pdf Box\\NAMC\\OE_Modeling\\NAMC_Supported_OEmodels\\CO\\Documentation\\EDAS2017\\Tutorial Guide to EDAS_Version 1.7.pdf",
            "to import bug and habitat data, harmonize taxa list, rarefy and compute MMI",
            "then read resulting excel file back into R to save results in the database.", sep="\n"))
  return(CObugs)
}
CO_pred_export<-function(prednew){
  #make empty df
  prednews<-as.data.frame(matrix(ncol=15,nrow=nrow(prednew)))
  names(prednews)<-c('StationID', 'WaterbodyName','Location',
                     'Lat_Dec', 'Long_Dec','NHDSLOPE',
                     'ELEV_SITE','ECO3','ECO4',
                     'SUMMER','WINTER','LOG_XP_PT',
                     'SQRT_TOPO','PRCPSHORTWS','DOY')
  #fill the df with predictors that are named how EDAS wants
  prednews$StationID=row.names(prednew)
  prednews$WaterbodyName=def_samples$siteName
  prednews$Location=rep('Location',nrow(prednew))
  prednews$Lat_Dec=prednew$Lat_Dec
  prednews$Long_Dec=prednew$Long_Dec
  prednews$NHDSLOPE=prednew$NHDSLOPE
  prednews$ELEV_SITE=prednew$ELEV_SITE
  prednews$ECO3=prednew$ECO3
  prednews$ECO4=prednew$ECO4
  prednews$SUMMER=prednew$SUMMER
  prednews$WINTER=rep('',nrow(prednews))
  prednews$LOG_XP_PT=rep('',nrow(prednews))
  prednews$SQRT_TOPO=rep('',nrow(prednews))
  prednews$PRCPSHORTWS=rep('',nrow(prednews))
  prednews$DOY=prednew$DOY
  write.csv(prednews,file = paste0(CO_path,"COpreds","boxId_",CObugs$Project[1],"_",Sys.Date(),".csv"),row.names=FALSE)
  cat(paste("csv with CO predictors has been written out to the CO EDAS imports folder.",
            "Convert this csv to excel 2003 and import into CO EDAS access database to compute the CSCI score.",
            "Follow instructions in this pdf Box\\NAMC\\OE_Modeling\\NAMC_Supported_OEmodels\\CO\\Documentation\\EDAS2017\\Tutorial Guide to EDAS_Version 1.7.pdf",
            "to import bug and habitat data, harmonize taxa list, rarefy and compute MMI",
            "then read resulting excel file back into R to save results in the database.", sep="\n"))
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
    sampleIds=210531#sampleIds
  )# unique unrarefied taxa NAMCr::query with pivoted taxonomy, and join translation name but not roll it up.... then summ in here

  AZsubsamp<-rarify(inbug=bugRaw, sample.ID="sampleId", abund="splitCount", subsiz=500)
  AZsubsamp<-AZsubsamp[AZsubsamp$splitCount>0,]

  bugsTranslation = NAMCr::query(
    "sampleTaxaTranslation",
    translationId = 26,
    sampleIds=sampleIds
  )
  samples = NAMCr::query(
    "samples",
    include = c("boxId","sampleId",'siteId','sampleDate',"siteName", "waterbodyName","sampleMethod","habitatName","area", "labSplit"),
    sampleIds=sampleIds
  )

  # join that data together into a single dataframe
  AZbugs=dplyr::left_join(AZsubsamp,bugsTranslation[,c('otuName','taxonomyId','sampleId')], by=c("taxonomyId", "sampleId"))
  AZbugs=dplyr::left_join(AZbugs,samples, by='sampleId')

  sampletax = NAMCr::query(
    "sampleTaxa",
    sampleIds=sampleIds
  )
  sampletax<-sampletax[!duplicated(sampletax$scientificName),]
  AZbugs<-merge(AZbugs,sampletax[,c('scientificName','lifeStageAbbreviation')],'scientificName')
  AZbugs2<-data.frame(matrix(nrow=nrow(AZbugs)));names(AZbugs2)<-'StationID'
  AZbugs2$StationID=AZbugs$siteName
  AZbugs2$WaterbodyName=AZbugs$waterbodyName
  AZbugs2$ActivityID=AZbugs$sampleId
  AZbugs2$RepNum=1
  AZbugs2$CommentsSample=NA
  AZbugs2$CollDate=AZbugs$sampleDate
  AZbugs2$CollMeth=rep('ADEQ Riffle bugs',nrow(AZbugs))
  #correction factor is 100 * 1/labsplit %. So we need to force decimal to %.
  AZbugs2$CorrectionFactor=100*(1/(AZbugs$labSplit*100))
  AZbugs2$FinalID=AZbugs$otuName
  AZbugs2$Individuals=AZbugs$splitCount
  AZbugs2$Stage=AZbugs$lifeStageAbbreviation
  AZbugs2$LargeRare='No'
  AZbugs2$Habitat=ifelse(AZbugs$habitatName=='Targeted Riffle','Riffle','Multi-Habitat')
  AZbugs2$Lab='NAMC'
  AZbugs2$LabID=NA



  AZbugs2<-AZbugs2[which(is.na(AZbugs2$FinalID)==F),]

  #write excel file to workspace

  write.csv(AZbugs2,file = paste0("C://Users//andrew.caudillo//Box//NAMC//OEModeling//NAMC_Supported_OEmodels//Arizona//InputFiles//AZbugs","boxId_",boxId,"_",Sys.Date(),".csv"),row.names=FALSE)
  cat(paste("csv with AZbugs has been written out to your current working directory.",
            "convert this csv to excel 2003 and import into AZ EDAS access database to compute the IBI score.",
            "/OE_Modeling/NAMC_Supported_OEmodels/AZ/Benthic Data Bulk Upload_SOP.doc",
            "to import bug and habitat data, harmonize taxa list, rarefy and compute MMI",
            "then read resulting excel file back into R to save results in the database.", sep="\n"))
  return(AZbugs2)

}
