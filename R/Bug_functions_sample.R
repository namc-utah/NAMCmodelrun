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

#' OR MMI  bug data export
#'
#' @param sampleIds
#'
#' @return raw bug data translated by OR translation and with column names formatted for BioMonTools metric function
#' @export
#'
#' @examples
OR_MMI_bug_export <- function(sampleIds){
  # get needed data from the APIs
  bugRaw = NAMCr::query(
    "sampleTaxaTranslationRarefied",
    sampleIds = sampleIds,
    translationId = translationId,
    fixedCount = fixedCount
  )# raw NAMCr::query with pivoted taxonomy, and join translation name but not roll it up.... then summ in here

  bugsTranslation = NAMCr::query(
    "sampleTaxaTranslation",
    translationId = 9,####edit
    sampleIds = sampleIds
  )
  sites = NAMCr::query(
    "samples",
    include = c("sampleId", 'siteName'),
    sampleIds = sampleIds
  )
  # join that data together into a single dataframe
  OR_MMIbugs=dplyr::left_join(bugRaw,bugsTranslation, by=c("taxonomyId", "sampleId"))
  OR_MMIbugs=dplyr::left_join(OR_MMIbugs,sites, by='sampleId')

  #rename columns for BiomonTools functions
   bug_tax_nhd <- get_NHD_info(df_bugs_taxa) |>
    dplyr::transmute(SampleID = Sample,
                     Area_mi2 = NA_integer_,
                     SurfaceArea = NA_integer_,
                     TaxaID = TAXAID,
                     N_Taxa = Count,
                     Index_Name ='MMI_metrics',
                     INDEX_CLASS = str_to_title(SITE_TYPE),
                     NonTarget,
                     SITE_TYPE,
                     Kingdom,
                     SubOrder,
                     SubFamily,
                     GenusGroup,
                     SpeciesGroup,
                     SpeciesSubGroup,
                     SpeciesComplex,
                     Phylum,
                     SubPhylum,
                     Class,
                     SubClass,
                     Order,
                     SuperFamily,
                     Family,
                     Tribe,
                     Genus,
                     SubGenus,
                     Species,
                     BCG_Attr  = BCG_attr,
                     FFG,
                     Habit,
                     Life_Cycle,
                     Thermal_Indicator = Thermal_indicator,
                     TolVal,
                     INFRAORDER = NA_character_,
                     HABITAT = Habitat,
                     ELEVATION_ATTR = NA_character_,
                     GRADIENT_ATTR = NA_character_,
                     WSAREA_ATTR = NA_character_,
                     HABSTRUCT = NA_character_,
                     UFC = NA_integer_,
                     Density_ft2. = NA_integer_,
                     DENSITY_M2 = NA_integer_,

    )

  #use biomontools markExcluded function to generate the exclude column
  bugs.excluded <- BioMonTools::markExcluded(bug_tax_nhd, TaxaLevels = c("Kingdom", "Phylum",
                                                                         "SubPhylum", "Class", "SubClass", "Order", "SubOrder", "SuperFamily",
                                                                         "Family", "SubFamily", "Tribe", "GenusGroup", "Genus", "SubGenus", "SpeciesGroup",
                                                                         "SpeciesSubGroup", "SpeciesComplex", "Species"))

  mets.keep <- c('pt_tv_intol', 'nt_habitat_rheo', 'pt_ti_stenocold_cold_cool', 'pi_EPTNoHydro')


  #Calculate metrics
  # The boo.Shiny argument prevents the code from stopping and asking permission to calculate metrics with missing parameters.
  # It will give a warning instead.

  # The function will return the metrics identified in mets.keep
  metricsdf <- BioMonTools::metric.values(bugs.excluded, "bugs",fun.MetricNames = mets.keep, boo.Shiny	= TRUE)


  bugnew= metricsdf
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
    "sampleTaxaTranslation",
    sampleIds=sampleIds,
    translationId=8
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
      "sampleTaxaTranslation",
      sampleIds=sampleIds,
      translationId=92,
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
  cat(paste("csv with CO bugs has been written out to the CO EDAS (2025) imports folder.",
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
AZ_bug_export<-function(sampleIds,AZ_traits){
  AZ_Traits=read.csv('C://Users//andrew.caudillo.BUGLAB-I9//Box//NAMC//OEModeling//NAMC_Supported_OEmodels//Arizona//Traits_and_taxa_NAMCcrosswalk_complete_clean.csv')
  #mats_taxa=read.csv("inputs/mats_taxa_species.csv")
  #load(file = "inputs/mats_site.rdata")

  # FFG and TOLVAL are by OTU not by taxa.  Add FFG and TOLVAL by OTU.

  # Correction for Operational Taxa Unit.  Ensures TolVal and FFG uses appropriate level from taxa table.  Note Taxa table has a lot of null values for tolerance and FFG.  Keep tolval and FFG if no OTU
  #getting trait info for OTU
  mats_taxa_otu <- AZ_Traits %>%
    mutate(filt = ifelse(TAXA_ID == OTU_TAXA_ID, "Y", "N")) %>% # Makes the OTU record the primary for FFG and TOLVAL.
    filter(filt == "Y") %>%
    distinct(OTU_TAXA_ID, TOLVAL, FFG)

  # Includes rolled up tolval from mats_taxa_otu (see chironomidae 42 all being tolval 6).  Also includes OTU that are missing but has tolvalues
  #applying OTU trait data to lower level taxa
  mats_taxa <- AZ_Traits %>%
    select(TAXA_ID, OTU_TAXA_ID, OTU_ADEQ, T_ORDER, FAMILY, GENUS, IBI_TAXA_STATUS_CD,NAMC_taxonomy_id, TOLVAL2 = TOLVAL, FFG2 = FFG) %>%
    left_join(mats_taxa_otu, by = c("OTU_TAXA_ID")) %>%
    mutate(TOLVAL = ifelse(is.na(TOLVAL) & !is.na(TOLVAL2), TOLVAL2, TOLVAL)) %>%
    mutate(FFG = ifelse(is.na(FFG) & !is.na(FFG2), FFG2, FFG)) %>%
    select(-TOLVAL2, -FFG2)

##get NAMC bug data
  bugRaw = NAMCr::query(
    "sampleTaxa",
    sampleIds=sampleIds
  )
  # sum taxa across all lifestages
  sumTaxa = bugRaw  %>%
    dplyr::group_by(sampleId, taxonomyId,scientificName) %>%
    dplyr::summarize(sumSplitCount = sum(splitCount))

  sumTaxa$NAMC_taxonomy_id= sumTaxa$taxonomyId
  sumTaxa=sumTaxa[,names(sumTaxa) != 'taxonomyId']

#join AZ trait data into NAMC bug data
  mats_individuals=left_join(sumTaxa,mats_taxa, by = "NAMC_taxonomy_id")%>%
    select(SampleID=sampleId,
           NAMC_taxonomy_id,
           scientificName,
           Individuals = sumSplitCount,
           OTU_ADEQ = OTU_ADEQ,
           Order = T_ORDER,
           Family = FAMILY,
           Genus = GENUS,
           TolVal = TOLVAL,
           FFG)

  # creating a column that has the heiarchy equivalent to a scientific name column
  # Join individuals and taxa
  mats_raw <- mats_individuals %>%
    mutate(Mark = "Y") %>%
    mutate(phylo = paste(Order, Family, Genus))


  ### using this dataset for richness metrics and filtering only taxa with OTUs
  # Ensure that one taxa per sample for multiple levels (family/order).  OTU null values use determined per metric.
  mats_no_na <- mats_raw %>%
    filter(!is.na(OTU_ADEQ)) %>% # we do not use taxa with no OTU
    group_by(SampleID, OTU_ADEQ) %>%
    mutate(Individuals = sum(Individuals, na.rm = TRUE)) %>%
    distinct(SampleID, OTU_ADEQ, .keep_all = TRUE) %>%
    ungroup()

  mats_just_fam <- mats_raw %>%
    filter(is.na(OTU_ADEQ)) %>%
    mutate(lowest = ifelse(!is.na(Family), "Family",
                           ifelse(!is.na(Order), "Order", "Reject"))) %>%
    filter(lowest == "Family") %>%
    group_by(SampleID, Family) %>% # add family
    mutate(Individuals = sum(Individuals)) %>%
    distinct(SampleID, Family, .keep_all = TRUE) %>% # added because lose if just distinct on OTU w/ NAs
    ungroup() %>%
    select(-lowest)

  mats_just_ord <- mats_raw %>%
    filter(is.na(OTU_ADEQ)) %>%
    mutate(lowest = ifelse(!is.na(Family), "Family",
                           ifelse(!is.na(Order), "Order", "Reject"))) %>%
    filter(lowest == "Order") %>%
    group_by(SampleID, Order) %>% # add order
    mutate(Individuals = sum(Individuals)) %>%
    distinct(SampleID, Order, .keep_all = TRUE) %>% # added because lose if just distinct on OTU w/ NAs
    ungroup() %>%
    select(-lowest)

  mats <- bind_rows(mats_no_na, mats_just_fam, mats_just_ord)


  # Figure out which NA's for OTU_ADEQ are present at order or family
  # Step 1 - Determine lowest level identified ----
  mats_lowest <- mats %>%
    filter(is.na(OTU_ADEQ)) %>%
    mutate(lowest = ifelse(!is.na(Family), "Family",
                           ifelse(!is.na(Order), "Order", "Reject"))) %>%
    distinct(SampleID, phylo, lowest,Mark)

  # Step 2 - Identify rows that already have taxa present at level 2 ----
  mats_exclude <- mats_lowest %>%
    left_join(mats, by = c("SampleID", "lowest" = "Family")) %>%
    mutate(family_flag = Mark.x) %>%
    left_join(mats, by = c("SampleID", "lowest" = "Order")) %>%
    mutate(order_flag = Mark.y) %>%
    dplyr::group_by(SampleID, lowest) %>% # look for multiples.  Determine situations where taxonomist able to identify some taxa to genus and some to family IN SAME GROUP (ex simulidae only and simulidae and genus should not count as 2 genera)
    dplyr::mutate(count = n()) %>% # 2's with a na in oTU = exclude flag
    mutate(OTU_Temp = ifelse(!is.na(OTU_ADEQ.x), OTU_ADEQ.x,
                             ifelse(!is.na(OTU_ADEQ.y), OTU_ADEQ.y, NA))) %>%
    mutate(exclude = ifelse(count > 1 & is.na(OTU_Temp), "Y", "N")) %>% # exclude flag for taxa not identified to lowest taxa for metrics like number of taxa but have valid information for other metrics like percent scraper.
    filter(exclude == "Y") %>%
    ungroup() %>%
    select(SampleID, phylo = phylo.x, exclude)

  # Add back to original
  mats <- mats %>%
    left_join(mats_exclude, by = c("SampleID", "phylo"))


    return(mats)
  # unam=NAMCr::query('sampleTaxaUnambiguous',boxId=3793)
  # mats$scientificName[mats$scientificName %in% unam$scientificName==F]
  # unam$scientificName[unam$scientificName %in% mats$scientificName==F]
  # unique(unam$scientificName)[unique(unam$scientificName) %in% unique(mats$scientificName)==F]
}




