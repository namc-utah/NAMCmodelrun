
ADEQ_bug_export<-function(sampleIds){
  bugRaw = NAMCr::query(
    "sampleTaxaUnambiguous",
    sampleIds=sampleIds
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
  AZbugs=dplyr::left_join(AZsubsamp,bugsTranslation, by=c("taxonomyId", "sampleId"))
  AZbugs=dplyr::left_join(AZsubsamp,samples, by='sampleId')

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
  AZbugs2$FinalID=AZbugs$scientificName
  AZbugs2$Individuals=AZbugs$splitCount
  AZbugs2$Stage=AZbugs$lifeStageAbbreviation
  AZbugs2$LargeRare='No'
  AZbugs2$Habitat=ifelse(AZbugs$habitatName=='Targeted Riffle','Riffle','Multi-Habitat')
  AZbugs2$Lab='NAMC'
  AZbugs2$LabID=NA



AZbugs2

  #write excel file to workspace
  write.csv(AZbugs2,file = paste0("AZbugs","boxId_",boxId,"_",Sys.Date(),".csv"),row.names=FALSE)
  cat(paste("csv with AZbugs has been written out to your current working directory.",
            "convert this csv to excel 2003 and import into AZ EDAS access database to compute the IBI score.",
            "/OE_Modeling/NAMC_Supported_OEmodels/AZ/Benthic Data Bulk Upload_SOP.doc",
            "to import bug and habitat data, harmonize taxa list, rarefy and compute MMI",
            "then read resulting excel file back into R to save results in the database.", sep="\n"))
  return(AZbugs2)

}
