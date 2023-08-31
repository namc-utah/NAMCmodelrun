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
    include = c("boxId","sampleId",'siteId','sampleDate',"siteName", "waterbodyName","sampleMethod","habitatName","area", "labSplit"),# possibly add waterbody name to the samples NAMCr::query
    sampleIds=sampleIds
  )

  # join that data together into a single dataframe
  AZbugs=dplyr::left_join(AZsubsamp,bugsTranslation, by=c("taxonomyId", "sampleId"))
  AZbugs=dplyr::left_join(AZsubsamp,samples, by='sampleId')

  AZbugs$StationID=AZbugs$siteName
  AZbugs$WaterbodyName=AZbugs$waterbodyName
  AZbugs$ActivityID=AZbugs$sampleId
  AZbugs$RepNum=1
  AZbugs$CommentsSample=NA
  AZbugs$CollDate=AZbugs$sampleDate
  AZbugs$CollMeth='ADEQ Riffle bugs'
  AZbugs$CorrectionFactor=100*(1/AZbugs$labSplit)
  AZbugs$FinalID=AZbugs$scientificName
  AZbugs$Individuals=AZbugs$splitCount
  AZbugs$Stage=AZbugs$lifeStageAbbreviation
  AZbugs$LargeRare='No'
  AZbugs$Habitat=ifelse(AZbugs$habitatName=='Targeted Riffle','Riffle','Multi-Habitat')
  AZbugs$Lab='NAMC'
  AZbugs$LabID=NA

  AZbugs2=AZbugs[,c("StationID","waterbodyName","ActivityID","RepNum","CollDate","CommentsSample","CollMeth",
                    "CorrectionFactor","FinalID","Individuals","Stage","LargeRare","Habitat","Lab","LabID")]
  #write excel file to workspace
  write.csv(AZbugs2,file = paste0("AZbugs","boxId_",boxId,"_",Sys.Date(),".csv"),row.names=FALSE)
  cat(paste("csv with AZbugs has been written out to your current working directory.",
            "convert this csv to excel 2003 and import into AZ EDAS access database to compute the IBI score.",
            "/OE_Modeling/NAMC_Supported_OEmodels/AZ/Benthic Data Bulk Upload_SOP.doc",
            "to import bug and habitat data, harmonize taxa list, rarefy and compute MMI",
            "then read resulting excel file back into R to save results in the database.", sep="\n"))
  return(AZbugs)
}
