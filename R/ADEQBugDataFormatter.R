
ADEQ_bug_export<-function(sampleIds){
  library(elevatr)
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
    include = c("boxId","sampleId",'siteId','sampleDate',"siteName", "waterbodyName","sampleMethod","habitatName","area", "labSplit",
                'siteLatitude','siteLongitude'),
    sampleIds=sampleIds
  )
  sites = NAMCr::query(
    "sites",
    sampleIds=sampleIds
  )
  # join that data together into a single dataframe
  AZbugs=dplyr::left_join(AZsubsamp,bugsTranslation[,c('otuName','taxonomyId','sampleId')], by=c("taxonomyId", "sampleId"))

  #getting invertReg
  spatdf<-st_as_sf(samples,coords=c("siteLongitude",'siteLatitude'))

  #assing lat/long projection
  spatdf<-st_set_crs(spatdf$geometry,"+proj=longlat")
  #get the elevation (returns a df)

  #get elevation, returns a DF
  invertreg_values<-get_elev_point(spatdf)
  #convert meters to feet (we could keep meters and just make a new
  #conditional statement, but this is fine.)
  invertreg_values$elevation<-invertreg_values$elevation*3.28084; invertreg_values$elev_units<-"feet"

  #assign the cold or warm
  invertreg_values$Type=ifelse(invertreg_values$elevation>5000,"cold","warm")
  #append sampleId (helps to find cold vs warm sites
  #on INSTAR, since we need to add the values individually on EDAS)
  invertreg_values$sampleId<-samples$sampleId
  plain_invertreg<-as.data.frame(invertreg_values[,c('sampleId','Type')])
  plain_invertreg<-plain_invertreg[,1:2]
  #get site name, too. This is what we will use when
  #assigning InvertReg on EDAS, since the ActivityID is less intuitive there.
  #(i.e., filter by the site names below and then add "warm" or "cold" accordingly)
  samples<-plyr::join(samples,plain_invertreg, by='sampleId')
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
  #correction factor is 100 * 1/labsplit %. So we need to force decimal to %, which is why we have split *100.
  AZbugs2$CorrectionFactor=100*(1/(AZbugs$labSplit*100))
  AZbugs2$FinalID=AZbugs$otuName
  AZbugs2$Individuals=AZbugs$splitCount
  AZbugs2$Stage=AZbugs$lifeStageAbbreviation
  AZbugs2$LargeRare='No'
  AZbugs2$Habitat=ifelse(AZbugs$habitatName=='Targeted Riffle','Riffle','Multi-Habitat')
  AZbugs2$Lab='NAMC'
  AZbugs2$LabID=NA
  AZbugs2$Latitude=AZbugs$siteLatitude
  AZbugs2$Longitude=AZbugs$siteLongitude
  AZbugs2$InvertReg=AZbugs$Type
  AZbugs2<-AZbugs2[which(is.na(AZbugs2$FinalID)==F),]

  #write excel file to workspace

  write.csv(AZbugs2,file = paste0("C://Users//andrew.caudillo.BUGLAB-I9//Box//NAMC//OEModeling//NAMC_Supported_OEmodels//Arizona//InputFiles//AZbugs","boxId_",boxId,"_",Sys.Date(),".csv"),row.names=FALSE)
  cat(paste("csv with AZbugs has been written out to your current working directory.",
            "convert this csv to excel 2003 and import into AZ EDAS access database to compute the IBI score.",
            "/OE_Modeling/NAMC_Supported_OEmodels/AZ/Benthic Data Bulk Upload_SOP.doc",
            "to import bug and habitat data, harmonize taxa list, rarefy and compute MMI",
            "then read resulting excel file back into R to save results in the database.", sep="\n"))
  return(AZbugs2)

}

#ADEQ_bug_export(sampleIds = sampleIds)


