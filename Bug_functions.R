
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
}


MMI_metrics<-function(sampleId,fixedCount){
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
}




CSCI_bug <- function(sampleId, translationId){
  
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
  bugRaw = query(
    "sampleRawTaxa",
    include = c(),
    translationId = def_models$translationId,
    sampleIds = def_model_results$SampleId
  )# raw query with pivoted taxonomy, and join translation name but not roll it up.... then summ in here
  taxonomy = query("taxonomy hearchy", sampleId)
  # can be
  CSCIbugs = bugsRaw
  CSCIbugs$SampleID =
    CSCIbugs$StationCode = sites$siteId#### need to get site info!!!
  CSCIbugs$FinalID = CSCIbugs$OTUName
  CSCIbugs$LifeStageCode = ifelse(CSCIbugs$lifeStage)
  # #LifeStageCode = case
  # when class <> 'insecta' then 'X'
  # when class is null then 'X'
  # else Lifestage
  # end,
  CSCIbugs$Distinct == 0
  CSCIbugs$BAResult = sum(CSCIbugs$splitCount, CSCIbugs$bigRare)# have query add this as a column
}


CO_bug_export<-function(boxId){
  CObugs=NAMCr::query("sampleTaxaTranslation",
                      translationId=,
                      sampleIds=)
  sites =query("samples", 
               include= c("sampleId",'siteId'),
               sampleIds = def_model_results$SampleId
  )
  
  CObugs$Project
  CObugs$Station
  CObugs$Name
  CObugs$Location
  CObugs$CollDate
  CObugs$Organism
  CObugs$Individuals
  CObugs$Stage
  CObugs$CommentsTaxa
  CObugs$RepNum
  CObugs$Grids
  CObugs$CommentsSample
  CObugs$CommentsRep
}

#PREDATOR NBR pseudo-O/E
#(The PREDitctive Assessment Tool for ORegon)
#based on FORTRAN Subsample.exe
#code by Andrew Caudillo, Lab Manager of NAMC
#modified from code ideas proposed by David Fowler
#Essentially a null O/E model (poor model; ref O/E SD is 0.29)

rm(list=ls())
setwd('Z://buglab//OE_Modeling//NAMC_Supported_OEmodels//OR//NBR//PREDATOR_Code')
ref_tax<-c('Baetis',
           'Brachycentrus',
           'Chironominae',
           'Diphetor_hageni',
           'Epeorus',
           'Optioservus',
           'Orthocladiinae',
           'Rhyacophila',
           'Trombidiformes',
           'Zaitzevia') #list of 10 "reference" taxa
E<-7.56 #null E is always the same. Does not account for environmental factors
nCount = 300 #count to sub-sample each sample's data by
countField = "Density" #the field name in data holding the counts 
sampleField = "SampleID"# the field name in data holding the sample identifier 
dataField= "OTUName" #the field name where the bugs are
#the data. this is ideally a query from NAMCr
#Where translations have already happened (I assummed they have with these data) 
#be sure to have that "_" for Diphetor_hageni!
sample_taxa<-read.csv('NBR_Predator_test.csv',stringsAsFactors = F)

sample_taxa<-sample_taxa[,c(2,4,7)] #isolate relevant variables

names(sample_taxa)<-c('SampleID','OTUName','Density') #consistent names


data = sample_taxa #the dataframe of data desired to subsample
data[[countField]]= floor( data[[countField]]) #get only integers out, not densities

# Explode the dataframe so each record is replicated by the number represented in the density column
data = data [rep( seq( nrow( data ) ), data[[countField]]), names(data) != countField]

library(dplyr) #dplyr will help us group by site and OTU
options(dplyr.summarise.inform = FALSE)
subData_list<-list() #create empty list for later

#this loop says for as many sites we have
#subsample to 300 for the bugs
#if a site has less than 300 bugs, just give us all the counts
#since we cannot subsample more than the population.
#then group by site and OTU for the final dataframe

for(i in 1:length(unique(data$SampleID))){
  site<-unique(data[[sampleField]])[i]
  X<-data[data[[sampleField]]==site,]
  
  nSampleIndividuals = nrow(X)
  indiv<-if(nSampleIndividuals < nCount) {
    X$OTUName
  } else {
    sample(X$OTUName, size=nCount, replace=F)
  }
  
  subData = data.frame(SampleID = rep(site, length(indiv)),
                       OTUName = indiv)%>%
    group_by(OTUName,SampleID) %>%
    summarize("Density"=n()) 
  subData_list[[i]]<-subData
}

big_sub_data<-do.call(rbind,subData_list) #rbind the list into one df
PREDATOR_subdata<-as.data.frame(big_sub_data) #we don't want a tibble
PREDATOR_subdata
PREDATOR_results<-data.frame(
  O=rep(NA, length(unique(sample_taxa$SampleID))),
  E=rep(E, length(unique(sample_taxa$SampleID))),
  OE=rep(NA, length(unique(sample_taxa$SampleID))),
  OE_cond=rep(NA, length(unique(sample_taxa$SampleID))),
  Station=rep(NA,length(unique(sample_taxa$SampleID))),
  Count=rep(NA, length(unique(sample_taxa$SampleID)))) #empty df for final O/Es

#this for loop makes the O/E scores
#for every site with at least 1 taxon
for(j in 1:length(unique(PREDATOR_subdata$SampleID))){
  X<-PREDATOR_subdata[PREDATOR_subdata$SampleID==
                        unique(PREDATOR_subdata$SampleID)[j],]
  X<-X[X$OTUName %in% ref_tax,] #only want to pull out reference taxa
  PREDATOR_results$Station[j]<-unique(PREDATOR_subdata$SampleID)[j]
  O<-nrow(X)
  PREDATOR_results$O[j]<-O
  PREDATOR_results$OE[j]<-O/E
  PREDATOR_results$Count[j]<-sum(X$Density)
}

#pick up the sites with 0s
PREDATOR_results$O[is.na(PREDATOR_results$OE)]<-0
PREDATOR_results$OE[is.na(PREDATOR_results$OE)]<-0
PREDATOR_results$Count[is.na(PREDATOR_results$Count)]<-0
PREDATOR_results$Station[is.na(PREDATOR_results$Station)]<-
  unique(sample_taxa$SampleID)[unique(sample_taxa$SampleID) %in% 
                                 PREDATOR_results$Station==F] #replace NAs with SampleIDs

PREDATOR_results$OE_cond<-ifelse(
  PREDATOR_results$OE>0.71,
  "Good",
  ifelse(PREDATOR_results$OE <0.42,
         "Poor",'Fair'))
PREDATOR_results<-PREDATOR_results[order(PREDATOR_results$Station),] #sort by SampleID
row.names(PREDATOR_results)<-NULL #Reset row names so they are correct (incorrect after sorting)


