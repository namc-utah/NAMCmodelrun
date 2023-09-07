
rm(list=ls())
library(randomForest)
library(NAMCr)
library(tidyverse)
library(devtools)
library(BMIMetrics)
library(CSCI)
library(factoextra)
#read in model app

source("R/ModelApplicabilityAll.R")
#read in calibration of model app
ecoregion_base_path="C://Users//andrew.caudillo//Box//NAMC//"
CalPredsModelApplicability=read.csv(paste0(ecoregion_base_path,"OEModeling/NAMC_Supported_OEModels/Model Applicability/CalPredsModelApplicability.csv"))
#read in the ADEQ IBI results
#remember that this needs to be calculated using the ADEQ IBI Access database,
#much like the CO EDAS.
AZ_IBI_results<-read.csv('C://Users//andrew.caudillo//Box//NAMC//OEModeling//NAMC_Supported_OEmodels//Arizona//Results//AZAIM2022_results.csv')
#calculate applicability preds
applicabilitypreds = NAMCr::query("samplePredictorValues",
                                  include = c(
                                    "sampleId",
                                    "predictorId",
                                    "status",
                                    "abbreviation",
                                    "predictorValue"
                                  ),
                                  sampleIds = AZ_IBI_results$ActivityID
) #need list of samples in database with values
applicabilitypreds = subset(applicabilitypreds, abbreviation %in% c('ElevCat','Tmean8110Ws','WsAreaSqKm','Precip8110Ws'))
applicabilitypreds$predictorValue=as.numeric(applicabilitypreds$predictorValue)
applicabilitypreds = tidyr::pivot_wider(applicabilitypreds,
                                        id_cols="sampleId",
                                        names_from = "abbreviation",
                                        values_from = "predictorValue")# add id_cols=sampleId once it gets added to end point
applicabilitypreds=as.data.frame(applicabilitypreds)

#make the model app for your samples
#maybe can make this more robust for cold vs warm water?
#but most are warm (less than 5000 feet in elevation)
#so no need to do this now as of 8/31/23
ModelApplicability_obj = ModelApplicability(CalPredsModelApplicability,
                                            modelId =236 ,
                                            applicabilitypreds)
#get invasive info
bugsOTU = NAMCr::query("sampleTaxaTranslationRarefied",
                       translationId = 26,
                       fixedCount = 500,
                       sampleIds=AZ_IBI_results$ActivityID)
sumrarefiedOTUTaxa = bugsOTU  %>%
  dplyr::group_by(sampleId) %>%
  dplyr::summarize(fixedCount = sum(splitCount))

################################################
###### get invasives ##### comment out this entire invasives section if running "National" AIM westwide reporting
# get raw bug data
bugRaw = NAMCr::query(
  "sampleTaxa",
  sampleIds=AZ_IBI_results$ActivityID
)
#subset taxa in samples to only invasives
bugraw = subset(bugRaw,taxonomyId %in% c(1330,1331,2633, 2671,4933,4934,4935,4936,4937,4938,4939,4940,4941,4942,1019,1994,5096,1515,1518,1604,2000,4074,1369,2013,1579))
#create list of invasives present at a site
invasives<-bugraw %>% dplyr::group_by(sampleId) %>% dplyr::summarize(InvasiveInvertSpecies=paste0(list(unique(scientificName)),collapse=''))
# remove list formatting
invasives$InvasiveInvertSpecies=gsub("^c()","",invasives$InvasiveInvertSpecies)
invasives$InvasiveInvertSpecies=gsub("\"","",invasives$InvasiveInvertSpecies)
invasives$InvasiveInvertSpecies=gsub("\\(","",invasives$InvasiveInvertSpecies)
invasives$InvasiveInvertSpecies=gsub("\\)","",invasives$InvasiveInvertSpecies)
# join to list of all samples with fixed counts
additionalbugmetrics=dplyr::left_join(sumrarefiedOTUTaxa,invasives, by="sampleId")
# if no invasives were present set to absent
additionalbugmetrics[is.na(additionalbugmetrics)]<-"Absent"
data.table::setnames(AZ_IBI_results,'ActivityID','sampleId')
#merge IBI with supplementary data (invasives, model app)
AZ_IBI_results<-merge(AZ_IBI_results,additionalbugmetrics,by='sampleId')
AZ_IBI_results<-merge(AZ_IBI_results,ModelApplicability_obj,by='sampleId')

#for loop to save the reults to the database
for(i in 1:nrow(AZ_IBI_results)){
dat_to_pass<-list(
  sampleId=AZ_IBI_results$sampleId[i],
  modelId=ifelse(AZ_IBI_results$InvertReg[i]=='warm',236,169),
  modelResult=AZ_IBI_results$IBI[i],
  fixedCount=AZ_IBI_results$TotalIndividualsRaw[i],
  modelApplicability=AZ_IBI_results$ModelApplicability[i],
  notes=AZ_IBI_results$InvasiveInvertSpecies[i])

NAMCr::save(
  api_endpoint = "setModelResult",
  args=dat_to_pass)

}

#Great! now you have saved the results.
#go to GeneralCode/NAMCSimpleReport and run the simple report
#after altering the Config file accordingly
#and send the results to Bill Wells and Janet Miller.
