# Load all needed packages
rm(list=ls())
#setwd('C://Users//andrew.caudillo//Documents//R_code//NAMCmodelrun')
library(randomForest)
library(NAMCr)
library(tidyverse)
library(devtools)
library(BMIMetrics)
library(CSCI)
library(factoextra)


# Load all needed functions
source("R/Bug_functions_sample.R")
source("R/Model_functions.R")
source("R/model.predict.RanFor.4.2.r")
source("R/model.predict.v4.1.r")
source("R/ModelApplicabilityAll.R")
#source("R/ADEQBugDataFormatter.R")

ecoregion_base_path="C:/Users/jenni/Box/NAMC (Trip Armstrong)/"
#ecoregion_base_path="C://Users//andrew.caudillo//Box//NAMC//"
# to run the scripts choose a boxId or projectId
boxId=8051
# test boxes
# 2141-UT,2065 OR WCCP and MCCP, null, 2152 PIBO, 2172 CSCI,2107 AREMP, 2054 CO, NV 2140, 1603 westwide and 2055,	2150 WY


#projects=NAMCr::query("projects")
#projectId=4182



# then input a modelID
models=NAMCr::query("models")

modelID=25

menu(c("Yes","No"),
     title=paste("Is this correct? Box number = ",boxId," Model(s) is/are ",modelID))
# all model results are always calculated
# overwrite controls which are saved in the database
#input overwrite='Y' if you want to replace existing model results in the database. Generally set this to 'N' unless you find an error you are trying to fix
overwrite='N'

# input file path for reference sites attributed with elevation, watershed area, and temperature for model applicability function
#remove this line once reference sites are all in database with stream cat data
CalPredsModelApplicability=read.csv(paste0(ecoregion_base_path,"OEModeling/NAMC_Supported_OEModels/Model Applicability/CalPredsModelApplicability.csv"))

#only need to run this part of CONFIG if AIM ID is the client

#This is ONLY for AIM Idaho
#The BLM has requested that we assign certain models
#to critical salmonid habitat, and another model for non-critical habitat.
ID_lookup<-read.csv(paste0(ecoregion_base_path,"OEModeling/NAMC_Supported_OEModels/Model Applicability/BLM_ID_COMID_metadata.csv"),stringsAsFactors = F)
names(ID_lookup)[2]<-'waterbodyCode'
ID_lookup<-ID_lookup[!duplicated(ID_lookup$waterbodyCode),]
#crit is 9, not crit 25

ID_lookup$modelId<-ifelse(ID_lookup$CritHab=='No',25,9)
IDsites= query(
  api_endpoint = "sites",
  args = list(boxIds = boxId) #change box number as needed
)
ID_lookup_table<-plyr::join(IDsites[,c('waterbodyCode','siteId')],ID_lookup[,c('waterbodyCode','CritHab','modelId')],by='waterbodyCode',type='inner')
CritHab<-ID_lookup_table[ID_lookup_table$modelId==9,]
NonCrit<-ID_lookup_table[ID_lookup_table$modelId==25,]



