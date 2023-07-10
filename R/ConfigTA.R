
rm(list=ls())


# Load all needed packages
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



# to run the scripts choose a boxId or projectId
boxId=3232
# test boxes
# 2141-UT,2065 OR WCCP and MCCP, null, 2152 PIBO, 2172 CSCI,2107 AREMP, 2054 CO, NV 2140, 1603 westwide and 2055,	2150 WY

projects=NAMCr::query("projects")
#projectId=2301


# then input a modelID
models=NAMCr::query("models")
modelID=c(2)

##TWA office
ecoregion_base_path="/Users/triparmstrong/Library/CloudStorage/Box-Box/NAMC/"
CalPredsModelApplicability=read.csv("/Users/triparmstrong/Library/CloudStorage/Box-Box/NAMC/OEModeling/NAMC_Supported_OEmodels/Model Applicability/CalPredsModelApplicability.csv")

##TWA home
#ecoregion_base_path="/Users/namc/Library/CloudStorage/Box-Box/NAMC/"
#CalPredsModelApplicability=read.csv("/Users/namc/Library/CloudStorage/Box-Box/NAMC/OE_Modeling/NAMC_Supported_OEmodels/Model Applicability/CalPredsModelApplicability.csv")


# all model results are always calculated
# overwrite controls which are saved in the database
#input overwrite='Y' if you want to replace existing model results in the database. Generally set this to 'N' unless you find an error you are trying to fix

overwrite='N'




