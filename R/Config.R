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
boxId=2190
# test boxes
# 2141-UT,2065 OR WCCP and MCCP, null, 2152 PIBO, 2172 CSCI,2107 AREMP, 2054 CO, NV 2140, 1603 westwide and 2055,	2150 WY

projects=NAMCr::query("projects")
#projectId=


# then input a modelID
models=NAMCr::query("models")
modelID=136


# all model results are always calculated
# overwrite controls which are saved in the database
#input overwrite='Y' if you want to replace existing model results in the database. Generally set this to 'N' unless you find an error you are trying to fix
overwrite='N'

# input file path for reference sites attributed with elevation, watershed area, and temperature for model applicability function
#remove this line once reference sites are all in database with stream cat data
CalPredsModelApplicability=read.csv("C:/Users/jenni/Box/NAMC (Trip Armstrong)/OE_Modeling/NAMC_Supported_OEModels/Model Applicability/CalPredsModelApplicability.csv")


