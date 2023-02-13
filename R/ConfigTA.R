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

##TWA office
ecoregion_base_path="/Users/triparmstrong/Library/CloudStorage/Box-Box/NAMC/"
CalPredsModelApplicability=read.csv("/Users/triparmstrong/Library/CloudStorage/Box-Box/NAMC/OE_Modeling/NAMC_Supported_OEmodels/Model Applicability/CalPredsModelApplicability.csv")

##TWA home
ecoregion_base_path="/Users/namc/Library/CloudStorage/Box-Box/NAMC/"
CalPredsModelApplicability=read.csv("/Users/namc/Library/CloudStorage/Box-Box/NAMC/OE_Modeling/NAMC_Supported_OEmodels/Model Applicability/CalPredsModelApplicability.csv")

# to run the scripts choose a boxId or projectId
boxId=4816

projects=NAMCr::query("projects")
#projectId=1641

# then input a modelID
models=NAMCr::query("models")

#modelID=25

modelID = c(4)



# all model results are always calculated
# overwrite controls which are saved in the database
#input overwrite='Y' if you want to replace existing model results in the database. Generally set this to 'N' unless you find an error you are trying to fix
overwrite='Y'

# input file path for reference sites attributed with elevation, watershed area, and temperature for model applicability function
#remove this line once reference sites are all in database with stream cat data




