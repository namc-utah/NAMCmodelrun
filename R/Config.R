# # Load all needed packages
# library(randomForest)
# library(NAMCr)
# library(tidyverse)
# library(CSCI)
#
#
# # Load all needed functions
# source("R/Bug_functions_sample.R")
# source("R/Model_functions.R")
# source("R/model.predict.RanFor.4.2.r")
# source("R/model.predict.v4.1.r")
# source("R/ModelApplicabilityAll.R")

library(randomForest)
library(NAMCr)
library(tidyverse)
library(DBI)
library(RSQLite)
library(nhdplusTools)
library(devtools)
library(BMIMetrics)
library(CSCI)
library("factoextra")
boxId=2141# 2141-UT,2065 OR WCCP and MCCP, null, 2152 PIBO, 2172 CSCI,2107 AREMP, 2054 CO, NV 2140, 1603 westwide and 2055,	2150 WY
projects=NAMCr::query("projects")
#projectId=
models=NAMCr::query("models")
modelId=8
modelId=3
#prednew=read.csv("C:/Users/jenni/Box/NAMC (Trip Armstrong)/OE_Modeling/NAMC_Supported_OEmodels/Wyoming/InputsAndResults/CurrentRun/WY_Habitat.csv", row.names="SampleID")
#prednew=read.csv("C:/Users/jenni/Box/NAMC (Trip Armstrong)/OE_Modeling/NAMC_Supported_OEmodels/West_Wide_Model/NewWestWIde_NONMidgeModel/InputsAndResults/AIM_NM_MT_2019/WW18_Preds.csv", row.names="SampleID")
#prednew=read.csv("C:/Users/jenni/Box/NAMC (Trip Armstrong)/OE_Modeling/NAMC_Supported_OEmodels/West_Wide_Model/NewWestWIde_NONMidgeModel/InputsAndResults/CurrentRun/WW18_Preds.csv", row.names="SampleID")
#prednew=read.csv("C:/Users/jenni/Box/NAMC (Trip Armstrong)/OE_Modeling/NAMC_Supported_OEmodels/AREMP2014/InputsAndResults/2020/Test_preds.csv",row.names="SampleID")
#prednew=read.csv("C:/Users/jenni/Box/NAMC (Trip Armstrong)/OE_Modeling/NAMC_Supported_OEmodels/CA/Hybrid_CA_Model/InputsAndResults/CurrentRun/habitat.csv")
#prednew=read.csv("C:/Users/jenni/Box/NAMC (Trip Armstrong)/OE_Modeling/NAMC_Supported_OEmodels/PIBO/InputsAndResults_PIBO2009oe/AIM_ID_2020/Habitat.csv",row.names = "SAMPLEID")
prednew=read.csv("C:/Users/jenni/Box/NAMC (Trip Armstrong)/OE_Modeling/NAMC_Supported_OEmodels/UTDEQ/AllSeasonsModel_2015/InputsAndResults/AIM2020/UTDEQ_Habitat.csv",row.names="SAMPLE")
#prednew=read.csv("C:/Users/jenni/Box/NAMC (Trip Armstrong)/OE_Modeling/NAMC_Supported_OEmodels/OR/InputsAndResults_PredatorORDEQ2005oe/AIM_OR_2019/MWCF/MWCF_Habitat.csv",row.names="SampleID")
#prednew=read.csv("C:/Users/jenni/Box/NAMC (Trip Armstrong)/OE_Modeling/NAMC_Supported_OEmodels/OR/InputsAndResults_PredatorORDEQ2005oe/AIM_OR_2019/WCCP/WCCP_Habitat.csv",row.names="SampleID")
#prednew=read.csv("C:/Users/jenni/Box/NAMC (Trip Armstrong)/OE_Modeling/NAMC_Supported_OEmodels/NV/NV_MMI/InputsAndResults/CurrentRun/NV_Habitat.csv",row.names="SampleID")
SQLite_file_path="C:/NAMC_S3/StreamCat/StreamCat2022.sqlite"
CalPredsModelApplicability=read.csv("C:/Users/jenni/Box/NAMC (Trip Armstrong)/OE_Modeling/NAMC_Supported_OEModels/Model Applicability/CalPredsModelApplicability.csv")




#29-37
#254-end

source("R/Bug_functions_box.R")
source("R/Bug_functions_sample.R")
source("R/Model_functions.R")
source("R/model.predict.RanFor.4.2.r")
source("R/model.predict.v4.1.r")
source("R/ModelApplicabilityAll.R")
