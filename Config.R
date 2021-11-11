library(randomForest)
library(labdsv)
library(vegan)
library(cluster)
library(Hmisc)
library(plyr)
library(dplyr)
library(gtools)
library(MASS)
#install.packages("devtools")#Install devtools from CRAN
#library(devtools)
#install_github("SCCWRP/BMIMetrics")
#install_github("SCCWRP/CSCI")
library(CSCI)

# load main John Vansickle functions
source("model_predict_function/model.predict.v4.1.r")#OR
source("model_predict_function/model.predict.RanFor.4.2.r")#WW, AREMP,UT
source("model_predict_function/model.predict.RanFor.r")#PIBO

# R objects
#AREMP
load("final_models/AREMP/AREMP_OE_standardized.rdata")
load("final_models/AREMP/AREMP_MMI.rdata")
#NV
load("final_models/NV/OE_MMI_models.rdata")
#OR MWCF
load("final_models/OR_MWCF/MWCF.rdata")
#load('final_models/OR_MWCF/Nov05model_MWCF_16jan13.Rdata')
#OR WCCP
load("final_models/OR_WCCP/WCCP.RData")
#load('Nov05model_WCCP_16jan13.RData')
#PIBO
load("final_models/PIBO/Benkendorf.RF.Model.Version1_standardized.Rdata")
#UT
load('final_models/UTDEQ15/UTDEQ_15_OE_model_standardized.rdata')
#WY
load('final_models/WYDEQ/WYRIVPACS2012_standardized.Rdata')
#WW
load("final_models/WW18/My.RF.Model.Version1_standardized.Rdata")
