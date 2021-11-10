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
source("model.predict.v4.1.r")
source("model.predict.RanFor.4.2.r")#WW
source("\\\\share1.bluezone.usu.edu\\miller\\\\buglab\\OE_Modeling\\NAMC_Supported_OEmodels\\PIBO\\PIBO_RF2020\\model.predict.RanFor.r")

# R objects
#AREMP
load("AREMP_OE.rdata")
load("AREMP_MMI.rdata")
#NV
load("OE_MMI_models.rdata")
#OR MWCF
load("MWCF.rdata")
load('Nov05model_MWCF_16jan13.Rdata')
#OR WCCP
load("WCCP.RData")
load('Nov05model_WCCP_16jan13.RData')
#PIBO
load("Benkendorf.RF.Model.Version1.Rdata")
#UT
load('UTDEQ_15_OE_model.rdata')
#WY
load('WYRIVPACS2012.Rdata')
#WW
load("My.RF.Model.Version1.Rdata")