library(randomForest)
library(labdsv)
library(vegan)
library(cluster)
library(Hmisc)
library(plyr)
library(dplyr)# is this used?
library(gtools)
library(MASS)

library(tidyverse)
#install.packages("devtools")#Install devtools from CRAN
#library(devtools)
#install_github("SCCWRP/BMIMetrics")
#install_github("SCCWRP/CSCI")
library(CSCI)

# create a new package that has the functions below and load that as a dependency
# load main John Vansickle functions
source("model_predict_function/model.predict.v4.1.r")#OR
source("model_predict_function/model.predict.RanFor.4.2.r")#WW, AREMP,UT
source("model_predict_function/model.predict.RanFor.r")#PIBO



# fn_args=list(
#   model.predict.RanFor.4.2=list(
#     bugcal.pa="bugcal.pa",
#     grps.final="grps.final",
#     preds.final, 
#     ranfor.mod,
#     prednew,
#     bugnew,
#     Pc=0.5,
#     Cal.OOB=FALSE
#     
#   )
  