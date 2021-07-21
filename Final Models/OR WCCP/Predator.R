############################

# From JVS original code:
# The 2 files (bugs and predictors) should have similar formats as the original taxa and predictor data sets used to build the model (see step 1 above);
# Notes on format --
#   A) The sample ID column in both files should be read into R as a row name (see Step 1 examples).
#   B) Predictor data set -- Must include columns with the same names, units, etc.,
#        as the model's predictor variables. All other columns will be ignored;
#        Column order does not matter;
#        Predictions and calculations of O/E will be made only for those samples that have;
#        complete data for all model predictors.;
#   C)  Sample-by-taxa matrix.  Sample ID's (row names) must match those of predictor data.
#        Model predictions will include only those taxa (columns) in this file that are also
#        in the calibration data (bugcal.pa);
#        To see a list of the calibration taxa, do the following:
#                                                         names(bugcal.pa)[colSums(bugcal.pa)>0];
#  Code modified by TWA (January 2019) in order to utilize NAMC file structure and output needs
##########;

#source in prediction script and load models
#source("//deqlab1/biomon/R Stats/Bio Tools_Upgrade with R/model.predict.v4.1.r") # bring in model predict function (JVS script)

# IF MWCF use this directory
setwd("\\\\share1.bluezone.usu.edu\\miller\\\\buglab\\OE_Modeling\\NAMC_Supported_OEmodels\\OR\\ModelFiles_PredatorORDEQ2005oe\\MWCF\\")
load("MWCF.rdata")
library(plyr)
source("model.predict.v4.1.R")
Boxdata <- ("\\\\share1.bluezone.usu.edu\\miller\\\\buglab\\OE_Modeling\\NAMC_Supported_OEmodels\\OR\\InputsAndResults_PredatorORDEQ2005oe\\CurrentRun")

###  Models can be run together all in single command; second section (WCCP) starts by clearing MWCF files

#################

# 5.1 MWCF model

################

load('Nov05model_MWCF_16jan13.Rdata');    # bring in MWCF model


#Drop all samples/sites from bug and predictor data that do not not have complete data for the model predictors;
#preds.mwcf.F<-rename(preds.mwcf.F, c("Long"="long"))

setwd(Boxdata)
getwd()
preds.mwcf.F <- read.csv("MWCF_Habitat.csv",row.names="SampleID",header=T)
bugs.MWCF.F <- read.csv("MWCF_Bugs.csv",row.names="SampleID",header=T)

preds.mwcf.F<-preds.mwcf.F[complete.cases(preds.mwcf.F[,preds.final]),]
bugs.MWCF.F<-bugs.MWCF.F[row.names(preds.mwcf.F),]
dim(preds.mwcf.F); dim(bugs.MWCF.F)   



#make predictions for test data;
OE.assess.test<-model.predict.v4.1(bugcal.pa,grps.final,preds.final, grpmns,covpinv,
                                   prednew=preds.mwcf.F,bugnew=bugs.MWCF.F,Pc=0.5);



oe.mwcf<-OE.assess.test$OE.scores; #create a d.f out of OE.scores
head(OE.assess.test$OE.scores)# look at O/E scores, for all samples;

head(OE.assess.test$Group.Occurrence.Probs) # Look at the predicted group membership probabilities;
head(OE.assess.test$Capture.Probs)  # Look at the predicted probabilties for each taxon at each site;

# Assign PREDATOR condition classes == MWCF benchmarks
oe.mwcf$OoverE<-round(oe.mwcf$OoverE, 2)
oe.mwcf$Condition<-(ifelse(oe.mwcf$OoverE <= 0.85, 'Most disturbed', 
                         ifelse(oe.mwcf$OoverE > 0.85 & oe.mwcf$OoverE < 0.92, 'Moderately disturbed', 
                                ifelse(oe.mwcf$OoverE >= 0.92 & oe.mwcf$OoverE < 1.25, 'Least disturbed', 
                                       ifelse(oe.mwcf$OoverE >= 1.25, 'Enriched', -999)))))


# calculate min - max for each condition class
ddply(oe.mwcf, .(oe.cond), summarize, min = min(OoverE), max = max(OoverE))
# verify that results are consistent with PREDATOR documentation benchmarks: <=0.85, 0.86 - 0.91, 0.92 - 1.24, > 1.24
# verify no '-999' values

## Write  tables of O/E outputs

write.csv(oe.mwcf,file = "MWFC_OEscores.csv")
write.csv( OE.assess.test$Group.Occurrence.Probs, "MWFC_grp.probs.csv", row.names=TRUE)
write.csv( OE.assess.test$Capture.Probs, "MWFC_capt.probs.csv", row.names=TRUE)




#############

# 5.2  WCCP model

#############

# need to remove MWCF objects, or WCCP predictions will be made on MWCF model data
rm(bugcal.pa,grps.final,preds.final, grpmns,covpinv, predcal) 



setwd("\\\\share1.bluezone.usu.edu\\miller\\\\buglab\\OE_Modeling\\NAMC_Supported_OEmodels\\OR\\ModelFiles_PredatorORDEQ2005oe\\WCCP\\")
load("WCCP.RData")
library(plyr)
source("model.predict.v4.1.R")
Boxdata <- ("\\\\share1.bluezone.usu.edu\\miller\\\\buglab\\OE_Modeling\\NAMC_Supported_OEmodels\\OR\\InputsAndResults_PredatorORDEQ2005oe\\CurrentRun")

# bring in WCCP model
load('Nov05model_WCCP_16jan13.RData')


setwd(Boxdata)
getwd()
preds.wccp.F <- read.csv("WCCP_Habitat.csv",row.names="SampleID",header=T)
bugs.wccp.F <- read.csv("WCCP_Bugs.csv",row.names="SampleID",header=T)

#Drop all samples/sites from bug and predictor data that do not not have complete data for the model predictors;
# preds.wccp.F<-rename(preds.wccp.F, c("Long"="long"))

preds.wccp.F<-preds.wccp.F[complete.cases(preds.wccp.F[,preds.final]),]
bugs.wccp.F<-bugs.wccp.F[row.names(preds.wccp.F),]
dim(preds.wccp.F); dim(bugs.wccp.F)     


#make predictions for test data;
OE.assess.test<-model.predict.v4.1(bugcal.pa,grps.final,preds.final, grpmns,covpinv,
                                   prednew=preds.wccp.F,bugnew=bugs.wccp.F,Pc=0.5);

oe.wccp<-OE.assess.test$OE.scores # create a d.f of OE.scores
head(oe.wccp)# look at O/E scores, for all samples;

head(OE.assess.test$Group.Occurrence.Probs)#Look at the predicted group membership probabilities;
head(OE.assess.test$Capture.Probs)#Look at the predicted probabilties for each taxon at each site;

# Assign PREDATOR condition classes == WCCP benchmarks
oe.wccp$OoverE<-round(oe.wccp$OoverE, 2)
oe.wccp$oe.cond<-(ifelse(oe.wccp$OoverE <= 0.78, 'Most disturbed', 
                         ifelse(oe.wccp$OoverE > 0.78 & oe.wccp$OoverE < 0.93, 'Moderately disturbed', 
                                ifelse(oe.wccp$OoverE >= 0.93 & oe.wccp$OoverE < 1.24, 'Least disturbed', 
                                       ifelse(oe.wccp$OoverE >= 1.24, 'Enriched', -999)))))

# calculate min - max for each condition class
ddply(oe.wccp, .(oe.cond), summarize, min = min(OoverE), max = max(OoverE))
# verify that results are consistent with PREDATOR documentation benchmarks: <=0.78, 0.79 - 0.92, 0.93 - 1.23, > 1.23


## Write  tables of O/E outputs
write.csv(oe.wccp,file = "WCCP_OEscores.csv")
write.csv( OE.assess.test$Group.Occurrence.Probs, "WCCP_grp.probs.csv", row.names=TRUE)
write.csv( OE.assess.test$Capture.Probs, "WCCP_capt.probs.csv", row.names=TRUE)

