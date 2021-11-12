############################################################################################################################
# Quick overview:
# This script takes python-extracted predictor variables (i.e., EC12 and TNTP12) and runs them through 3 rf models resulting in predictions of EC, TN, and TP.
############################################################################################################################

############################################################################################################################
# More detailed summary:
# THE BIG IDEA IN THE NEXT CHUNK OF CODE IS THAT 3 RANDOM FOREST MODELS HAVE ALREADY BEEN CREATED.  THE 3 RESPONSE VARIABLES OF THESE MODELS ARE: 
#     1. ELECTRICAL CONDUCTANCE (EC) (CODE TALK FOR SPECIFIC CONDUCTANCE OR SALINITY), 
#     2. TOTAL PHOSPHOROUS (TP), AND 
#     3. TOTAL NITROGEN (TN). 
# THOSE MODELS HAVE BEEN SAVED OFF AND WE'LL RETRIEVE THEM LATER DOWN BELOW. 
#     IN THE MEANTIME, WE WANT TO PREPARE TO RUN THESE MODELS ON "NEW" AIM DATASETS.  
#     RECALL THAT THESE MODELS ARE BUILT USING REFERENCE SITE PREDICTORS.  
#         FOR EXAMPLE, THE EC MODEL REQUIRES N=19 PREDICTORS
#         THE TP MODEL REQUIRES N=15 PREDICTORS, AND
#         THE TN MODEL REQUIRES N=12 PREDICTORS.  
# OUR FIRST STEP IS TO BRING UP THOSE REFERENCE SITE PREDICTORS AND USE BOXPLOT COMPARISONS TO SEE IF OUR AIM DATA PREDICTORS FALL WITHIN THEIR RANGES.  
#     FOR EXAMPLE, THE REFERENCE DATA RANGES BETWEEN 0.700 AND 32.470 FOR THE "CaO_Mean" PREDICTOR IN THE EC MODEL, BUT FOR THE SAME PREDICTOR, THE AIM DATA RANGES BETWEEN 0.9766 AND 31.1929. SO, YOU HAVE TO DECIDE IF THE AIM PREDICTOR DATA IS SUFFICIENTLY COUCHED WITHIN THE EXPERIENCE OF THE REFERENCE SITE DATA, WHICH WOULD SUGGEST THAT THE MODEL CAN APPROPRIATELY BE RUN FOR THE AIM DATA. 
#     WE RUN 3 SEPARATE BOXPLOT COMPARISONS OF REFERENCE SITE PREDICTOR VALUES VS. AIM SITE PREDICTOR VALUES - ONE PANEL OF BOXPLOT COMPARISONS FOR EACH OF THE 3 DIFFERENT WATER QUALITY MODEL PREDICTOR SETS.
############################################################################################################################

library(randomForest)

# I need to bring UID into this somehow!!!!!!!!!!!!!!!!

# REF sites (EC, TN, TP)
directory="\\\\share1.bluezone.usu.edu\\miller\\"#to allow for multiple machines to run from share
WD=sprintf("%sbuglab\\OE_Modeling\\NAMC_Supported_OEmodels\\WQ models\\",directory)
ECref=read.csv(paste0(WD,'/Rfiles/EC_Model_ReferenceData_cleaned_JC_28April2016.csv'))
TNTP_ref=read.csv(paste0(WD,'Rfiles\\TPTN_data_cleaned_JC_28April2016.csv'))
TNTP_ref=TNTP_ref[,1:28]

# # AIM sites (EC, TN, TP)
# WDcp = "C:/Users/Christian/999_trash/PredictorRequests/20171214_Jennifer_WQ_Request/20180122_utah_shapes/"
setwd("\\\\share1.bluezone.usu.edu\\miller\\buglab\\OE_Modeling\\NAMC_Supported_OEmodels\\WQ models\\")

EC_AIM=read.csv('Z:\\buglab\\Research Projects\\AIM\\Analysis\\WQ_Bug_Model_GIS_stats_and_results\\2019\\EC_stats_from_Alex2140.csv')
rownames(EC_AIM) = EC_AIM$UID
names(EC_AIM)[names(EC_AIM) == 'SITE_ID'] = 'SiteCode'
EC_AIM$SiteCode=as.factor(EC_AIM$SiteCode)
EC_AIM$Type = "AIM"
EC_AIM = EC_AIM[,c("SiteCode", "AtmCa", "AtmMg", "AtmSO4","BDH_AVE", "CaO_Mean", "EVI_MaxAve", "KFCT_AVE", "LPREM_mean", "LST32AVE", "MAXWD_WS", "MEANP_WS", "MgO_Mean", "MINP_WS", "MINWD_WS", "PRMH_AVE", "S_Mean", "SumAve_P", "TMAX_WS","UCS_Mean", "XWD_WS","Type")]

TNTP_AIM=read.csv('Z:\\buglab\\Research Projects\\AIM\\Analysis\\WQ_Bug_Model_GIS_stats_and_results\\2020\\TNTP_NVAIM2020_Priority_Predictors.csv')
rownames(TNTP_AIM) = TNTP_AIM$UID
names(TNTP_AIM)[names(TNTP_AIM) == 'SITE_ID'] = 'SiteCode'
TNTP_AIM$SiteCode=as.factor(TNTP_AIM$SiteCode)
TNTP_AIM$Type = "AIM"
TNTP_AIM = TNTP_AIM[,c("SiteCode", "CaO_Mean", "TP_Mean", "Vol_ave", "AtmCa", "AtmNa", "AtmNO3", "AtmSO4", "TMIN_WS", "XWD_WS", "RH_WS", "PT_Tmin", "AWC_soil", "Db3rdbar", "Kfact", "Pct_Alfi", "SOC", "EVI_AveAve", "Evergr_ave", "alru_dom", "Wb_mx_area", "slpavg", "GW_P_Sp_Mx", "ER13", "PPT_ACCUM", "PPT_2MoAvg", "DOY", "Type")]


# Combine and rearrange order of REF and AIM sets:
ECboxplot=rbind(ECref,EC_AIM)
ECboxplot=ECboxplot[,c(6,13,17,21,9,2,3,4,14,12,18,19,22,11,10,5,8,16,7,20)] # In the order in which they are listed in the saved RandomForest model 
TNTPboxplot=rbind(TNTP_ref,TNTP_AIM)
refpred=refpred[,c(1,17,2,3,4,5,25,12,11,18,15,13,16,21,22,24)] # In the order in which they are listed in the saved RandomForest model 
refpred=refpred[,c(1,6,7,8,9,10,23,27,19,14,26,18,20)] # In the order in which they are listed in the saved RandomForest model 


# Boxplots: These will write as .png's directly to the current working directory
getwd()

# EC 
png('AIM2020_EC12_boxplots.png',height=1000,width=2000,units="px")
par(mfrow=c(3,7))
for (i in 1:(length(ECboxplot)-1)) {
  boxplot(ECboxplot[,i]~ECboxplot$Type, main=names(ECboxplot[i]),col=c(7,4))
}
dev.off()

# TP
png('AIM2019_TP_boxplots.png',height=1000,width=2000,units="px")
par(mfrow=c(3,5))
for (i in 1:(length(TPboxplot)-2)) {
  boxplot(TPboxplot[,i]~TPboxplot$Type, main=names(TPboxplot[i]),col=c(7,4))
}
dev.off()

#TN:
png('AIM2019_TN_boxplots.png',height=1000,width=2000,units="px")
par(mfrow=c(3,4))
for (i in 1:(length(TNboxplot)-1)) {
  boxplot(TNboxplot[,i]~TNboxplot$Type, main=names(TNboxplot[i]),col=c(7,4))
}
dev.off()


############################################################################################################################
# THIS CHUNK OF CODE RETRIEVES THE 3 EXISTING RANDOM FOREST MODELS FOR EC, TP, AND TN.  
# YOU WILL RUN THE TEST DATA (E.G., "EC_AIM") THROUGH EACH OF THE 3 MODELS RESULTING IN PREDICTIONS OF EC, TP, AND TN.
############################################################################################################################

#1. CHRISTIAN: THIS STEP IS NECESSARY SINCE THE AIM DATA DOES NOT HAVE THE "Y" LEVEL PRESENT FOR THE "ER13" VARIABLE IN ANY ROWS, BUT THE TRAINING DATA DOES:
str(TNTP_AIM); levels(TNTP_AIM$ER13)
levels(TNTP_AIM$ER13) = c(levels(TNTP_AIM$ER13),"Y")
levels(TNTP_AIM$ER13); str(TNTP_AIM)

#2. load R objects (package of formula, RF original input, RF model, etc)
load(sprintf("%sRfiles\\rf17bCnd9",WD)); rf17bCnd9; #EC
load(sprintf("%sRfiles\\rf11accTP12m2c",WD)); rf11accTP12m2c; #TP
load(sprintf("%sRfiles\\rf24wmTN12a",WD)); rf24wmTN12a; #TN:

#3. run RF in predict mode:
EC_AIM_RF_out=as.data.frame(predict(rf17bCnd9,EC_AIM))
colnames(EC_AIM_RF_out) = c("EC_predicted")
write.csv(EC_AIM_RF_out,'AIM2020_2140EC12_model_predictions.csv')

TN_AIM_RF_out=as.data.frame(predict(rf24wmTN12a,TNTP_AIM)) 
colnames(TN_AIM_RF_out) = c("TN_predicted")
write.csv(TN_AIM_RF_out,"AIM2020TN_model_predictions.csv")

TP_AIM_RF_out=as.data.frame(predict(rf11accTP12m2c,TNTP_AIM))
colnames(TP_AIM_RF_out) = c("TP_predicted")
write.csv(TP_AIM_RF_out,"AIM2020TP_model_predictions.csv")





# Now, we wish to create a single file with columns for EC_predicted, TN_predicted, and TP_predicted as well as additional columns for all the predictor variables needed to run these 3 models.

# First, we'll bring in the predictor variables:

EC_AIM_predictors=EC_AIM
rownames(EC_AIM_predictors) = EC_AIM_predictors$UID
names(EC_AIM_predictors)[names(EC_AIM_predictors) == 'SITE_ID'] = 'SiteCode'
EC_AIM_predictors$SiteCode=as.factor(EC_AIM_predictors$SiteCode)
EC_AIM_predictors$Type = "AIM"
EC_AIM_predictors = EC_AIM_predictors[,c("UID","SiteCode", "AtmCa", "AtmMg", "AtmSO4","BDH_AVE", "CaO_Mean", "EVI_MaxAve", "KFCT_AVE", "LPREM_mean", "LST32AVE", "MAXWD_WS", "MEANP_WS", "MgO_Mean", "MINP_WS", "MINWD_WS", "PRMH_AVE", "S_Mean", "SumAve_P", "TMAX_WS", "Type","UCS_Mean", "XWD_WS")]

TNTP_AIM_predictors=read.csv(paste0(WDcp,"AIM2019_sites_TNTP12_predictors.csv"))
rownames(TNTP_AIM_predictors) = TNTP_AIM_predictors$UID
names(TNTP_AIM_predictors)[names(TNTP_AIM_predictors) == 'SITE_ID'] = 'SiteCode'
TNTP_AIM_predictors$SiteCode=as.factor(TNTP_AIM_predictors$SiteCode)
TNTP_AIM_predictors$Type = "AIM"
TNTP_AIM_predictors = TNTP_AIM_predictors[,c("UID","SiteCode", "CaO_Mean", "TP_Mean", "Vol_ave", "AtmCa", "AtmNa", "AtmNO3", "AtmSO4", "TMIN_WS", "XWD_WS", "RH_WS", "PT_Tmin", "AWC_soil", "Db3rdbar", "Kfact", "Pct_Alfi", "SOC", "EVI_AveAve", "Evergr_ave", "alru_dom", "Wb_mx_area", "slpavg", "GW_P_Sp_Mx", "ER13", "PPT_ACCUM", "PPT_2MoAvg", "DOY", "Type")]

# Now, we'll merge all the necessary files together and eliminate any duplicate predictor variables (e.g., "AtmCa" is required for both the EC and TP models).
ALL_WQ_predictions_and_predictors = cbind(EC_AIM_RF_out,EC_AIM)
ALL_WQ_predictions_and_predictors = cbind(ALL_WQ_predictions_and_predictors,TNTP_AIM)
ALL_WQ_predictions_and_predictors = cbind(ALL_WQ_predictions_and_predictors,TN_AIM_RF_out)
ALL_WQ_predictions_and_predictors = cbind(ALL_WQ_predictions_and_predictors,TP_AIM_RF_out)
ALL_WQ_predictions_and_predictors = ALL_WQ_predictions_and_predictors[,c(3,2,1,54,55,4:53)]
ALL_WQ_predictions_and_predictors = ALL_WQ_predictions_and_predictors[,c(1:23,25:26,30:31,33:34,36,38:44)]
names(ALL_WQ_predictions_and_predictors)
write.csv(ALL_WQ_predictions_and_predictors,"AIM2019_ECTNTP_predictions_dec2020.csv")
