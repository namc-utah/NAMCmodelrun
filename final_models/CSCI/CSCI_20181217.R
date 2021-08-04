#install.packages("devtools")#Install devtools from CRAN
#library(devtools)
#install_github("SCCWRP/BMIMetrics")
#install_github("SCCWRP/CSCI")
#Load the CSCI library
library(CSCI)

#Upload your data

setwd("Z:\\buglab\\OE_Modeling\\NAMC_Supported_OEmodels\\CA\\Hybrid_CA_Model\\InputsAndResults\\CurrentRun")
bugs<-read.csv("bug_input.csv")
bugs.reqd<-c("StationCode","SampleID","FinalID","BAResult","LifeStageCode","Distinct")
bugs.reqd %in% names(bugs) #Are all your required fields in "bugs"?
bugs.reqd[!bugs.reqd %in% names(bugs)] #Which predictors are missing?


stations<-read.csv("habitat.csv")
stations$LogWSA<-log10(stations$AREA_SQKM)

predictors.reqd<-c("StationCode","AREA_SQKM","New_Long","New_Lat","SITE_ELEV","ELEV_RANGE",
                   "TEMP_00_09","PPT_00_09","SumAve_P","KFCT_AVE","BDH_AVE","MgO_Mean","P_MEAN",
                   "CaO_Mean","PRMH_AVE","S_Mean","PCT_SEDIM","LPREM_mean","N_MEAN")


predictors.reqd %in% names(stations) #Are all required predictors in "stations?
predictors.reqd[!predictors.reqd %in% names(stations)] # List of missing fields


#stations<-stations[,predictors.reqd] 

#Create the CSCI function
#This will be incorporated into the final CSCI package

report<-CSCI(bugs,stations)
#SWJ 4/5/2013: will get an error if there are mismatches in final id and lifestage that prevents the function from running. Correct these and rerun, report the mismatches to Sarah Judson and Raphael Mazor so they can be updated in the buglab database and script taxa list respectively.
ls(report) #see all the components of the report.

Core<-report$core #Core results for each sample
Supp1.mmi<-report$Suppl1_mmi #mean metrics for each sample
Supp1.grps<-report$Suppl1_grps #group probabilities for each site
Supp1.OE<-report$Suppl1_OE #Capture probabilities and mean abundance of each OTU
Supp2.mmi<-report$Suppl2_mmi #Metric for each iteration of each sample, plus the means again for some reason
Supp2.OE<-report$Suppl2_OE #Abundance of each OTU for each iteration for each sample

head(report$core)

write.csv(Core,"core.csv", row.names=F)
write.csv(Supp1.mmi,"Supp1.mmi.csv", row.names=F)
write.csv(Supp1.grps,"Supp1.grps.csv", row.names=F)
write.csv(Supp1.OE,"Supp1.OE.csv", row.names=F)
write.csv(Supp2.mmi,"Supp2.mmi.csv", row.names=F)
write.csv(Supp2.OE,"Supp2.OE.csv", row.names=F)


