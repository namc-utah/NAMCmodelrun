# Code for getting Bugs and WQ from NAMC for BLM AIM lotic indicator table or for national reporting #

# Step 1 get an account with NAMC by signing up here

#https://namc-usu.org/    then to SampleProcessing/Online Sample Submission (INSTAR) to create an account.
# Contact Trip or Jennifer to get your role set to customer and external collaborator


# Step 2 install package to connect to NAMC database - only need to do once
remotes::install_github("namc-utah/NAMCr", force=TRUE)

#test your credentials by running the code below. If you receive an error contact Trip or Jennifer
NAMCr::query("auth")



#### Step 3 determine what results exist for all AIM sites ####
# all AIM sites = projectId=49 but may take a bit to load.
# each field season/ state has its own boxId
# results may not exist but bugs make have been processed. To view if bugs have been processed review list of AIM boxes by running code below
boxes=NAMCr::query("boxes",include = c('boxId','alias','boxState','sampleCount','boxReceivedDate','processingCompleteDate','projectedCompleteDate','sampleCount'),
                   entityIds=614)
#all AIM data
samples=NAMCr::query("samples",projectId=49)
#Only data for a certain box or boxes of interest
#samples=NAMCr::query("samples",boxId=c()#input list of boxes of interest from above boxes query
Report=NAMCr::query("modelResults", sampleIds=samples$sampleId)
Report=subset(Report, is.na(modelResult)==FALSE)


#### Step 4 format output into format needed for AIM database ####

# bug data
# subset and order result columns
Report2=Report[,c("sampleId",'visitId','customerSiteCode','modelId','modelAbbr','fixedCount','oResult','eResult','modelApplicability','modelResult','condition','notes')]
# rename columns for the AIM database
Report2=setNames(Report2,c('SampleID','EvaluationID','PointID','modelId','OE_MMI_ModelUsed','MacroinvertebrateCount','ObservedInvertRichness','ExpectedInvertRichness','OE_MMI_ModelApplicability','modelResult','condition','InvasiveInvertSpecies'))
Report2$OE_MMI_ModelApplicability=ifelse(Report2$OE_MMI_ModelApplicability=='TRUE','Pass',ifelse(Report2$OE_MMI_ModelApplicability=='FALSE','Fail',Report2$OE_MMI_ModelApplicability))
Report2$MMI_Macroinvertebrate=ifelse(Report2$modelId %in% c(3,4,5,6,8,24),Report2$modelResult,NA)
Report2$OE_Macroinvertebrate=ifelse(Report2$modelId %in% c(1,2,7,9,10,11,12,13:23,25:26),Report2$modelResult, NA)
# get standard field office level results for AIM database
Report3=subset(Report2, modelId %in% c(1:7,9:26) & InvasiveInvertSpecies!='National')
###If you want results for a national report get westwide OE scores comment out above line and uncomment this line... invasive species should be gotten only from field office level query above
##Report3=subset(Report2, modelId %in% c(25,26))


# WQ data (predicted WQ values only) raw WQ data back from the Baker lab is not stored by NAMC and rather just excel spreadsheet is forwarded on
WQ=subset(Report2,modelId %in% c(27,28,29), select=c('SampleID','EvaluationID','PointID','OE_MMI_ModelUsed','modelResult'))
WQ$OE_MMI_ModelUsed=ifelse(WQ$OE_MMI_ModelUsed=='TN','PredictedTotalNitrogen',
                           ifelse(WQ$OE_MMI_ModelUsed=='TP','PredictedTotalPhosphorous',
                                  ifelse(WQ$OE_MMI_ModelUsed=='EC12','PredictedSpecificConductance',
                                         WQ$OE_MMI_ModelUsed)))
WQ_pivot=tidyr::pivot_wider(WQ,id_cols=c('SampleID','EvaluationID','PointID'),names_from=c('OE_MMI_ModelUsed'), values_from=c("modelResult"))

#join bug and WQ data
Final=dplyr::full_join(Report3,WQ_pivot, by=c('SampleID','EvaluationID','PointID'))
write.csv(Final,paste0('AIM_Bug_WQ_results_for_import_',Sys.Date(),'.csv'))




#future modifications if table structure changes
#notescol=setNames(as.data.frame(stringr::str_split_fixed(Report$notes,";",2)),c("InvasiveInverSpecies","Purpose"))
#Report=cbind(Report,notescol)
#Report=subset(Report,Pupose=="StateBasedAIM")
