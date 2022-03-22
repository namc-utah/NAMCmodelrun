

Report=NAMCr::query("modelResults", projectId=49)
Report=subset(Report, is.na(modelResult)==FALSE)
notescol=setNames(as.data.frame(stringr::str_split_fixed(t$notes,";",2)),c("InvasiveInverSpecies","Purpose"))
Report=cbind(Report,notescol)
Report=subset(Report,Pupose=="StateBasedAIM")

Report2=Report[,c("sampleId",'visitId','customerSiteCode','modelId','modelAbbr','fixCount','oResult','eResult','modelApplicability','modelResult','condition','notes')]

Report2=setNames(Report2,c('SampleID','EvaluationID','PointID','modelId','OE_MMI_ModelUsed','MacroinvertebrateCount','ObservedInvertRichness','ExpectedInvertRichness','OE_MMI_ModelApplicability','modelResult','condition','InvasiveInvertSpecies'))
Report3=subset(Report2, modelId %in% c(1:7,9:24))
#test=as.data.frame(tidyr::pivot_wider(Report2,id_cols=c('SampleID','EvaluationID','PointID'),names_from=c('OE_MMI_ModelUsed'), values_from=c("modelResult",'ObservedInvertRichness')))
oddbugs=subset(Report2,modelId %in% c(8,25,26))
join=dplyr::full_join(Report3,oddbugs, by=c('SampleID','EvaluationID','PointID'))
# first merge na rows then join in WQ after getting rid of x

# in the notes field put "NationalReporting" and filter those samples out or put "AIMStatebased" and filter only those samples. that way any research uses of data are filtered out
# if doing this then invasives needs made a metric, put in the notes as InvasiveInvertSpecies: or made its own model (preferable)

WQ=subset(Report2,modelId %in% c(27,28,29), select=c('SampleID','EvaluationID','PointID','OE_MMI_ModelUsed','modelResult'))
WQ_pivot=tidyr::pivot_wider(WQ,id_cols=c('SampleID','EvaluationID','PointID'),names_from=c('OE_MMI_ModelUsed'), values_from=c("modelResult"))

Final=full_join(Report3,WQ_pivot, by=c('SampleID','EvaluationID','PointID'))


Report2$MMI_Macroinvertebrate=ifelse(Report2$modelId %in% c(3,4,5,6,8,24),Report2$modelResult,NA)
Report2$OE_Macroinvertebrate=ifelse(Report2$modelId %in% c(1,2,7,9,10,11,12,13:23,25:26),Report2$modelResult, NA)
Report2$PredictedTotalNitrogen=ifelse(modelId=27,Report2$modelResult,NA)
Report2$PredictedTotalPhosphorous=ifelse(modelId=28,Report2$modelResult,NA)
Report2$PredictedSpecificConductance=ifelse(modelId=29,Report2$modelResult,NA)

