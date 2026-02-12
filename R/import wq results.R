WQ=read.csv("C:/Users/A01964518/Box/NAMC/Research Projects/AIM/WQ data/WQ_Bug_Model_GIS_stats_and_results/2025_predictions/2025_predictedWQimport.csv")
t=tidyr::pivot_longer(WQ,cols=c(EC_predict,	EC_ModelApplicability,	TN_predict,	TN_ModelApplicability	,TP_predict	,TP_ModelApplicability,	TN_Observed	,TP_Observed,	EC_Observed
))
library(tidyverse)
df_separated_delim <- t %>%
  separate_wider_delim(
    cols = name, # Column to split
    delim = "_",     # Delimiter character (a space)
    names = c("Model", "param") # Names for the new columns
  )
write.csv(df_separated_delim,'wq.csv')
WQ=read.csv('wq.csv')
finalResults=tidyr::pivot_wider(WQ,names_from=param,values_from=value)
finalResults=na.omit(finalResults)
finalResults$modelResult=finalResults$oResult-finalResults$eResult

finalResults=read.csv("C:/Users/A01964518/Box/NAMC/OEModeling/NAMC_Supported_OEmodels/WaterQuality/Updated WQ Models/import into database.csv")
for (i in 1:nrow(finalResults)){

  tryCatch({
 NAMCr::save(
  api_endpoint = "setModelResult",
  sampleId = finalResults$sampleId[i],
  modelId = finalResults$modelId[i],
  oResult=finalResults$oResult[i],
  eResult=finalResults$eResult[i],
  modelResult = finalResults$modelResult[i],
  modelApplicability = finalResults$ModelApplicability[i],
  fixedCount=0
)
  }, error =function(e){
    cat(paste0("\n\tSAMPLE ERROR: results may already exist in database and were not overwritten",finalResults$sampleId[i],"\n"))
    str(e,indent.str = "   "); cat("\n")
  })

}
