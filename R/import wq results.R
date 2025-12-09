finalResults=read.csv("C:/Users/A01964518/Box/NAMC/OEModeling/NAMC_Supported_OEmodels/WaterQuality/Updated WQ Models/import into database.csv")
for (i in 1:nrow(finalResults)){

  tryCatch({
 NAMCr::save(
  api_endpoint = "setModelResult",
  sampleId = finalResults$sampleId[i],
  modelId = finalResults$modelId[i],
  oResult=finalResults$o_result[i],
  eResult=finalResults$e_result[i],
  modelResult = finalResults$modelResult[i],
  modelApplicability = finalResults$ModelApplicability[i],
  fixedCount=0
)
  }, error =function(e){
    cat(paste0("\n\tSAMPLE ERROR: results may already exist in database and were not overwritten",finalResults$sampleId[i],"\n"))
    str(e,indent.str = "   "); cat("\n")
  })

}
