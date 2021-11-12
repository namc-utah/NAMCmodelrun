translationId

bugsTranslation = query("sampleTaxaTranslation", translationId=def_models$translationId,sampleIds = def_model_results$SampleId)
sites =query("samples", 
             include= c("sampleId",'siteId'),
             sampleIds = def_model_results$SampleId
)
bugRaw=query("sampleRawTaxa", include=c(), translationId=def_models$translationId,sampleIds = def_model_results$SampleId)# raw query with pivoted taxonomy, and join translation name but not roll it up.... then summ in here
taxonomy=query("taxonomy hearchy")
if (def_models$modelId == 1) {# can be 
  CSCIbugs = bugsRaw
  CSCIbugs$SampleID = 
    CSCIbugs$StationCode = sites$siteId#### need to get site info!!!
  CSCIbugs$FinalID = CSCIbugs$OTUName
  CSCIbugs$LifeStageCode = ifelse(
    CSCIbugs$lifeStage)
  # #LifeStageCode = case 
  # when class <> 'insecta' then 'X'
  # when class is null then 'X'
  # else Lifestage
  # end,
  CSCIbugs$Distinct == 0
  CSCIbugs$BAResult = sum(CSCIbugs$splitCount, CSCIbugs$bigRare)# have query add this as a column
  