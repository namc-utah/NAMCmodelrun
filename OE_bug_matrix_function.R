
translationId=def_models$translationId
fixedCount
sampleId

bugsOTU = query(
  "sampleMetrics",
  translationId = def_models$translationId,
  fixedCount = def_models$fixedCount,
  sampleIds = def_model_results$SampleId
)
rarefiedOTUTaxa = subset(bugsOTU, metricName == "Rarefied Taxa")
rarefiedOTUTaxa = NAMCr::json.expand(rarefiedOTUTaxa, "metricValue")
sumrarefiedOTUTaxa = rarefiedOTUTaxa  %>%
  dplyr::group_by(sampleId, OTUName) %>%
  dplyr::summarize(sumSplitCount = sum(splitCount)) # why are multiple records exported here per OTU???

sumrarefiedOTUTaxa$presence = ifelse(sumrarefiedOTUTaxa$sumSplitCount >=1, 1, 0)
bugnew = tidyr::pivot_wider(sumrarefiedOTUTaxa,id_cols = "sampleId", names_from = "OTUName",values_from = "presence")
