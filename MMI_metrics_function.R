fixCount
sampleId
# need to add list of required metrics to model table as well as random forest input file
#get needed list of metrics from new end point
# NV and AREMP relative abundances were OTU standardized too!
bugsMetrics = query("sampleMetrics",
                    def_models$fixedCount,
                    sampleIds = def_model_results$SampleId)
modelInfo = query("modelInfo", modelId = def_model_results$modelId)
#MMI_metrics = subset(bugsMetrics, metricId %in% ()) # store translations to metric names in database
# need to replace metric name with a model specific metric abbreviation
bugnew = tidyr::pivot_wider(MMI_metrics,id_cols = "sampleId",names_from = "metricName",values_from = "metricValue")
