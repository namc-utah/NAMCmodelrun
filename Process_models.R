# get model id for what model you are running
models=query("models")

# get what predictors are needed for needed for that model
predictors = query("predictors",modelId=3,expand_metadata = FALSE)

# get all predictor values needed for a box or project # note this either needs a loop written over it or a different API end point 
samplePredictorValues=query("samplePredictorValues",modelId=155612) #need list of samples in database with values

#subset the predictor values in the database to get only those needed for that model
samplePredictorValues=subset(samplePredictorValues,predictorId %in% predictors$predictorId)

# pivot the predictor values into the habitat data input
library(tidyr)
predictorvalues=tidyr::pivot_wider(samplePredictorValues,names_from="abbreviation",values_from="predictorValue")# add id_cols=sampleId once it gets added to end point

modelInfo=query("modelInfo",modelId=4,expand_metadata = FALSE)
# if modelType=bug MMI or bug OE then

    # get the bug data needed for the model
    sampleMetrics=query("sampleMetrics",translationId=2,fixedCount=300,sampleIds=c(115217))
    
    # if modelType= bug OE get OTU taxa matrix
    rarefiedOTUTaxa=subset(sampleMetrics,metricName=="Rarefied Taxa")
    rarefiedOTUTaxa = NAMCr::json.expand(rarefiedOTUTaxa, "metricValue")
    sumrarefiedOTUTaxa=rarefiedOTUTaxa  %>% dplyr::group_by(sampleId,OTUName) %>% dplyr::summarize(sumSplitCount=sum(splitCount)) # why are multiple records exported here per OTU???
    sumrarefiedOTUTaxa$presence=ifelse(sumrarefiedOTUTaxa$sumSplitCount>=1,1,0)
    bugInput=tidyr::pivot_wider(sumrarefiedOTUTaxa, id_cols="sampleId",names_from="OTUName",values_from="presence")
    
    # if modelType= bug MMI get 
    # need to add list of required metrics to model table as well as random forest input file
    #get needed list of metrics from new end point
    modelInfo=query("modelInfo",modelId=27)
    MMI_metrics=subset(sampleMetrics,metricId %in% ())
    # need to replace metric name with a model specific metric abbreviation
    bugInput=tidyr::pivot_wider(MMI_metrics, id_cols="sampleId",names_from="metricName",values_from="metricValue")

# otherwise do nothing

# create a function for each model that has habitat/predictorvalues, bug inputs, and random forest input



# WQ models will have predictor and R object inputs

#run model applicability function

#save output in database as model result

#query the model result table to get conditions automatically applied
modelConditions=query("modelConditions",modelId=3)
modelResults=query("modelResults", sampleIds=150807)

