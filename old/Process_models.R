# R objects
#AREMP
load("final_models/AREMP/AREMP_OE_standardized.rdata")
load("final_models/AREMP/AREMP_MMI.rdata")
#NV
load("final_models/NV_MMI/OE_MMI_models.Rdata")
#OR MWCF
load("final_models/OR_MWCF/MWCF.rdata")
#load('final_models/OR_MWCF/Nov05model_MWCF_16jan13.Rdata')
#OR WCCP
load("final_models/OR_WCCP/WCCP.RData")
#load('Nov05model_WCCP_16jan13.RData')
#PIBO
load("final_models/PIBO/Benkendorf.RF.Model.Version1_standardized.Rdata")
#UT
load('final_models/UTDEQ15/UTDEQ_15_OE_model_standardized.rdata')
#WY
load('final_models/WYDEQ/WYRIVPACS2012_standardized.Rdata')
#WW
load("final_models/WW18/My.RF.Model.Version1_standardized.Rdata")

load("model data/TN.Rdata")
load("model data/EC12.Rdata")
load("model data/TP.Rdata")

library(NAMCr)
# get model id for what model you are running
models=NAMCr::query("models")

# get what predictors are needed for needed for that model
predictors = query("predictors",expand_metadata = FALSE)

# for all samples in a box or project

    # get all predictor values needed for a box or project # note this either needs a loop written over it or a different API end point 
    samplePredictorValues=query("samplePredictorValues",sampleId=155612) #need list of samples in database with values
    
    #subset the predictor values in the database to get only those needed for that model
    samplePredictorValues=subset(samplePredictorValues,predictorId %in% predictors$predictorId)
    
    # if predictor alkalinity or conductivity query model result table for predictor
    
    # if not null append this with other predictors
    
    # pivot the predictor values into the habitat data input
    library(tidyr)
    prednew=tidyr::pivot_wider(samplePredictorValues,names_from="abbreviation",values_from="predictorValue")# add id_cols=sampleId once it gets added to end point
    
    modelInfo=query("modelInfo",modelId=4,expand_metadata = FALSE)
    # if modelType=bug MMI or bug OE then
    
        # if modelType= bug OE get OTU taxa matrix
        bugsOTU=query("sampleMetrics",translationId=2,fixedCount=300,sampleIds=c(115217))
        rarefiedOTUTaxa=subset(bugsOTU,metricName=="Rarefied Taxa")
        rarefiedOTUTaxa = NAMCr::json.expand(rarefiedOTUTaxa, "metricValue")
        sumrarefiedOTUTaxa=rarefiedOTUTaxa  %>% dplyr::group_by(sampleId,OTUName) %>% dplyr::summarize(sumSplitCount=sum(splitCount)) # why are multiple records exported here per OTU???
        sumrarefiedOTUTaxa$presence=ifelse(sumrarefiedOTUTaxa$sumSplitCount>=1,1,0)
        bugnew=tidyr::pivot_wider(sumrarefiedOTUTaxa, id_cols="sampleId",names_from="OTUName",values_from="presence")
        
        # if modelType= bug MMI get 
        # need to add list of required metrics to model table as well as random forest input file
        #get needed list of metrics from new end point
        # NV and AREMP relative abundances were OTU standardized too!
        bugsMetrics=query("sampleMetrics",fixedCount=300,sampleIds=c(115217))
        modelInfo=query("modelInfo",modelId=27)
        MMI_metrics=subset(bugsMetrics,metricId %in% ())
        # need to replace metric name with a model specific metric abbreviation
        test_bugs_metrics=tidyr::pivot_wider(MMI_metrics, id_cols="sampleId",names_from="metricName",values_from="metricValue")
    
        # if model_id=co or CSCI
        bugsRaw=query("sampleTaxa",sampleIds=c(115217))
        
    # otherwise do nothing
    
    # create a function for each model that has habitat/predictorvalues, bug inputs, and random forest input
    
    
    
    # WQ models will have predictor and R object inputs
    
    #run model applicability function
    
    # get only predictors needed for model applicability 
    predictors = query("predictors",modelId=3,expand_metadata = FALSE)
        
    # get all predictor values needed for a box or project # note this either needs a loop written over it or a different API end point 
    samplePredictorValues=query("samplePredictorValues",modelId=) #need list of samples in database with values
        
    source(ModelApplicabilityAll.R)
        
    #save output in database as model result
    
        
    #query the model result table to get conditions automatically applied
    modelConditions=query("modelConditions",modelId=3)
    modelResults=query("modelResults", sampleIds=150807)

