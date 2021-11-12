
###### function to run all samples in the database at once. API endpoint still needs developed
#' run all predictors for all samples at once
#'
#' @return
#' @export
#'
#' @examples
process_models = function() {
  logger = Logger$new(
    logPath = "./",
    fileName = "",
    enabled = TRUE,
    consoleOutput = TRUE,
    appendDate = TRUE
  )
  logger$startLog()


  def_samples = NAMCr::query(api_endpoint = "samples2process",#need separate end point for models
                             include = c("sampleId", "modelId"))
  for (sample in def_samples) {
    process_sample_models(sample$sampleId, sample$modelId) 
  }
  logger$stopLog()
}


####### run predictors for a box
#' predictors for each box
#'
#' @param boxId
#'
#' @return
#' @export
#'
#' @examples
process_box_models = function(boxId,modelId) {
  logger = Logger$new(
    logPath = "/",
    fileName = "",
    enabled = TRUE,
    consoleOutput = TRUE,
    appendDate = TRUE
  )
  logger$startLog()

  def_boxes = NAMCr::query(
    api_endpoint = "samples",
    include = c("sampleId"),
    boxIds = boxId
  )

  # for (i in seq_len(nrow(def_boxes))) {
  #   process_sample_predictors(def_boxes$sampleId[i])
  # }

  by(def_boxes, seqlen(nrow(def_boxes)), function(sample) {
    process_sample_models(sample$sampleId,modelId)
  })

  logger$stopLog()
}


####### run models for one sample at a time
#' Process sample models
#' @description
#' @details saving each predictor for each sample one at a time in the database
#'
#' @param sampleId
#' @param config
#'
#' @return none
#' @export
#'
#' @examples
process_sample_models = function(sampleId, modelId, config = config) {
  tryCatch({
    # ---------------------------------------------------------------
    # determine if the needed model has already been run for the sample
    # ---------------------------------------------------------------
   # consider adding a loop through models here 
    
    # getting a list of samples and associated models that do not already have results in the table
    # has site location changed or predictor values changed if not dont rerun , handle with valid flag
    # not ready status if predictors are not there
    def_model_results = NAMC::query(
      api_endpoint = "modelResults",
      include = c(
        "modelId",
        "status",
        "modelVersion",
        "abbreviation",
        #"O",
        #"E",
        "modelResult",
        #"modelApplicability"
        #invasive species
        #is.Active ---look for only one record
        #modelRunId - query knows if metrics are out of date
         ),
      sampleId = sampleId,
      modelId = modelId
    )
    def_models = def_models[def_models$status = !"Valid", ]
    # if any predictors not valid look for them in model results (predicted coductivity and alalinity
    
    # ---------------------------------------------------------------
    # get model metadata needed to run the model
    # ---------------------------------------------------------------
    def_models=NAMC::query(
      api_endpoint = "ModelInfo",
      include = c(
        "modelId",
        "modelType",
       "translationId",
        "fixed count",
        ),
      modelId = def_model_results$modelId
    )
    
    # ---------------------------------------------------------------
    # Get predictor values needed for the model, if they dont exist yet stop here
    # ---------------------------------------------------------------
           
    # getting predictor values associated with those samples and models coming out of the def_models query above
    def_predictors = NAMCr::query(
      api_endpoint = "samplePredictorValues",
      include = c("predictorId",
                  "status",
                  "abbreviation",
                  "predictorValue"),
      sampleId = def_models$SampleId,
      modelId = def_models$modelId
    )
    if (any(def_predictors$status==!"Valid")) {
      
    }
     
    #if not in model results write out error calculated predictor values first and stop!

    # get predictors into wide format needed for model functions
    prednew=tidyr::pivot_wider(def_predictors,names_from="abbreviation",values_from="predictorValue")# add id_cols=sampleId once it gets added to end point

    # ---------------------------------------------------------------
    # Get bug data for model functions
    # ---------------------------------------------------------------
    #
    # if modelType= bug OE get OTU taxa matrix
    if (def_models$modelType == "OE") {
      bugnew=OE_bug_matrix(sampleId=sampleId,translationId=def_models$translationId,fixedCount=def_models$fixedCount)
        }
    # if modelType= bug MMI get
    else if (def_models$modelType == "MMI") {
      bugnew=MMI_metrics(sampleId=sampleId,fixedCount=def_models$fixedCount)
     }
    # CSCI requires just the raw taxa list translated for misspelling
    else if (def_models$modelId %in% c(1)) {# add in model names in comments
      CSCIbugs=CSCI_bug(sampleId=sampleId, translationId=def_models$translationId)
      }
    # CO model must be written out as an excel file using a separate bank of code and function
    else if (def_models$modelId %in% c(4,5,6)){
      print(paste0("CO model must be written out as an excel file using a separate bank of code and function: ",sampleIds))
    }
    else {}
    
    # ---------------------------------------------------------------
    # load model specific R objects which include reference bug data and predictors RF model objects
    # ---------------------------------------------------------------   
    # every model has an R object that stores the random forest model and reference data
    # the R objects are named with the model abbreviation
    load(paste0(def_models$model_abbreviation,".Rdata"))
        
    # ---------------------------------------------------------------
    # Run models
    # ---------------------------------------------------------------   
    # models using latest version of van sickle function include:
    if (def_models$modelId %in% c(7,2,25,26)){
       OE<-model.predict.RanFor.4.2(bugcal.pa,grps.final,preds.final, ranfor.mod,prednew,bugnew,Pc=0.5,Cal.OOB=FALSE)#....
    }
    # models using version 4.1 of van sickle code includes: 
    else if (def_models$modelId %in% c(10,11,13:23)){
     OE<-model.predict.v4.1(bugcal.pa,grps.final,preds.final, grpmns=grpmns,covpinv,prednew,bugnew,Pc=0.5)# add elpsis...
    } 
    # PIBO model was one of the earliest models built and used a version of van sickle's function that didnt have a version number
    else if (def_models$modelId==9){
      OE<-PIBO_model(bug.otu,bugall, grps.final, preds.final, ranfor.mod, prednew)
    } 
    # CSCI has its own package and function    
    else if (def_models$modelId==1){
    report<-CSCI(bugs=CSCIbugs,stations=prednew) 
    OE=report$core
    }
    # all MMIs will need their own function added here because there is a rf model for each metric
    else if (def_models$modelId==8){
      AREMP_MMI_model(bugnew,prednew,rf_model) #seperate out all the models here 
    }
    # all MMIs will need their own function added here because there is a rf model for each metric
    else if (def_models$modelId==3){
     NV_MMI_model(bugnew,prednew,rf_model)  #seperate out all the models here 
    }
     
    else{}
   
       
    )
    
   
    # ---------------------------------------------------------------
    # Save model results
    # ---------------------------------------------------------------
   #has permission to save then spit out result to console
     # pass Nas for anything not used
      NAMCr::save(
        api_endpoint = "newModelResult",
        sampleId = def_samples$sampleId[1],
        modelId = def_model_results=modelId,
        O= ifelse (not NA),
        E=,
        model_result= ,
        modelApplicability=      
          )
    } 
      ) 
      }
   
    }
  },
  error = function(e) {
    cat(paste0("\n\tPREDICTOR ERROR: ", predictor$abbreviation, "\n"))
    str(e, indent.str = "   ")
    cat("\n")
  })
})
},
error = function(e) {
  cat(paste0("\n\tSAMPLE ERROR: ", sampleId, "\n"))
  str(e, indent.str = "   ")
  cat("\n")
})
}


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


