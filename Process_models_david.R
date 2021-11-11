
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


####### run predictors for one sample at a time
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
    def_model_results = NAMC::query(
      api_endpoint = "ModelResults",
      include = c(
        "modelId",
        "status",
        "modelVersion",
        "abbreviation",
        #"O",
        #"E",
        "modelResult",
        #"modelApplicability"
         ),
      sampleId = sampleId,
      modelId = modelId
    )
    def_models = def_models[def_models$status = !"Valid", ]
    
    # ---------------------------------------------------------------
    # get model metadata neede to run the model
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
    #def_predictors = def_predictors[def_predictors$status = "Valid", ]
    
    if (any(def_predictors$status==!"Valid")) {
      
    }
    # if any predictors not valid look for them in model results (predicted coductivity and alalinity 
    #if not in model results write out error calculated predictor values first and stop!

    # get predictors into wide format needed for model functions
    prednew=tidyr::pivot_wider(def_predictors,names_from="abbreviation",values_from="predictorValue")# add id_cols=sampleId once it gets added to end point

    # ---------------------------------------------------------------
    # Get bug data for model functions
    # ---------------------------------------------------------------
    #
    # if modelType= bug OE get OTU taxa matrix
    if (def_models$modelType == "OE") {
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
    }
    
    # if modelType= bug MMI get
    else if (def_models$modelType == "MMI") {
      # need to add list of required metrics to model table as well as random forest input file
      #get needed list of metrics from new end point
      # NV and AREMP relative abundances were OTU standardized too!
      bugsMetrics = query("sampleMetrics",
                          def_models$fixedCount,
                          sampleIds = def_model_results$SampleId)
      modelInfo = query("modelInfo", modelId = def_model_results$modelId)
      #MMI_metrics = subset(bugsMetrics, metricId %in% ())
      # need to replace metric name with a model specific metric abbreviation
      bugnew = tidyr::pivot_wider(MMI_metrics,id_cols = "sampleId",names_from = "metricName",values_from = "metricValue")
    }
    
    else if (def_models$modelId %in% c(1, 4, 5, 6)) {
      bugsRaw = query("sampleTaxa", sampleIds = def_model_results$SampleId)
      sites =query(
        "samples", 
        include= c("sampleId",'siteId'),
        sampleIds = def_model_results$SampleId
        )
      if (def_models$modelId == 1) {
        CSCIbugs = bugsRaw
        CSCIbugs$SampleID = sites$siteId
        CSCIbugs$StationCode = #### need to get site info!!!
          CSCIbugs$FinalID = CSCIbugs$OTUName
        CSCIbugs$LifeStageCode = ifelse(
          CSCIbugs$lifeStage)
        # #LifeStageCode = case 
        # when class <> 'insecta' then 'X'
        # when class is null then 'X'
        # else Lifestage
        # end,
        CSCIbugs$Distinct == 0
        CSCIbugs$BAResult = sum(CSCIbugs$splitCount, CSCIbugs$bigRare)
        }
        else {
          CObugs= bugsRaw
          CObugs$Project
          CObugs$Station
          CObugs$Name
          CObugs$Location
          CObugs$CollDate
          CObugs$Organism
          CObugs$Individuals
          CObugs$Stage
          CObugs$CommentsTaxa
          CObugs$RepNum
          CObugs$Grids
          CObugs$CommentsSample
          CObugs$CommentsRep
          }
    }
    else {}
    
    # ---------------------------------------------------------------
    # load model specific R objects which include reference bug data and predictors RF model objects
    # ---------------------------------------------------------------   
    #AREMP
    load("final_models/AREMP/AREMP_OE_standardized.rdata")
    load("final_models/AREMP/AREMP_MMI.rdata")
    #NV
    load("final_models/NV/OE_MMI_models.rdata")
    #OR MWCF
    load("final_models/OR_MWCF/MWCF.rdata")
    #OR WCCP
    load("final_models/OR_WCCP/WCCP.RData")
    #PIBO
    load("final_models/PIBO/Benkendorf.RF.Model.Version1_standardized.Rdata")
    #UT
    load('final_models/UTDEQ15/UTDEQ_15_OE_model_standardized.rdata')
    #WY
    load('final_models/WYDEQ/WYRIVPACS2012_standardized.Rdata')
    #WW
    load("final_models/WW18/My.RF.Model.Version1_standardized.Rdata")
    # ---------------------------------------------------------------
    # Run models
    # ---------------------------------------------------------------   
    # models using latest version of van scickle function
    
    if (def_models$modelId %in% c(7,2,25,26)){
       OE<-model.predict.RanFor.4.2(bugcal.pa,grps.final,preds.final, ranfor.mod,prednew,bugnew,Pc=0.5,Cal.OOB=FALSE)
    }
    else if (def_models$modelId %in% c(10,11,13:23)){
     OE<-model.predict.v4.1(bugcal.pa,grps.final,preds.final, grpmns,covpinv,prednew,bugnew,Pc=0.5)
    } 
    else if (def_models$modelId==9){
      OE<-PIBO_model(bug.otu,bugall, grps.final, preds.final, ranfor.mod, prednew)
    } 
    else if (def_models$modelId==1){
    OE<-CSCI(bugs=CSCIbugs,stations=prednew)  
    }
    else{}
   
       
    # uses environment[[ function_name ]]() syntax to call each predictor function
    # Data needs to be in json format
  
    predictor_value = model_fns[[def_models$calculationScript]](
      polygon2process = polygon2process ,
      point2process =  geojsonsf:geojson_sf(def_sites$location[1]) ,
      predictor_name = predictor$abbreviation,
      predictor_geometry = pred_geometries[[predictor$abbreviation]],
      geometry_input_path <-
        paste0(pred_geometry_base_path, predictor$geometry_file_path),
      CurrentYear = lubridate::year(def_samples$sampleDate[1]),
      JulianDate = lubridate::yday(def_samples$sampleDate[1]),
      USGS_NED = USGS_NED
    )
    
    
    # ---------------------------------------------------------------
    # Save predictors
    # ---------------------------------------------------------------
    if (predictor$isTemporal) {
      NAMCr::save(
        api_endpoint = "newSamplePredictorValue",
        sampleId = def_samples$sampleId[1],
        predictorId = predictor$predictorId,
        value = predictor_value
      )
    } else{
      NAMCr::save(
        api_endpoint = "newSitePredictorValue",
        siteId = def_samples$siteId[1],
        predictorId = predictor$predictorId,
        value = predictor_value
      )
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
# -----------------------------------------------------------------------
# -----------------------------------------------------------------------
# Alternative ways of getting a list of predictors needed for each sample

# def_sites_models = NAMCr::query(
#   api_endpoint = "siteModels",
#   include = c("siteId", "modelId"),
#   siteId = def_samples$siteId
# )
#
# # ---------------------------------------------------------------
# # Query for the model / predictor definitions
# # ---------------------------------------------------------------
# #
# def_models = NAMCr::query(
#   api_endpoint = "models",
#   include = c("modelId","abbreviation",""),
#
# )
#
# def_models = def_models[ def_models$modelId %in% unique(def_sites_models$modelId), ]
# #**** Need missing API endpoint to join models to predictors
# # def_models_predictors = NAMCr::query(...)
#
# #modify "predictors" endpoint to accept multiple model Ids
# def_predictors = NAMCr::query(
#   api_endpoint = "predictors",
#   include = c("predictorId","abbreviation","calculationScript"),
#   modelId = unique( def_models_predictors$modelId )
# )
#
#
#
# # ---------------------------------------------------------------
# # Query for previously calculated predictor values
# # ---------------------------------------------------------------
# #
# # If state is dirty allow the script to recalculate site predictors
# if( def_samples$predictorState[1] != "dirty" ) {
#   def_sitePredictorValues = NAMCr::query(
#     api_endpoint = "sitePredictorValues",
#     include = c("predictorId"),
#     siteId = def_samples$siteId
#   )
#   # Removed already calculated predictors from list to process
#   def_predictors = def_predictors %>%
#     filter( !(predictorId %in% def_sitePredictorValues$predictorId) )
# }


