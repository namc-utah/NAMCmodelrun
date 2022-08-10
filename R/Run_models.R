#'
#' ####### run models for a list of samples in a box
#' #' Run models
#' #' @description
#' #' @details run models for a list of samples (i.e. box) and save each model result sample by sample
#' #'
#' #' @param boxIds
#' #' @param modelIds
#' #'
#' #' @return none
#' #' @export
#' #'
#' #' @examples
#' run_models = function(boxIds, modelIds ) {
  # tryCatch({

#   # ---------------------------------------------------------------
# get a list of samples in a box or project
# ---------------------------------------------------------------
if (exists("boxId")){
  def_samples=NAMCr::query("samples",boxId=boxId)
}else {def_samples=NAMCr::query("samples",projectId=projectId)
  }

  #   # ---------------------------------------------------------------
  # get a list of samples if the needed model has already been run for the sample
  # ---------------------------------------------------------------

  # getting a list of samples and associated models that do not already have results in the table
  # has site location changed or predictor values changed if not dont rerun by excluding status=current
  # not ready status if predictors are not there
  def_model_results = NAMCr::query(
    api_endpoint = "modelResults",
    sampleIds=def_samples$sampleId
  )

  def_model_results=subset(def_model_results,modelId==modelID)


  # ---------------------------------------------------------------
  # get model metadata needed to run the model - philip said he would change apis so that model id would just be provided and translation id and fixed count wouldnt be needed
  # ---------------------------------------------------------------
  def_models = NAMCr::query(
    api_endpoint = "modelInfo",
    include = c("modelId",
                "modelTypeAbbreviation",
                "abbreviation",
                "translationId",
                "fixedCount"),
    modelId = def_model_results$modelId[1]
  )

  # ---------------------------------------------------------------
  # Get predictor values needed for the model, if they dont exist yet stop here
  # ---------------------------------------------------------------

  # getting predictor values associated with those samples and models coming out of the def_models query above
  def_predictors = NAMCr::query(
    api_endpoint = "samplePredictorValues",
    include = c("sampleId",
                "predictorId",
                "status",
                "abbreviation",
                "predictorValue"
                ),
    sampleIds = def_model_results$sampleId)
  modelpred=NAMCr::query("predictors",modelId=modelID)

  if (modelID %in% c(4,5,6,28)){ #CO and TP models have predictors that are categorical but all other models need predictors converted from character to numeric after pulling from database
    def_predictors_categorical=subset(def_predictors,predictorId %in% c(111,75))
    prednew1 = tidyr::pivot_wider(def_predictors_categorical,
                                 id_cols="sampleId",
                                 names_from = "abbreviation",
                                 values_from = "predictorValue")# add id_cols=sampleId once it gets added to end point
    prednew1=as.data.frame(prednew1)
     '%notin%' <- Negate('%in%')
    def_predictors=subset(def_predictors,predictorId %notin% c(111,75))
    def_predictors$predictorValue=as.numeric(def_predictors$predictorValue)
    prednew2 = tidyr::pivot_wider(def_predictors,
                                 id_cols="sampleId",
                                 names_from = "abbreviation",
                                 values_from = "predictorValue")# add id_cols=sampleId once it gets added to end point

    prednew2=as.data.frame(prednew2)
    prednew=rbind(prednew1,pred2)

    }else {def_predictors=subset(def_predictors,predictorId %in% modelpred$predictorId)
            def_predictors$predictorValue=as.numeric(def_predictors$predictorValue)
            # if (any(def_predictors$status!="Valid")) {
            #   print(paste0("predictors need calculated"))
            # } else{
              # get predictors into wide format needed for model functions
              prednew = tidyr::pivot_wider(def_predictors,
                                           id_cols="sampleId",
                                           names_from = "abbreviation",
                                           values_from = "predictorValue")# add id_cols=sampleId once it gets added to end point
              prednew=as.data.frame(prednew)
            }
    rownames(prednew)<-prednew$sampleId
    prednew<-prednew[,-1]
    # ---------------------------------------------------------------
    # Get bug data for model functions
    # ---------------------------------------------------------------
    #
    # if modelType= bug OE get OTU taxa matrix
    if (def_models$modelId == 12) {
      bugnew = OR_NBR_bug(
        sampleIds = def_model_results$sampleId,
        translationId = def_models$translationId,
        fixedCount = def_models$fixedCount
      )
      #need a way to distinguish this model from others.. call NULL OE?
    } else if (def_models$modelTypeAbbreviation == "OE") {
      bugnew = OE_bug_matrix(
        sampleIds = def_model_results$sampleId,
        translationId = def_models$translationId,
        fixedCount = def_models$fixedCount
      )
     # CO model must be written out as an excel file using a separate bank of code and function
    } else if (def_models$modelId %in% c(4, 5, 6)) {
      #write get bugs from database, write out as a csv and save as CObugs object
      CObugs=CO_bug_export(sampleIds=def_model_results$sampleId)
      #write out predictors as a csv
      write.csv(prednew,file = paste0("COpredictors","boxId_",CObugs$Project[1],"_",Sys.Date(),".csv"),row.names=FALSE)
      cat(paste("csv with COpredictors has been written out to your current working directory.",
                "Convert this csv to excel 2003 and import into CO EDAS access database to compute the CSCI score.",
                "Follow instructions in this pdf Box\\NAMC\\OE_Modeling\\NAMC_Supported_OEmodels\\CO\\Documentation\\EDAS2017\\Tutorial Guide to EDAS_Version 1.7.pdf",
                "to import bug and habitat data, harmonize taxa list, rarefy and compute MMI",
                "then read resulting excel file back into R to save results in the database.", sep="\n"))

      # CSCI requires just the raw taxa list translated for misspelling
    } else if (def_models$modelId %in% c(1)) {
      bugnew = CSCI_bug(sampleIds = def_model_results$sampleId)
    } else if (def_models$modelId %in% c(169)){
      bugnew=AZ_bug_export(sampleIds = def_model_results$sampleId)
    }else if (def_models$modelTypeAbbreviation == "MMI") {# if modelType= bug MMI get
      bugnew = MMI_metrics(sampleIds = def_model_results$sampleId, translationId=def_models$translationId, fixedCount = def_models$fixedCount,modelId=def_models$modelId)
 }else {

    }
    if (def_models$modelId %in% c(12)){#no predictors are needed for OR null model

    }else {bugnew<-subset(bugnew,rownames(bugnew) %in% rownames(prednew))
    prednew<-subset(prednew,rownames(prednew) %in% rownames(bugnew))
    #reorder them the same just in case model functions dont already do this
    bugnew = bugnew[order(rownames(bugnew)),];
    prednew = prednew[order(rownames(prednew)),];
}


    # ---------------------------------------------------------------
    # load model specific R objects which include reference bug data and predictors RF model objects
    # ---------------------------------------------------------------
    # every model has an R object that stores the random forest model and reference data
    # the R objects are named with the model abbreviation
    # instead of all these if statements the R file name could be stored in the database... but WY and NV require two models and R file names
    #if CO, CSCI, or OR null model no R data file needs loaded in
    if (def_models$modelId %in% c(1,4,5,6,12,169)){
      print("no R object needs loaded")
      #if WY model only one Rdata file needs loaded and not one for each "model" but Alkalinity also needs added
    } else if (def_models$modelId %in% c(13:23)){
      load("sysdata.rda/WY2018.Rdata")
      load("sysdata.rda/Alkalinity.Rdata")### objects named the same so they will be overwritten.... how do we deal with
      #if westwide model only one R data file needs loaded in and not one for each model
    }else if (def_models$modelId %in% c(25:26)){
      load(paste0("sysdata.rda/Westwide2018.Rdata"))
      # all other models should have R data files named identical to model name
    }else{
      load(paste0("sysdata.rda/",def_models$abbreviation, ".Rdata"))
    }

    # ---------------------------------------------------------------
    # Run models
    # ---------------------------------------------------------------

    # ------------------------------
    # OE models
    # ------------------------------
    # models using john vansickles RIVPACS random forest code : AREMP, UTDEQ15, Westwide, PIBO
    if (def_models$modelId %in% c(7, 2, 25, 26,9)) {
      OE <-model.predict.RanFor.4.2(
          bugcal.pa,
          grps.final,
          preds.final,
          ranfor.mod,
          prednew,
          bugnew,
          Pc = 0.5,
          Cal.OOB = FALSE
        )#....
      modelResults<-OE$OE.scores

    # models using John VanSickles RIVPACS discriminant function code: OR_WCCP, OR_MWCF
    } else if (def_models$modelId %in% c(10, 11)) {
      OE <-model.predict.v4.1(bugcal.pa,
                           grps.final,
                           preds.final,
                           grpmns,
                           covpinv,
                           prednew,
                           bugnew,
                           Pc = 0.5)# add elpsis...
      modelResults<-OE$OE.scores

    # WY also uses John vansickles discriminant function code but requires alkalinity model as a dependency
    } else if (def_models$modelId %in% (13:23)) {
      ALK_LOG = setNames(as.data.frame(
        log10(predict(ranfor.mod, prednew, type = "response"))
      ), c("ALK_LOG"))
      prednew = cbind(prednew, ALK_LOG) ## need to subset to only model predictors or maybe it doesnt matter??
      OE <-model.predict.v4.1(bugcal.pa,
                           grps.final,
                           preds.final,
                           grpmns,
                           covpinv,
                           prednew,
                           bugnew,
                           Pc = 0.5)
      modelResults<-OE$OE.scores

    # OR eastern region is a null model and no predictors are used
    } else if (def_models$modelId == 12) {
      modelResults <- OR_NBR_model(bugnew)

    # CSCI has its own package and function
    } else if (def_models$modelId == 1) {
      report <- CSCI::CSCI(bugs = bugnew, stations = prednew)
      modelResults = report$core
      rownames(modelResults)=modelResults$SampleID

      # ------------------------------
      # MMI models
      # ------------------------------
      # all MMIs will need their own function added here because there is a rf model for each metric
    } else if (def_models$modelId == 8) {# AREMP MMI
      modelResults <-
        AREMP_MMI_model(
          bugnew,
          prednew,
          CLING_rich.rf,
          DIPT_rich.rf,
          LLT_rich.rf,
          NON_INSECT_rich.rf,
          PER_EPT.rf,
          PER_INTOL.rf,
          rf_models,
          mdeg_metrics_adj_cal,
          ref_metrics_adj
        )
    } else if (def_models$modelId == 3) {# NV MMI
    # need to call conductivity model first before calling the NV model because predicted conductivity is a predictor for the NV model
      load(file="sysdata.rda/EC12.Rdata")
      PrdCond = setNames(as.data.frame(
        predict(ranfor.mod, prednew, type = "response")# need to subset to only model predictors or maybe it doesnt matter??
      ), c('PrdCond'))
      prednew = cbind(prednew, PrdCond)
      modelResults <-
        NV_MMI_model(
          bugnew,
          prednew,
          CLINGER.rf,
          INSET.rf,
          NONSET.rf,
          PER_CFA.rf,
          PER_EPHEA.rf,
          PER_PLECA.rf
        )
    } else if (def_models$modelId == 136){#modeled insect richness
      UniqueRichness_Insecta_pred=predict(rfmod_UniqueRichness_Insecta, prednew, type="response")
      # join predicted insect richness to intial data
      modelResults=cbind(bugnew, UniqueRichness_Insecta_pred)
      # convert to O/E ratio
      modelResults$modeledInsectRichness=modelResults$UniqueRichness_Insecta/modelResults$UniqueRichness_Insecta_pred

      # ------------------------------
      # WQ models
      # ------------------------------
      #conductivity, tp, tn,temperature
    }else if (def_models$modelId %in% c(27, 28, 29, 30)) {
      modelResults = as.data.frame(predict(ranfor.mod, prednew, type = "response"))# make sure prednew has sampleIds as the rows
    }else{

    }

    # ---------------------------------------------------------------
    # Always run model applicability test
    # ---------------------------------------------------------------
    # get all predictor values needed for a box or project # note this either needs a loop written over it or a different API end point
    applicabilitypreds = NAMCr::query("samplePredictorValues",
                                include = c(
                                  "sampleId",
                                  "predictorId",
                                  "status",
                                  "abbreviation",
                                  "predictorValue"
                                  ),
                                sampleIds = def_model_results$sampleId
                                ) #need list of samples in database with values
    applicabilitypreds = subset(applicabilitypreds, abbreviation %in% c('ElevCat','Tmean8110Ws','WsAreaSqKm','Precip8110Ws'))
    applicabilitypreds$predictorValue=as.numeric(applicabilitypreds$predictorValue)
    applicabilitypreds = tidyr::pivot_wider(applicabilitypreds,
                                 id_cols="sampleId",
                                 names_from = "abbreviation",
                                 values_from = "predictorValue")# add id_cols=sampleId once it gets added to end point
    applicabilitypreds=as.data.frame(applicabilitypreds)
    # run model applicability function
    ModelApplicability = ModelApplicability(CalPredsModelApplicability,
                                            modelId = modelID,
                                            applicabilitypreds) # add to config file or add an R object with calpreds

    finalResults=merge(modelResults,ModelApplicability,by="row.names")

    # ---------------------------------------------------------------
    # Get additional bug metrics (fixed count and invasives)
    # ---------------------------------------------------------------
    ##### get fixed count column #####
    bugsOTU = NAMCr::query("sampleTaxaTranslationRarefied",
                           translationId = def_models$translationId,
                           fixedCount = def_models$fixedCount,
                           sampleIds=def_model_results$sampleId
    )
    sumrarefiedOTUTaxa = bugsOTU  %>%
      dplyr::group_by(sampleId) %>%
      dplyr::summarize(fixedCount = sum(splitCount))


    ###### get invasives #####
    # get raw bug data
    bugRaw = NAMCr::query(
      "sampleTaxa",
      sampleIds=def_model_results$sampleId
    )
    #subset taxa in samples to only invasives
    bugraw = subset(bugRaw,taxonomyId %in% c(1330,1331,2633, 2671,4933,4934,4935,4936,4937,4938,4939,4940,4941,4942,1019,1994,5096,1515,1518,1604,2000,4074,1369,2013,1579))
    #create list of invasives present at a site
    invasives<-bugraw %>% dplyr::group_by(sampleId) %>% dplyr::summarize(InvasiveInvertSpecies=paste0(list(unique(scientificName)),collapse=''))
    # remove list formatting
    invasives$InvasiveInvertSpecies=gsub("^c()","",invasives$InvasiveInvertSpecies)
    invasives$InvasiveInvertSpecies=gsub("\"","",invasives$InvasiveInvertSpecies)
    invasives$InvasiveInvertSpecies=gsub("\\(","",invasives$InvasiveInvertSpecies)
    invasives$InvasiveInvertSpecies=gsub("\\)","",invasives$InvasiveInvertSpecies)
    # join to list of all samples with fixed counts
    additionalbugmetrics=dplyr::left_join(sumrarefiedOTUTaxa,invasives, by="sampleId")
    # if no invasives were present set to absent
    additionalbugmetrics[is.na(additionalbugmetrics)]<-"Absent"

    finalResults=dplyr::left_join(finalResults,additionalbugmetrics,by="sampleId")


    # ---------------------------------------------------------------
    # Save model results
    # ---------------------------------------------------------------
if (overwrite=='N'){
  def_model_results2 = NAMCr::query(
    api_endpoint = "modelResults",
    sampleIds=def_samples$sampleId
  )

  def_model_results2=subset(def_model_results2,modelId==modelID & is.na(modelResult)==TRUE)
  finalResults=subset(finalResults,sampleId %in% def_model_results2$sampleId)
} else{}

for (i in 1:nrow(finalResults) ){# need to add invasives and extra metrics to the notes field in some easy fashion???
    #has permission to save then spit out result to console
    # pass Nas for anything not used
  tryCatch({
    if (def_models$modelTypeAbbreviation == "OE") {
      NAMCr::save(
        api_endpoint = "setModelResult",
        sampleId = finalResults$sampleId[i],
        modelId = def_model_results$modelId[1],
        oResult = finalResults$O[i],
        eResult = finalResults$E[i],
        modelResult = finalResults$OoverE[i] ,
        fixedCount = finalResults$fixedCount[i],
        modelApplicability = finalResults$ModelApplicability[i],
        notes=finalResults$InvasiveInvertSpecies[i]
      )
    }else if (def_models$modelTypeAbbreviation == "Both Observed Versus Expected and Multimetric Index") {
      NAMCr::save(
        api_endpoint = "setModelResult",
        sampleId = finalResults$sampleId[i],
        modelId = def_model_results$modelId[1],
        modelResult = finalResults$CSCI[i],
        fixedCount = finalResults$fixedCount[i],
        modelApplicability = finalResults$ModelApplicability[i],
        notes=finalResults$InvasiveInvertSpecies[i]
      )
    }else if (def_models$modelTypeAbbreviation == "MMI") {
      NAMCr::save(
        api_endpoint = "setModelResult",
        sampleId = finalResults$sampleId[i],
        modelId = def_model_results$modelId[1],
        modelResult = finalResults$MMI[i],
        fixedCount = finalResults$fixedCount[i],
        modelApplicability = finalResults$ModelApplicability[i],
        notes=finalResults$InvasiveInvertSpecies[i]
      )
    }else if (def_models$modelTypeAbbreviation == "WQ") {
      NAMCr::save(
        api_endpoint = "setModelResult",
        sampleId = finalResults$sampleId[i],
        modelId = def_model_results$modelId[1],
        modelResult = finalResults$WQ[i],###need to fix....
        modelApplicability = finalResults$ModelApplicability[i]
      )
    }else{
    }

  }, error =function(e){
    cat(paste0("\n\tSAMPLE ERROR: results may already exist in database and were not overwritten",finalResults$sampleId[i],"\n"))
    str(e,indent.str = "   "); cat("\n")
  })
}

#  }
#}




# query the model result table to get conditions automatically applied
Report=query("modelResults", sampleIds=def_model_results$sampleId)
Report=subset(Report,modelId==modelID)

#use the following to see what thresholds were applied
modelConditions=NAMCr::query("modelConditions",modelId=modelID)
