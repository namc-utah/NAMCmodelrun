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

  def_model_results=subset(def_model_results,modelId %in% modelID)


  # ---------------------------------------------------------------
  # get model metadata needed to run the model - philip said he would change apis so that model id would just be provided and translation id and fixed count wouldnt be needed
  # ---------------------------------------------------------------
  arbitrary_list<-list()

  if(length(modelID)>1){
    print("There are multiple modelIDs in this set")
  for(i in 1:length(modelID)){
  x = NAMCr::query(
    api_endpoint = "modelInfo",
    include = c("modelId",
                "modelTypeAbbreviation",
                "abbreviation",
                "translationId",
                "fixedCount"),
    modelId = unique(def_model_results$modelId)[i]
  )
  arbitrary_list[[i]]<-x
  }
  def_models<-bind_rows(arbitrary_list)
  }else{
  print("There is only one modelID for this set")
    def_models = NAMCr::query(
      api_endpoint = "modelInfo",
      include = c("modelId",
                  "modelTypeAbbreviation",
                  "abbreviation",
                  "translationId",
                  "fixedCount"),
      modelId =def_model_results$modelId[1]
    )
    def_models=as.data.frame(do.call('rbind',def_models))
  }


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

  if(length(modelID)>1){
    for(i in 1:length(modelID)){
        modelpred=NAMCr::query("predictors",modelId=modelID[i])
        arbitrary_list[[i]]<-modelpred
    }
    modelpred<-do.call('rbind',arbitrary_list)
  }else{
    modelpred=NAMCr::query("predictors",modelId=modelID[i])
  }



  def_predictors=subset(def_predictors,predictorId %in% modelpred$predictorId)
  if (length(modelID[modelID%in%c(4,5,6,28)]==T)>=1){ #CO and TP models have predictors that are categorical but all other models need predictors converted from character to numeric after pulling from database
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
    prednew=cbind(prednew1,prednew2)

    }else {def_predictors$predictorValue=as.numeric(def_predictors$predictorValue)
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
    if (length(def_models$modelId[def_models$modelId %in% 12]==T)>=1) {
      bugnew = OR_NBR_bug(
        sampleIds = def_model_results$sampleId,
        translationId = def_models$translationId[1],
        fixedCount = def_models$fixedCount[1]
      )
      #need a way to distinguish this model from others.. call NULL OE?
    } else if (length(def_models$modelTypeAbbreviation[def_models$modelTypeAbbreviation=='OE'])>=1) {
      bugnew = OE_bug_matrix(
        sampleIds = def_model_results$sampleId,
        translationId = def_models$translationId[1],
        fixedCount = def_models$fixedCount[1]
      )
     # CO model must be written out as an excel file using a separate bank of code and function
    } else if (length(def_models$modelId[def_models$modelId %in% c(4,5,6)]==T)>=1) {
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
    } else if (length(def_models$modelId[def_models$modelId %in% 1]==T)>=1) {
      bugnew = CSCI_bug(sampleIds = def_model_results$sampleId)
    } else if (length(def_models$modelId[def_models$modelId %in% 169]==T)>=1){
      bugnew=AZ_bug_export(sampleIds = def_model_results$sampleId)
    }else if (def_models$modelTypeAbbreviation == "MMI") {# if modelType= bug MMI get
      bugnew = MMI_metrics(sampleIds = def_model_results$sampleId, translationId=def_models$translationId, fixedCount = def_models$fixedCount,modelId=def_models$modelId)
 }else {

    }
    if (length(def_models$modelId[def_models$modelId %in% 12]==T)>=1){#no predictors are needed for OR null model

    } else if (length(def_models$modelId[def_models$modelId %in% c(1,4:6)]==T)>=1){#CSCI bug file doesnt have row names and has multiple rows for a given sampleID, so does CO but CO is written out to disk
      bugnew<-subset(bugnew,bugnew$SampleID %in% rownames(prednew))
      prednew<-subset(prednew,rownames(prednew) %in% row.names(bugnew))
      #reorder them the same just in case model functions dont already do this
      bugnew = bugnew[order(row.names(bugnew)),];
      prednew = prednew[order(rownames(prednew)),];
    } else {bugnew<-subset(bugnew,rownames(bugnew) %in% rownames(prednew))
    prednew<-subset(prednew,rownames(prednew) %in% rownames(bugnew))
    #reorder them the same just in case model functions dont already do this
    bugnew = bugnew[order(rownames(bugnew)),];
    prednew = prednew[order(rownames(prednew)),];
    }

    ## special handling of AREMP predictor names is needed here.
    #Tmax_WS is used for NV WQ and other models but for AREMP this same predictor is called TMAX_WS
    #database is not case sensitive so cant add as unique predictor in database so have to handle here.
if (length(def_models$modelId[def_models$modelId %in% c(7,8)]==T)>=1){
  prednew$TMAX_WS=prednew$Tmax_WS
} else{}
# same issue with RH_WS being included in TP model but different capitalization. However, revisit if we revise TP model
if (length(def_models$modelId[def_models$modelId %in% 8]==T)>=1){
      prednew$rh_WS=prednew$RH_WS
    } else{}



    # ---------------------------------------------------------------
    # load model specific R objects which include reference bug data and predictors RF model objects
    # ---------------------------------------------------------------
    # every model has an R object that stores the random forest model and reference data
    # the R objects are named with the model abbreviation
    # instead of all these if statements the R file name could be stored in the database... but WY and NV require two models and R file names
    #if CO, CSCI, or OR null model no R data file needs loaded in
    if (length(def_models$modelId[def_models$modelId %in% c(1,4,5,6,12,169)]==T)>=1){
      print("no R object needs loaded")
      #if WY model only one Rdata file needs loaded and not one for each "model" but Alkalinity also needs added
    } else if (length(def_models$modelId[def_models$modelId %in% 13:23]==T)>=1){
      load("sysdata.rda//WY2018.Rdata")
      load("sysdata.rda//Alkalinity.Rdata")### objects named the same so they will be overwritten.... how do we deal with
      #if westwide model only one R data file needs loaded in and not one for each model
    }else if (length(def_models$modelId[def_models$modelId %in% 25:26]==T)>=1){
      load(paste0("sysdata.rda//Westwide2018.Rdata"))
      # all other models should have R data files named identical to model name
    }else{
      load(paste0("sysdata.rda//",def_models$abbreviation, ".Rdata"))
    }

    # ---------------------------------------------------------------
    # Run models
    # ---------------------------------------------------------------

    # ------------------------------
    # OE models
    # ------------------------------
    # models using john vansickles RIVPACS random forest code : AREMP, UTDEQ15, Westwide, PIBO
    if (length(def_models$modelId[def_models$modelId %in% c(2,7,25,26,29)]==T)>=1) {
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
    } else if (length(def_models$modelId[def_models$modelId %in% 10:11]==T)>=1) {
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
    } else if (length(def_models$modelId[def_models$modelId %in% 13:23]==T)>=1) {
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
    } else if (length(def_models$modelId[def_models$modelId %in% 12]==T)>=1) {
      modelResults <- OR_NBR_model(bugnew)

    # CSCI has its own package and function
    } else if (length(def_models$modelId[def_models$modelId %in% 1]==T)>=1) {
      prednew$sampleId=as.numeric(row.names(prednew))
      prednew=left_join(prednew,def_samples[,c('siteName','sampleId')],by='sampleId')
      prednew$StationCode<-prednew$sampleId
      report <- CSCI::CSCI(bugs = bugnew, stations = prednew)
      modelResults = report$core
      rownames(modelResults)=modelResults$SampleID

      # ------------------------------
      # MMI models
      # ------------------------------
      # all MMIs will need their own function added here because there is a rf model for each metric
    } else if (length(def_models$modelId[def_models$modelId %in% 8]==T)>=1) {# AREMP MMI
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
    } else if (length(def_models$modelId[def_models$modelId %in% 3]==T)>=1) {# NV MMI
    # need to call conductivity model first before calling the NV model because predicted conductivity is a predictor for the NV model
      load(file="sysdata.rda//EC12.Rdata")
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
    } else if (length(def_models$modelId[def_models$modelId %in% 136]==T)>=1){#modeled insect richness
      UniqueRichness_Insecta_pred=predict(rfmod_UniqueRichness_Insecta, prednew, type="response")
      # join predicted insect richness to intial data
      modelResults=cbind(bugnew, UniqueRichness_Insecta_pred)
      # convert to O//E ratio
      modelResults$modeledInsectRichness=modelResults$UniqueRichness_Insecta/modelResults$UniqueRichness_Insecta_pred

      # ------------------------------
      # WQ models
      # ------------------------------
      #conductivity, tp, tn,temperature
    }else if (length(def_models$modelId[def_models$modelId %in% 27:30]==T)>=1) {
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
    if(length(modelID)>1){
      for(i in 1:length(modelID)){
    ModelApplicability_obj = ModelApplicability(CalPredsModelApplicability,
                                            modelId = modelID[i],
                                            applicabilitypreds) # add to config file or add an R object with calpreds

    arbitrary_list[[i]]<-ModelApplicability_obj
      }
      ModelApplicability_obj<-do.call('rbind',arbitrary_list)
      finalResults=merge(modelResults,ModelApplicability_obj,by="row.names")
      }else{
      ModelApplicability_obj = ModelApplicability(CalPredsModelApplicability,
                                              modelId = modelID,
                                              applicabilitypreds) # add to config file or add an R object with calpreds

      finalResults=merge(modelResults,ModelApplicability_obj,by="row.names")
}
    # ---------------------------------------------------------------
    # Get additional bug metrics (fixed count and invasives)
    # ---------------------------------------------------------------
    ##### get fixed count column #####
    bugsOTU = NAMCr::query("sampleTaxaTranslationRarefied",
                           translationId = def_models$translationId[1],
                           fixedCount = def_models$fixedCount[1],
                           sampleIds=def_model_results$sampleId
    )
    sumrarefiedOTUTaxa = bugsOTU  %>%
      dplyr::group_by(sampleId) %>%
      dplyr::summarize(fixedCount = sum(splitCount))

    ################################################
    ###### get invasives ##### comment out this entire invasives section if running "National" AIM westwide reporting
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
    #################################################

    #IF NATIONAL COMMENT OUT THIS LINE OF CODE AND UNCOMMENT OUT THE FOLLOWING TWO LINES
    finalResults=dplyr::left_join(finalResults,additionalbugmetrics,by="sampleId")
    # finalResults=dplyr::left_join(finalResults,sumrarefiedOTUTaxa,by="sampleId")
    # finalResults$InvasiveInvertSpecies='National'


    # ---------------------------------------------------------------
    # Create modelId column to save results by appropriate modelId and write results out to csv
    # ---------------------------------------------------------------
    # join in sample info and coordinates
    finalResults=dplyr::left_join(finalResults,def_samples[,c('sampleId','siteLongitude','siteLatitude')],by='sampleId')
    finalResults_sf=sf::st_as_sf(finalResults,coords=c('siteLongitude','siteLatitude'),crs=4269)
    # create modelId column specific to geography
    if (length(modelID[modelID %in% 25:26]==T)>=1){
    ecoregion=sf::st_read(paste0(ecoregion_base_path,"GIS//GIS_Stats//CONUS//ecoregion//hybrid_level_III_ecoregions.shp"))
    ecoregion=sf::st_make_valid(ecoregion)
    finalResults_sf=sf::st_transform(finalResults_sf,5070)
    finalResults_sf=sf::st_intersection(finalResults_sf,ecoregion)
    finalResults=dplyr::left_join(finalResults,sf::st_drop_geometry(finalResults_sf[,c('sampleId','modelId')]),by='sampleId')
    } else if (length(modelID[modelID %in% 13:23]==T)>=1){
    ecoregion=sf::st_read(paste0(ecoregion_base_path,"GIS//GIS_Stats//Wyoming//ecoregion//BIOREGIONS_2011_modifiedCP.shp"))
    ecoregion=sf::st_make_valid(ecoregion)
    finalResults_sf=sf::st_transform(finalResults_sf,5070)
    finalResults_sf=sf::st_intersection(finalResults_sf,ecoregion)
    finalResults=dplyr::left_join(finalResults,sf::st_drop_geometry(finalResults_sf[,c('sampleId','modelId')]),by='sampleId')
    } else{
    finalResults$modelId=modelID
    }

    write.csv(finalResults,paste0('finalresults_',paste(modelID,collapse='_'),"_",Sys.Date(),'.csv'))

    # ---------------------------------------------------------------
    # Save model results
    # ---------------------------------------------------------------
if (overwrite=='N'){
  def_model_results2 = NAMCr::query(
    api_endpoint = "modelResults",
    sampleIds=def_samples$sampleId
  )
  def_model_results2=subset(def_model_results2,modelId%in%modelID & is.na(modelResult)==TRUE)
  #def_model_results2=subset(def_model_results2,modelId==modelID & is.na(modelResult)==TRUE)
  finalResults=subset(finalResults,sampleId %in% def_model_results2$sampleId)
  names(finalResults)[1]<-'sampleId' #changing the name from row.names to sampleId
  finalResults$sampleId<-as.integer(finalResults$sampleId) #making it a good class, not "AsIs"
} #else{}


for (i in 1:nrow(finalResults) ){# need to add invasives and extra metrics to the notes field in some easy fashion???
    #has permission to save then spit out result to console
    # pass Nas for anything not used
  tryCatch({
    if (length(def_models$modelTypeAbbreviation[def_models$modelTypeAbbreviation == "OE"]>=1)) {
      print('O/E')
      dat_to_pass<-list(sampleId = finalResults$sampleId[i],
                        modelId = finalResults$modelId[i],
                        oResult = finalResults$O[i],
                        eResult = finalResults$E[i],
                        modelResult = finalResults$OoverE[i] ,
                        fixedCount = finalResults$fixedCount[i],
                        modelApplicability = finalResults$ModelApplicability[i],
                        notes=finalResults$InvasiveInvertSpecies[i])
      NAMCr::save(
        api_endpoint = "setModelResult",
        args=dat_to_pass)
    }else if (length(def_models$modelTypeAbbreviation[def_models$modelTypeAbbreviation == "Hybrid"]>=1)) {
      print('Hybrid')
      NAMCr::save(
        api_endpoint = "setModelResult",
        sampleId = finalResults$sampleId[i],
        modelId = finalResults$modelId[i],
        modelResult = finalResults$CSCI[i],
        fixedCount = finalResults$fixedCount[i],
        modelApplicability = finalResults$ModelApplicability[i],
        notes=finalResults$InvasiveInvertSpecies[i]
      )
    }else if (length(def_models$modelTypeAbbreviation[def_models$modelTypeAbbreviation == "MMI"]>=1)) {
      print('MMI')
      NAMCr::save(
        api_endpoint = "setModelResult",
        sampleId = finalResults$sampleId[i],
        modelId = finalResults$modelId[i],
        modelResult = finalResults$MMI[i],
        fixedCount = finalResults$fixedCount[i],
        modelApplicability = finalResults$ModelApplicability[i],
        notes=finalResults$InvasiveInvertSpecies[i]
      )
    }else if (length(def_models$modelTypeAbbreviation[def_models$modelTypeAbbreviation == "WQ"]>=1)) {
      print('WQ')
      NAMCr::save(
        api_endpoint = "setModelResult",
        sampleId = finalResults$sampleId[i],
        modelId = finalResults$modelId[i],
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
Report=subset(Report,modelId%in%modelID)

#use the following to see what thresholds were applied
if(length(modelID)>1){
  for(i in 1:length(modelID)){
    modelConditions=NAMCr::query("modelConditions",modelId=modelID[i])
    arbitrary_list[[i]]<-modelConditions
  }
  modelConditions<-bind_rows(arbitrary_list)
}else{
modelConditions=NAMCr::query("modelConditions",modelId=modelID)
}
write_model<-paste(modelID,collapse='_')
write.csv(Report,paste('Report_',write_model,"_",Sys.Date(),'.csv'))

