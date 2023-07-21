#Mutliple modelID modelrun
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

#OE_list equates to bugnew for OE indices not null

#model types
#1 - OE
#2 - MMI
#3 - CSCI
#4- WQ
basic_models = query(
  api_endpoint = "models"
)

OEs<-basic_models$modelId[basic_models$modelTypeId==1 & basic_models$modelId !=12]
MMIs<-basic_models$modelId[basic_models$modelTypeId==2]
CSCIs<-basic_models$modelId[basic_models$modelTypeId==3]
WQs<-basic_models$modelId[basic_models$modelTypeId==4]
AREMP<-c(7,8)
NullOE<-12
MIR<-136 #change
PIBO <-9 #change


if (exists("boxId")){
  def_samples=NAMCr::query("samples",boxId=boxId)
}else {def_samples=NAMCr::query("samples",projectId=projectId)
}

sampleIds = def_samples$sampleId
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

#subset only necessary modelIDs. for example,
#box 2162 has 4 models assigned to it, but only 3 were needed
#based on geography
modelID<-modelID[modelID %in% def_model_results$modelId]


sampleOEs<-def_model_results[def_model_results$modelId %in% OEs,]
sampleMMIs<-def_model_results[def_model_results$modelId %in% MMIs,]
sampleCSCIs<-def_model_results[def_model_results$modelId %in% CSCIs,]
sampleWQs<-def_model_results[def_model_results$modelId %in% WQs,]
sampleNullOEs<-def_model_results[def_model_results$modelId %in% NullOE,]
sampleAREMP<-def_model_results[def_model_results$modelId %in% AREMP,]
sampleMIRs<-def_model_results[def_model_results$modelId %in% MIR,]
samplePIBO<-def_model_results[def_model_results$modelId %in% PIBO,]
# ---------------------------------------------------------------
# get model metadata needed to run the model - philip said he would change apis so that model id would just be provided and translation id and fixed count wouldnt be needed
# ---------------------------------------------------------------
arbitrary_list<-list()

if(length(unique(def_model_results$modelId))>1){
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
  def_models=as.data.frame(t(as.data.frame(do.call('rbind',def_models))))
  def_models[,c(1,4,5)]<-as.integer(def_models[,c(1,4,5)])
}

# ---------------------------------------------------------------
# load model specific R objects which include reference bug data and predictors RF model objects
# ---------------------------------------------------------------
# every model has an R object that stores the random forest model and reference data
# the R objects are named with the model abbreviation
# instead of all these if statements the R file name could be stored in the database... but WY and NV require two models and R file names
#if CO, CSCI, or OR null model no R data file needs loaded in

#for loop that loops through the model Ids and brings in only necessary support.
for(i in 1:nrow(def_models)){
  print(paste('iteration:',i,sep=' '))
  if (def_models$modelId[i] %in%  c(1,4,5,6,12,169)) {
    print(modelID[modelID %in% c(1,4,5,6,12,169)])
    print("no R object needs loading")
    #if WY model only one Rdata file needs loaded and not one for each "model" but Alkalinity also needs added
  } else if (def_models$modelId[i] %in% 13:23){
    print(modelID[modelID %in% 13:23])
    load("sysdata.rda//WY2018.Rdata")
    load("sysdata.rda//Alkalinity.Rdata")### objects named the same so they will be overwritten.... how do we deal with
    #if westwide model only one R data file needs loaded in and not one for each model
  } else if (def_models$modelId[i] %in% 25:26){
    print(modelID[modelID %in% 25:26])
    load(paste0("sysdata.rda//Westwide2018.Rdata"))
    # all other models should have R data files named identical to model name
  }else{
    load(paste0("sysdata.rda//",def_models$abbreviation[i], ".Rdata"))
    print('reading in support fxn that was not listed')
  }

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


def_predictors <- def_predictors[!duplicated(def_predictors), ]

  if(length(modelID)>1){
    for(i in 1:length(modelID)){
      print('there are multiple modelIDs')
      modelpred=NAMCr::query("predictors",modelId=modelID[i])
      arbitrary_list[[i]]<-modelpred
    }
    modelpred<-do.call('rbind',arbitrary_list)
  }else{
    print("there is only 1 modelID here")
    modelpred=NAMCr::query("predictors",modelId=modelID)
  }


def_predictors=subset(def_predictors,predictorId %in% modelpred$predictorId)
if (nrow(sampleMMIs)>=1 | nrow(sampleWQs)>=1){ #CO and TP models have predictors that are categorical but all other models need predictors converted from character to numeric after pulling from database
  print("CO or TP model")
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
  rownames(prednew)<-prednew$sampleId
  prednew<-prednew[,-1]

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
rownames(prednew)<-prednew$sampleId
prednew<-prednew[,-1]
}



#get bug info for the models

if(nrow(sampleNullOEs)>=1){
  Nullbugnew = OR_NBR_bug(
    sampleIds = sampleNullOEs$sampleId,
    translationId = sampleNullOEs$translationId[1],
    fixedCount = sampleNullOEs$fixedCount[1]
  )
  #run the model right here!
  NullOR_modelResults <- OR_NBR_model(Nullbugnew)
}

if(nrow(sampleOEs)>=1){
  print('Running O/E')
  OE_list<-list()
  for(i in 1:length(unique(sampleOEs$modelId))){
    y<-sampleOEs[sampleOEs$modelId==unique(sampleOEs$modelId)[i],]
  bugnew = OE_bug_matrix(
    sampleIds =y$sampleId,
    translationId = y$translationId[1],
    fixedCount = y$fixedCount[1])
  OE_list[[i]]<-bugnew
  names(OE_list)[i]<-unique(sampleOEs$modelId)[i]
  }

  if (nrow(sampleOEs[sampleOEs$modelId %in% c(2,7,9,25,26,29),])>=1) {
    print('running O/Es using 4.2RF')
    RunOE_list<-list()
    bug_sub_list<-sampleOEs$modelId[sampleOEs$modelId %in% c(2,7,9,25,26,29),]
    n_unique_OE_mods<-length(unique(bug_sub_list$modelId))
    for(i in 1:n_unique_OE_mods){
      model_id_burn<-unique(bug_sub_list$modelId)[1]
      oe_bug_burn<-as.data.frame(OE_list[as.character(model_id_burn)])
      colnames(oe_bug_burn)<-sub(c('X2.','X7.','X9.','X25.','X26.','X29.'),colnames(oe_bug_burn))
      OE <-model.predict.RanFor.4.2(
        bugcal.pa,
        grps.final,
        preds.final,
        ranfor.mod,
        prednew,
        bugnew=oe_bug_burn,
        Pc = 0.5,
        Cal.OOB = FALSE)
      modelResults<-OE$OE.scores
      RunOE_list[[i]]<-modelResults
      names(RunOE_list)[i]<-model_id_burn
    }
    #if only one model met that condition above
    #force the list to a df
    if(length(OE_list)==1)
      OEmodelResults<-as.data.frame(RunOE_list[[1]])



  }
}else if (sampleMMIs$modelId %in% c(4,5,6)) {
    print('CO MMI, this is run in Access. Only predictors and bugs needed')
    #write get bugs from database, write out as a csv and save as CObugs object
    CObugs=CO_bug_export(sampleIds=def_model_results$sampleId)
    #write out predictors as a csv
    write.csv(prednew,file = paste0("COpredictors","boxId_",CObugs$Project[1],"_",Sys.Date(),".csv"),row.names=FALSE)
    cat(paste("csv with COpredictors has been written out to your current working directory.",
              "Convert this csv to excel 2003 and import into CO EDAS access database to compute the CSCI score.",
              "Follow instructions in this pdf Box\\NAMC\\OE_Modeling\\NAMC_Supported_OEmodels\\CO\\Documentation\\EDAS2017\\Tutorial Guide to EDAS_Version 1.7.pdf",
              "to import bug and habitat data, harmonize taxa list, rarefy and compute MMI",
              "then read resulting excel file back into R to save results in the database.", sep="\n"))
} else if (nrow(sampleCSCIs)>=1) {
  bugnew = CSCI_bug(sampleIds = def_model_results$sampleId)
} else if (sampleMMIs$modelId %in% 169){
  bugnew=AZ_bug_export(sampleIds = def_model_results$sampleId)
}else if (nrow(sampleMMIs)>=1) {
  MMI_list<-list()
  for(i in 1:length(unique(sampleMMIs$modelId))){
    m<-sampleMMIs[sampleMMIs$modelId==unique(sampleMMIs$modelId)[i],]
  bugnew = MMI_metrics(sampleIds = m$sampleId,
                       translationId=m$translationId[1],
                       fixedCount = m$fixedCount[1],
                       modelId=m$modelId[1])
  }
}
def_predictors=subset(def_predictors,predictorId %in% modelpred$predictorId)

if (nrow(sampleWQs)>=1){ #TP models have predictors that are categorical but all other models need predictors converted from character to numeric after pulling from database
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
  rownames(prednew)<-prednew$sampleId
  prednew<-prednew[,-1]

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
rownames(prednew)<-prednew$sampleId
prednew<-prednew[,-1]
}

if (nrow(sampleCSCIs)>=1 | sampleMMIs$modelId %in% 4:6){#CSCI bug file doesnt have row names and has multiple rows for a given sampleID, so does CO but CO is written out to disk
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

if (nrow(sampleAREMP)>=1){
  prednew$TMAX_WS=prednew$Tmax_WS
}
# same issue with RH_WS being included in TP model but different capitalization. However, revisit if we revise TP model
if (nrow(sampleAREMP[sampleAREMP$modeId==8])>=1){
  prednew$rh_WS=prednew$RH_WS
}



# ---------------------------------------------------------------
# Run models
# ---------------------------------------------------------------

# ------------------------------
# OE models
# ------------------------------
# models using john vansickles RIVPACS random forest code : AREMP, UTDEQ15, Westwide, PIBO


if (nrow(sampleOEs[sampleOEs$modelId %in% c(2,7,9,25,26,29),])>=1) {
  print('running O/Es using 4.2RF')
  RunOE_list<-list()
  bug_sub_list<-sampleOEs$modelId[sampleOEs$modelId %in% c(2,7,9,25,26,29),]
  n_unique_OE_mods<-length(unique(bug_sub_list))
  for(i in 1:n_unique_OE_mods){
    model_id_burn<-unique(bug_sub_list$modelId)[i]
    oe_bug_burn<-OE_list[as.character(model_id_burn)]
    OE <-model.predict.RanFor.4.2(
      bugcal.pa,
      grps.final,
      preds.final,
      ranfor.mod,
      prednew,
      bugnew=oe_bug_burn,
      Pc = 0.5,
      Cal.OOB = FALSE)
    modelResults<-OE$OE.scores
    RunOE_list[[i]]<-modelResults
    names(RunOE_list)[i]<-model_id_burn
  }
  #if only one model met that condition above
  #force the list to a df
  if(length(OE_list)==1)
    OEmodelResults<-as.data.frame(RunOE_list[[1]])


} if (nrow(sampleOEs[sampleOEs$modelId %in% 10:11,])>=1) {
  print('running O/E using 4.1RF')
  bug_sub_list<-sampleOEs[sampleOEs$modelId %in% 10:11,]
  RunOE_L<-list()
  for(i in 1:length(bug_sub_list)){
    model_id_burn<-unique(bug_sub_list$modelId)[i]
    oe_bug_burn<-OE_list[as.character(model_id_burn)]
    OE <-model.predict.v4.1(bugcal.pa,
                            grps.final,
                            preds.final,
                            grpmns,
                            covpinv,
                            prednew,
                            bugnew=oe_bug_burn,
                            Pc = 0.5)
    modelResults<-OE$OE.scores
    RunOE_L[[i]]<-modelResults
  }
  if(length(OE_L)==1)
    OE_L<-as.data.frame(OE_L)

} else if (nrow(sampleOEs[sampleOEs$modelId %in% 13:23,])>=1) {
    print ('WY models')
    WY_sub_list<-sampleOEs$modelId[sampleOEs$modelId %in% 13:23]
    n_unique_WYmods<-length(unique(WY_sub_list))
    ALK_LOG = setNames(as.data.frame(
      log10(predict(ranfor.mod, prednew, type = "response"))
    ), c("ALK_LOG"))
    prednew = cbind(prednew, ALK_LOG)
    WY_list<-list()
    ## need to subset to only model predictors or maybe it doesnt matter??
    for(i in 1:length(n_unique_WYmods)){
      model_id_burn<-unique(WY_sub_list$modelId)[i]
      oe_bug_burn<-OE_list[as.character(model_id_burn)]
    OE <-model.predict.v4.1(bugcal.pa,
                            grps.final,
                            preds.final,
                            grpmns,
                            covpinv,
                            prednew,
                            bugnew=oe_bug_burn,
                            Pc = 0.5)
    modelResults<-OE$OE.scores
    WY_list[[i]]<-modelResults
    names(WY_list)[i]<-model_id_burn

    }
  #prep data for CSCI, if needed
} else if (nrow(sampleCSCIs)>=1) {
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
 } else if (nrow(sampleAREMP[sampleAREMP$modelId==8])>=1) {# AREMP MMI
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

   }else if (nrow(sampleMMIs[sampleMMIs$modelId==3,])>=1) {# NV MMI
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

} else if (nrow(sampleMIRs)>=1){#modeled insect richness
     UniqueRichness_Insecta_pred=predict(rfmod_UniqueRichness_Insecta, prednew, type="response")
     # join predicted insect richness to intial data
     modelResults=cbind(bugnew, UniqueRichness_Insecta_pred)
     # convert to O//E ratio
     modelResults$modeledInsectRichness=modelResults$UniqueRichness_Insecta/modelResults$UniqueRichness_Insecta_pred

     # ------------------------------
     # WQ models
     # ------------------------------
     #conductivity, tp, tn,temperature
}else if (nrow(sampleWQs)>=1) {
  modelResults = as.data.frame(predict(ranfor.mod, prednew, type = "response"))# make sure prednew has sampleIds as the rows
}

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


#It looks like model 12 breaks the applicabilit test!
#is this something that we know of, or a weird bug?
#removing it from the model app makes it work just fine.
modelID_for_app<-modelID[modelID!=12]
model_app_list<-list()
    for(i in 1:length(modelID)){
      ModelApplicability_obj = ModelApplicability(CalPredsModelApplicability,
                                                  modelId = modelID_for_app[i],
                                                  applicabilitypreds) # add to config file or add an R object with calpreds

      model_app_list[[i]]<-ModelApplicability_obj
      names(model_app_list)[i]<-modelID_for_app[i]
    }
model_app_list[[1]]
OE_final_results<-list()
for(i in 1:length(modelID_for_app)){
  model_app_burn<-model_app_list[as.character(modelID_for_app[i])]
  OE_final_results[[i]]<-plyr::join(OE_L[[i]],model_app_burn,'row.names','left')
}
  #arbitrary_list[[i]]<-ModelApplicability_obj
  ModelApplicability_obj<-do.call('rbind',arbitrary_list)
  finalResults=merge(modelResults,ModelApplicability_obj,by="row.names")
  #else{
    #ModelApplicability_obj = ModelApplicability(CalPredsModelApplicability,
    #                                            modelId = modelID,
    #                                            applicabilitypreds) # add to config file or add an R object with calpreds

    #finalResults=merge(modelResults,ModelApplicability_obj,by="row.names")

#}


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
    if (nrow(sampleOEs[sampleOEs$modelId %in% 25:26,])>=1){
      ecoregion=sf::st_read(paste0(ecoregion_base_path,"GIS//GIS_Stats//CONUS//ecoregion//hybrid_level_III_ecoregions.shp"))
      ecoregion=sf::st_make_valid(ecoregion)
      finalResults_sf=sf::st_transform(finalResults_sf,5070)
      finalResults_sf=sf::st_intersection(finalResults_sf,ecoregion)
      finalResults=dplyr::left_join(finalResults,sf::st_drop_geometry(finalResults_sf[,c('sampleId','modelId')]),by='sampleId')
    }else if (nrow(sampleOEs[sampleOEs$modelId %in% 13:23,])>=1){
      ecoregion=sf::st_read(paste0(ecoregion_base_path,"GIS//GIS_Stats//Wyoming//ecoregion//BIOREGIONS_2011_modifiedCP.shp"))
      ecoregion=sf::st_make_valid(ecoregion)
      finalResults_sf=sf::st_transform(finalResults_sf,5070)
      finalResults_sf=sf::st_intersection(finalResults_sf,ecoregion)
      finalResults=dplyr::left_join(finalResults,sf::st_drop_geometry(finalResults_sf[,c('sampleId','modelId')]),by='sampleId')
    } else{
      finalResults$modelId=modelID
    }

    #[need to edit for new code]
    #little if statements that subset out criticial or not critical salmonid habitats
    #for AIM ID

    if(modelID==25 & unique(def_samples$customerAbbreviation)=='BLM-AIM' & unique(def_samples$usState=="Idaho")){
      finalResults_tosub<-finalResults
      intermediate<-plyr::join(def_samples[,c('sampleId','siteId'),],IDsites[,c('siteId','waterbodyCode')],by='siteId','left')
      finalResults_tosub<-plyr::join(finalResults_tosub,intermediate,by='sampleId','left')
      finalResults<-finalResults_tosub[which(finalResults_tosub$siteId %in% NonCrit$siteId),]
    }
    if(modelID==9 & unique(def_samples$customerAbbreviation)=='BLM-AIM'){
      finalResults_tosub<-finalResults
      intermediate<-plyr::join(def_samples[,c('sampleId','siteId'),],IDsites[,c('siteId','waterbodyCode')],by='siteId','left')
      finalResults_tosub<-plyr::join(finalResults_tosub,intermediate,by='sampleId','left')
      finalResults<-finalResults_tosub[which(finalResults_tosub$siteId %in% CritHab$siteId),]
    }
    write.csv(finalResults,paste0('finalresults_',paste(modelID,collapse='_'),"_",Sys.Date(),'.csv'))
    finalResults
