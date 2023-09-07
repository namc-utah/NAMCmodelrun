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

##Run models code, originally developed by Jennifer Courtwright
#and edited by Andrew Caudillo for multiple modelID ingestion

#the workflow of this code is as follows:

#1. Read in basic model info (IDs, etc)
#2. subset the sampleIds via their respective models
#there may be multiple models assigned to a sampleID and that is ok!
#3. Generate base Model applicability info for easy access later
#4. Generate model metadata for each model to be run
#5. Read in support functions for RF objects etc.
#6. Get predictor values for the model(s) needed
#7. Start getting the bug information ("bugnew" in the original code)
#AND run the indices here, too. This step also joins spatial data
#and runs model applicability
#and stores the information in a list, with each element
#referring to a certain modelId
#if there is only one modelId, the list is coerced into a data frame.
#Section 7 breaks the models out by type. OE, MMI, Hybrid, etc.

#model types
#1 - OE
#2 - MMI
#3 - CSCI
#4- WQ
basic_models = query(
  api_endpoint = "models")

OEs<-basic_models$modelId[basic_models$modelTypeId==1 & basic_models$modelId !=12]
MMIs<-basic_models$modelId[basic_models$modelTypeId==2]
CSCIs<-basic_models$modelId[basic_models$modelTypeId==3]
WQs<-basic_models$modelId[basic_models$modelTypeId==4]
AREMP<-c(7,8)
NullOE<-12
MIR<-136
PIBO <-9


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
  sampleIds=sampleIds
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

All_OE_results<-list()
ALL_MMI_results<-list()
ALL_CSCI_results<-list()
ALL_MIR_results<-list()
ALL_WQ_Results<-list()

#get base model applicability info

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


#def_predictors <- def_predictors[!duplicated(def_predictors), ]

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
#edit this for CO only
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

# get predictors into wide format needed for model functions
prednew = tidyr::pivot_wider(def_predictors,
                             id_cols="sampleId",
                             names_from = "abbreviation",
                             values_from = "predictorValue")# add id_cols=sampleId once it gets added to end point
prednew=as.data.frame(prednew)
rownames(prednew)<-prednew$sampleId
prednew<-prednew[,-1]
}


# ---------------------------------------------------------------
# Run models
# ---------------------------------------------------------------

#quick export section (CO MMI and ADEQ IBI)
#no analyses performed by NAMC for these indices
#so we will just export the appropriate file
#for MS Access ingestion.
 if (nrow(sampleMMIs[sampleMMIs$modelId %in% c(169),])>=1){
   setwd('C://Users//andrew.caudillo//Box//NAMC//OEModeling//NAMC_Supported_OEmodels//Arizona//InputFiles')
  ADEQ_bug_export(sampleIds = def_model_results$sampleId)
}
# ------------------------------
# OE models
# ------------------------------
# models using john vansickles RIVPACS random forest code : AREMP, UTDEQ15, Westwide, PIBO

#AND collect bug into simultaneously

if(nrow(sampleNullOEs)>=1){
  Nullbugnew = OR_NBR_bug(
    sampleIds = sampleNullOEs$sampleId,
    translationId = 14,
    fixedCount = 300
  )
  #run the model right here!
  #NBR PREDATOR doest not get model applicability.
  NullOR_modelResults <- OR_NBR_model(Nullbugnew)
NullOR_modelResults$modelId<-rep(12,nrow(NullOR_modelResults))
  for(n in 1:nrow(NullOR_modelResults)){
    #this yields split count. Should be rarefied to 300, but if we have
    #a sample with fewer than 300, it will yield that number instead.
    taxa_counts = query(
      api_endpoint = "sampleTaxaTranslationRarefied",
      args = list(translationId = 16, fixedCount = 300, sampleIds=NullOR_modelResults$sampleId[n])
    )
    agged_tax_counts<-aggregate(splitCount~sampleId,data=taxa_counts,FUN=sum)
    #getting invasives
    bugsOTU = NAMCr::query("sampleTaxaTranslationRarefied",
                           translationId = 14,
                           fixedCount = 300,
                           sampleIds=sampleNullOEs$sampleId
    )
    sumrarefiedOTUTaxa = bugsOTU  %>%
      dplyr::group_by(sampleId) %>%
      dplyr::summarize(fixedCount = sum(splitCount))

    ################################################
    ###### get invasives ##### comment out this entire invasives section if running "National" AIM westwide reporting
    # get raw bug data
    bugRaw = NAMCr::query(
      "sampleTaxa",
      sampleIds=sampleNullOEs$sampleId
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
    NullOR_modelResults_burn=dplyr::left_join(NullOR_modelResults,additionalbugmetrics,by="sampleId")
    n_dat_to_pass<-list(sampleId = NullOR_modelResults_burn$sampleId[n],
                      modelId = NullOR_modelResults_burn$modelId[n],
                      oResult = NullOR_modelResults_burn$O[n],
                      eResult = NullOR_modelResults_burn$E[n],
                      modelResult = NullOR_modelResults_burn$OoverE[n],
                      fixedCount=agged_tax_counts$splitCount,
                      notes=NullOR_modelResults_burn$InvasiveInvertSpecies[n])

    NAMCr::save(
     api_endpoint = "setModelResult",
      args=n_dat_to_pass)
    print('results saved!')
    print(n)
  }

}

#get bug info for the non-null models
#we will call this data in each specific loop
#for each modelID.
if(nrow(sampleOEs)>=1){
  print('Running O/E')
  OE_list<-list()
  for(i in 1:length(unique(sampleOEs$modelId))){
    uniq_mod<-unique(sampleOEs$modelId)[i]
    mod_val<-def_models[def_models$modelId==uniq_mod,]
    y<-sampleOEs[sampleOEs$modelId==unique(sampleOEs$modelId)[i],]
  bugnew = OE_bug_matrix(
    sampleIds =y$sampleId,
    translationId = mod_val$translationId[1],
    fixedCount = mod_val$fixedCount[1])
  OE_list[[i]]<-bugnew
  names(OE_list)[i]<-unique(sampleOEs$modelId)[i]
  }
}

#non OR O/E indices (WW, PIBO, etc.)
  if (nrow(sampleOEs[sampleOEs$modelId %in% c(2,7,9,25,26,29),])>=1) {
    print('running O/Es using 4.2RF')
    print(unique(sampleOEs$modelId[which(sampleOEs$modelId %in% c(2,7,9,25,26,29))]))
    bug_sub_list<-sampleOEs[which(sampleOEs$modelId %in% c(2,7,9,25,26,29)),]
    n_unique_OE_mods<-length(unique(bug_sub_list$modelId))
    print(n_unique_OE_mods)
    General_OE_results<-list()


    for(i in 1:n_unique_OE_mods){
      print('getting OE results')
      model_id_burn<-as.character(unique(bug_sub_list$modelId)[i])
      model_sub<-bug_sub_list[bug_sub_list$modelId==as.integer(model_id_burn),]
      oe_bug_burn<-OE_list[[model_id_burn]]
      #colnames(oe_bug_burn)<-sub(c('X2.','X7.','X9.','X25.','X26.','X29.'),colnames(oe_bug_burn))
      OE <-model.predict.RanFor.4.2(
        bugcal.pa,
        grps.final,
        preds.final,
        ranfor.mod,
        prednew,
        bugnew=oe_bug_burn,
        Pc = 0.5,
        Cal.OOB = FALSE)


      modelResults<-OE$OE.scores[(grepl("NA", row.names(OE$OE.scores), fixed = TRUE))==F,]
      #modelResults$modelID<-as.integer(model_id_burn)
      General_OE_results[[i]]<-modelResults
      names(General_OE_results)[i]<-model_id_burn
      print('model app time')

      ModelApplicability_obj = ModelApplicability(CalPredsModelApplicability,
                                                  modelId = as.integer(model_id_burn),
                                                  applicabilitypreds)

      General_OE_results[[i]]<-merge(General_OE_results[[i]],
                                          ModelApplicability_obj,
                                     by='row.names')

      General_OE_results[[i]]<-dplyr::left_join(General_OE_results[[i]],def_samples[,c('sampleId','siteLongitude','siteLatitude')],by='sampleId')

      GeneralOEResults_sf<-sf::st_as_sf(General_OE_results[[i]],coords=c('siteLongitude','siteLatitude'),crs=4269)
      if (model_id_burn %in% c("25","26")){
        print('Generating ModelId for a WestWide O/E model')
        ecoregion=sf::st_read(paste0(ecoregion_base_path,"GIS//GIS_Stats//CONUS//ecoregion//hybrid_level_III_ecoregions.shp"))
        ecoregion=sf::st_make_valid(ecoregion)
        GeneralOEResults_sf=sf::st_transform(GeneralOEResults_sf,5070)
        GeneralOEResults_sf=sf::st_intersection(GeneralOEResults_sf,ecoregion)
        General_OE_results[[i]]=dplyr::left_join(General_OE_results[[i]],sf::st_drop_geometry(GeneralOEResults_sf[,c('sampleId','modelId')]),by='sampleId')
      } else{
      print('not a westwide model, going to other options (WY, OR)')
      }
      General_OE_results[[i]]$modelId=as.integer(model_id_burn)}

      bugsOTU = NAMCr::query("sampleTaxaTranslationRarefied",
                             translationId = model_sub$translationId[i],
                             fixedCount = model_sub$fixedCount[i],
                             sampleIds=model_sub$sampleId
      )
      sumrarefiedOTUTaxa = bugsOTU  %>%
        dplyr::group_by(sampleId) %>%
        dplyr::summarize(fixedCount = sum(splitCount))

      ################################################
      ###### get invasives ##### comment out this entire invasives section if running "National" AIM westwide reporting
      # get raw bug data
      bugRaw = NAMCr::query(
        "sampleTaxa",
        sampleIds=model_sub$sampleId
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
      General_OE_results[[i]]=dplyr::left_join(General_OE_results[[i]],additionalbugmetrics,by="sampleId")
      # finalResults=dplyr::left_join(finalResults,sumrarefiedOTUTaxa,by="sampleId")
      # finalResults$InvasiveInvertSpecies='National'

      print('writing O/E results')
      for(j in 1:nrow(model_sub)){
      dat_to_pass<-list(sampleId = General_OE_results[[i]]$sampleId[j],
                        modelId = General_OE_results[[i]]$modelId[j],
                        oResult = General_OE_results[[i]]$O[j],
                        eResult = General_OE_results[[i]]$E[j],
                        modelResult = General_OE_results[[i]]$OoverE[j] ,
                        fixedCount = General_OE_results[[i]]$fixedCount[j],
                        modelApplicability = General_OE_results[[i]]$ModelApplicability[j],
                        notes=General_OE_results[[i]]$InvasiveInvertSpecies[j])

      NAMCr::save(
        api_endpoint = "setModelResult",
        args=dat_to_pass)
      }
    }

#if only one model met that condition above
#force the list to a df. maybe keeping it a list is better... for looping
#the output would work with the latter option.

    #if(length(General_OE_results)==1)
    #  General_OE_results<-General_OE_results[[1]]


#this section if for non PREDATOR O/E indices
    if (nrow(sampleOEs[sampleOEs$modelId %in% 10:11,])>=1) {
      print('running O/E using 4.1RF')
      bug_sub_list<-sampleOEs[sampleOEs$modelId %in% 10:11,]
      model_id_burn<-as.character(unique(bug_sub_list$modelId)[i])
      model_sub<-bug_sub_list[bug_sub_list$modelId==as.integer(model_id_burn),]
      print(unique(sampleOEs$modelId[sampleOEs$modelId %in% 10:11]))
      n_unique_OE_mods<-length(unique(bug_sub_list$modelId))
      print(n_unique_OE_mods)
      OR_OE_result<-list()
      for(i in 1:n_unique_OE_mods){
        print('getting OR OE results')
        uniq_mod<-unique(sampleOEs$modelId)[i]
        mod_val<-def_models[def_models$modelId==uniq_mod,]
        model_id_burn<-as.character(unique(bug_sub_list$modelId)[i])
        oe_bug_burn<-OE_list[[model_id_burn]]
        pred_burn<-prednew[row.names(prednew) %in% row.names(oe_bug_burn),]
        OE <-model.predict.v4.1(bugcal.pa,
                                grps.final,
                                preds.final,
                                grpmns,
                                covpinv,
                                prednew=pred_burn,
                                bugnew=oe_bug_burn,
                                Pc = 0.5)
        modelResults<-OE$OE.scores
        modelResults<-OE$OE.scores[(grepl("NA", row.names(OE$OE.scores), fixed = TRUE))==F,]
        OR_OE_result[[i]]<-modelResults

        names(OR_OE_result)[i]<-model_id_burn
        print('model app time')

        ModelApplicability_obj = ModelApplicability(CalPredsModelApplicability,
                                                    modelId = as.integer(model_id_burn),
                                                    applicabilitypreds)
        print("Model app done")
        OR_OE_result[[i]]<-merge(OR_OE_result[[i]],
                                       ModelApplicability_obj,
                                       by='row.names')

        OR_OE_result[[i]]<-dplyr::left_join(OR_OE_result[[i]],def_samples[,c('sampleId','siteLongitude','siteLatitude')],by='sampleId')



        OR_OE_result[[i]]$modelId=as.integer(model_id_burn)}

      bugsOTU = NAMCr::query("sampleTaxaTranslationRarefied",
                             translationId = mod_val$translationId,
                             fixedCount = mod_val$fixedCount,
                             sampleIds=model_sub$sampleId
      )
      sumrarefiedOTUTaxa = bugsOTU  %>%
        dplyr::group_by(sampleId) %>%
        dplyr::summarize(fixedCount = sum(splitCount))

      ################################################
      ###### get invasives ##### comment out this entire invasives section if running "National" AIM westwide reporting
      # get raw bug data
      bugRaw = NAMCr::query(
        "sampleTaxa",
        sampleIds=model_sub$sampleId
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
      OR_OE_result[[i]]=dplyr::left_join(OR_OE_result[[i]],additionalbugmetrics,by="sampleId")
      # finalResults=dplyr::left_join(finalResults,sumrarefiedOTUTaxa,by="sampleId")
      # finalResults$InvasiveInvertSpecies='National'

      print('writing O/E results')


      for(j in 1:nrow(model_sub)){
        dat_to_pass<-list(sampleId = OR_OE_result[[i]]$sampleId[j],
                          modelId = OR_OE_result[[i]]$modelId[j],
                          oResult = OR_OE_result[[i]]$O[j],
                          eResult = OR_OE_result[[i]]$E[j],
                          modelResult = OR_OE_result[[i]]$OoverE[j] ,
                          fixedCount = OR_OE_result[[i]]$fixedCount[j],
                          modelApplicability = OR_OE_result[[i]]$ModelApplicability[j],
                          notes=OR_OE_result[[i]]$InvasiveInvertSpecies[j])

        NAMCr::save(
          api_endpoint = "setModelResult",
          args=dat_to_pass)
      }
    }



##WY O/Es
#these will never have just one model ID, so no need to
#coerce to a single dataframe, unlike the O/Es before.

 if (nrow(sampleOEs[sampleOEs$modelId %in% 13:23,])>=1) { #if T
    print ('Running OE for WY models')
    WY_sub_list<-sampleOEs[sampleOEs$modelId %in% 13:23,]
    n_unique_WYmods<-length(unique(WY_sub_list$modelId))
    ALK_LOG = setNames(as.data.frame(
    log10(predict(ranfor.mod, prednew, type = "response"))
      ), c("ALK_LOG"))
      prednew = cbind(prednew, ALK_LOG)
      WY_OE_results<-list()
        ## need to subset to only model predictors or maybe it doesnt matter??
        for(i in 1:n_unique_WYmods){ #for i
          print('getting WY model id')
          model_id_burn<-as.character(unique(WY_sub_list$modelId)[i])
          print(model_id_burn)
          model_sub<-WY_sub_list[WY_sub_list$modelId==as.integer(model_id_burn),]
          oe_bug_burn<-OE_list[[model_id_burn]]
          pred_burn<-prednew[row.names(prednew) %in% row.names(oe_bug_burn),]
          OE <-model.predict.v4.1(bugcal.pa,
                                  grps.final,
                                  preds.final,
                                  grpmns,
                                  covpinv,
                                  prednew=pred_burn,
                                  bugnew=oe_bug_burn,
                                  Pc = 0.5)
          print('OE done')
          modelResults<-OE$OE.scores
          WY_OE_results[[i]]<-modelResults
          names(WY_OE_results)[i]<-model_id_burn
          if(model_id_burn!="13")
            if(0){

          print('model app time')

          ModelApplicability_obj = ModelApplicability(CalPredsModelApplicability,
                                                      modelId = as.integer(model_id_burn),
                                                      applicabilitypreds)

         WY_OE_results[[i]]<-merge(WY_OE_results[[i]],
                                         ModelApplicability_obj,
                                         by='row.names')
            }
          WY_OE_results[[i]]$sampleId<-as.integer(row.names(WY_OE_results[[i]]))

          WY_OE_results[[i]]<-dplyr::left_join(WY_OE_results[[i]],def_samples[,c('sampleId','siteLongitude','siteLatitude')],by='sampleId')
          #using the extraction method in sf returns too many
          #modelIds for WY because of the strong overlap in the regions
          #use what INSTAR has assigned to it and skip over this section.
          #assign modelId via model_id_burn
if(0){
          WYOEResults_sf<-sf::st_as_sf(WY_OE_results[[i]],coords=c('siteLongitude','siteLatitude'),crs=4269)
          #fix to work with above code
            ecoregion=sf::st_read(paste0(ecoregion_base_path,"GIS//GIS_Stats//Wyoming//ecoregion//BIOREGIONS_2011_modifiedCP.shp"))
            ecoregion=sf::st_make_valid(ecoregion)
            WY_OE_results_sf=sf::st_transform(WYOEResults_sf,5070)
            WYOEresults_sf=sf::st_intersection(WY_OE_results_sf,ecoregion)
            WYOEresults_sf<-WYOEresults_sf[WYOEresults_sf$modelId==as.integer(model_id_burn),]

            WY_OE_results[[i]]=dplyr::left_join(WY_OE_results[[i]],sf::st_drop_geometry(WYOEresults_sf[,c('sampleId','modelId')]),by='sampleId')
}
            bugsOTU = NAMCr::query("sampleTaxaTranslationRarefied",
                                   translationId = mod_val$translationId,
                                   fixedCount = mod_val$fixedCount,
                                   sampleIds=model_sub$sampleId
            )
            sumrarefiedOTUTaxa = bugsOTU  %>%
              dplyr::group_by(sampleId) %>%
              dplyr::summarize(fixedCount = sum(splitCount))

            ################################################
            ###### get invasives ##### comment out this entire invasives section if running "National" AIM westwide reporting
            # get raw bug data
            bugRaw = NAMCr::query(
              "sampleTaxa",
              sampleIds=model_sub$sampleId
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
            WY_OE_results[[i]]=dplyr::left_join(WY_OE_results[[i]],additionalbugmetrics,by="sampleId")
            # finalResults=dplyr::left_join(finalResults,sumrarefiedOTUTaxa,by="sampleId")
            # finalResults$InvasiveInvertSpecies='National'

            print('writing O/E results')



            for(j in 1:nrow(WY_OE_results[[i]])){ #for j
              if(model_id_burn=="13"){
                print('only WY model with applicability attached')
              dat_to_pass<-list(sampleId = WY_OE_results[[i]]$sampleId[j],
                                modelId = WY_OE_results[[i]]$modelId[j],
                                oResult = WY_OE_results[[i]]$O[j],
                                eResult = WY_OE_results[[i]]$E[j],
                                modelResult = WY_OE_results[[i]]$OoverE[j] ,
                                fixedCount = WY_OE_results[[i]]$fixedCount[j],
                                modelApplicability = WY_OE_results[[i]]$ModelApplicability[j],
                                notes=WY_OE_results[[i]]$InvasiveInvertSpecies[j])
              }else{ #else
                print('all other WY models here')
                dat_to_pass<-list(sampleId = WY_OE_results[[i]]$sampleId[j],
                                  modelId = as.integer(model_id_burn),
                                  oResult = WY_OE_results[[i]]$O[j],
                                  eResult = WY_OE_results[[i]]$E[j],
                                  modelResult = WY_OE_results[[i]]$OoverE[j] ,
                                  fixedCount = WY_OE_results[[i]]$fixedCount[j],
                                  notes=WY_OE_results[[i]]$InvasiveInvertSpecies[j])
              }#else

              NAMCr::save(
                api_endpoint = "setModelResult",
                args=dat_to_pass)
            }# for j


        } #for i
      print(i)
} #if T

  #Starting off the MMI section is the CO MMI, a quick and easy export
if (nrow(sampleMMIs[sampleMMIs$modelId %in% c(4,5,6),])>=1) {
    print('CO MMI, this is run in Access. Only predictors and bugs needed')
    #write get bugs from database, write out as a csv and save as CObugs object
    CObugs=CO_bug_export(sampleIds=sampleIds)
    #write out predictors as a csv
    write.csv(prednew,file = paste0("COpredictors","boxId_",CObugs$Project[1],"_",Sys.Date(),".csv"),row.names=FALSE)
    cat(paste("csv with COpredictors has been written out to your current working directory.",
              "Convert this csv to excel 2003 and import into CO EDAS access database to compute the CSCI score.",
              "Follow instructions in this pdf Box\\NAMC\\OE_Modeling\\NAMC_Supported_OEmodels\\CO\\Documentation\\EDAS2017\\Tutorial Guide to EDAS_Version 1.7.pdf",
              "to import bug and habitat data, harmonize taxa list, rarefy and compute MMI",
              "then read resulting excel file back into R to save results in the database.", sep="\n"))
    #throwing in CSCI and MIR here, just because they are small
} else if (nrow(sampleCSCIs)>=1) {
  bugnew = CSCI_bug(sampleIds = def_model_results$sampleId)
} else if (nrow(sampleMMIs)>=1) {

#Now the MMI section starts.
#Make MMI list, like OE list
#and fill with modelIDs

  MMI_list<-list()
  for(i in 1:length(unique(sampleMMIs$modelId))){
  m<-sampleMMIs[sampleMMIs$modelId==unique(sampleMMIs$modelId)[i],]
  bugnew = MMI_metrics(sampleIds = m$sampleId,
                       translationId=m$translationId[1],
                       fixedCount = m$fixedCount[1],
                       modelId=m$modelId[1])
  MMI_list[[i]]<-bugnew
  names(MMI_list)<-as.character(m$modelId[1])
  }
}
#AREMP MMI is a little special because of metrics
#create its results here
if (nrow(sampleAREMP[sampleAREMP$modelId==8])>=1) {# AREMP MMI
  AREMP_bugnew<-MMI_list[["8"]]
  modelResults <-
    AREMP_MMI_model(
      AREMP_bugnew,
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

}else if (nrow(sampleMMIs[sampleMMIs$modelId==3,])>=1) {
  # NV MMI
  # need to call conductivity model first before calling the NV model because predicted conductivity is a predictor for the NV model
  NV_bugnew<-MMI_list[["3"]]
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

}
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
  #getting results  here, since this is the last time WQ is queried
  WQmodelResults = as.data.frame(predict(ranfor.mod, prednew, type = "response"))# make sure prednew has sampleIds as the rows
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

#renaming columns to match database
if (nrow(sampleAREMP)>=1){
  prednew$TMAX_WS=prednew$Tmax_WS
}
# same issue with RH_WS being included in TP model but different capitalization. However, revisit if we revise TP model
if (nrow(sampleAREMP[sampleAREMP$modeId==8])>=1){
  prednew$rh_WS=prednew$RH_WS
}

  #prep data for CSCI, if needed
if (nrow(sampleCSCIs)>=1) {
  prednew$sampleId=as.numeric(row.names(prednew))
  prednew=left_join(prednew,def_samples[,c('siteName','sampleId')],by='sampleId')
  prednew$StationCode<-prednew$sampleId
  report <- CSCI::CSCI(bugs = bugnew, stations = prednew)
  modelResults = report$core
  rownames(modelResults)=modelResults$SampleID


#Modeled Insect Richness (an effort largely created by Jennifer Courtwright)

} else if (nrow(sampleMIRs)>=1){#modeled insect richness
     UniqueRichness_Insecta_pred=predict(rfmod_UniqueRichness_Insecta, prednew, type="response")
     # join predicted insect richness to intial data
     modelResults=cbind(bugnew, UniqueRichness_Insecta_pred)
     # convert to O//E ratio
     modelResults$modeledInsectRichness=modelResults$UniqueRichness_Insecta/modelResults$UniqueRichness_Insecta_pred

}



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


    ####

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

