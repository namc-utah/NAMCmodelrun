# if any predictors not valid look for them in model results (predicted coductivity and alalinity

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
  modelId = modelId
)

# ---------------------------------------------------------------
# Get bug data for model functions
# ---------------------------------------------------------------
#
# if modelType= bug OE get OTU taxa matrix
if (def_models$modelId == 12) {#need a way to distinguish this model from others.. call NULL OE?
  bugnew = OR_NBR_bug_box(
    boxId = boxId,
    translationId = def_models$translationId,
    fixedCount = def_models$fixedCount)
} else if(def_models$modelTypeAbbreviation == "OE") {
  bugnew = OE_bug_matrix_box(
    boxId = boxId,
    translationId = def_models$translationId,
    fixedCount = def_models$fixedCount
  )
} else if (def_models$modelId %in% c(4, 5, 6)) {# CO model must be written out as an excel file using a separate bank of code and function
  CObugs=CO_bug_export_box(boxId = boxId)
} else if (def_models$modelId %in% c(1)) {# CSCI requires just the raw taxa list translated for misspelling
  bugnew = CSCI_bug_box(boxId = boxId)
}else if (def_models$modelTypeAbbreviation == "MMI") {# if modelType= bug MMI get
  bugnew = MMI_metrics_box(boxId = boxId, translationId=def_models$translationId, fixedCount = def_models$fixedCount)
} else {

}
bugnew<-subset(bugnew,rownames(bugnew) %in% rownames(prednew))
prednew<-subset(prednew,rownames(prednew) %in% rownames(bugnew))
#reorder them the same just in case model functions dont already do this
bugnew = bugnew[order(rownames(bugnew)),];
prednew = prednew[order(rownames(prednew)),];

# ---------------------------------------------------------------
# load model specific R objects which include reference bug data and predictors RF model objects
# ---------------------------------------------------------------
# every model has an R object that stores the random forest model and reference data
# the R objects are named with the model abbreviation
# instead of all these if statements the R file name could be stored in the database... and should be!!
#if CO, CSCI, or OR null model no R data file needs loaded in
if (def_models$modelId %in% c(1,4,5,6,12)){
print("no R object needs loaded")
  #if WY model only one Rdata file needs loaded and not one for each "model" but Alkalinity also needs added
} else if (def_models$modelId %in% c(13:23)){
  load("sysdata.rda/WY2018.Rdata")
  load("sysdata.rda/Alkalinity.Rdata")
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
# models using latest version of van sickle function include: AREMP, UTDEQ15, Westwide, PIBO
if (def_models$modelId %in% c(7, 2, 25, 26, 9)) {
  OE <-
    model.predict.RanFor.4.2(
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
} else if (def_models$modelId %in% c(10, 11)) {# models using version 4.1 of van sickle code include: OR_WCCP, OR_MWCF
  OE <-
    model.predict.v4.1(bugcal.pa,
                       grps.final,
                       preds.final,
                       grpmns,
                       covpinv,
                       prednew,
                       bugnew,
                       Pc = 0.5)# add elpsis...
  modelResults<-OE$OE.scores
}else if (def_models$modelId %in% (13:23)) {# WY also uses version 4.1 of van sickle code but requires alkalinity model as a dependency
  ALK_LOG = setNames(as.data.frame(
    predict(ranfor.mod, prednew, type = "response")), c("ALK_LOG"))# need to log value
  prednew = cbind(prednew, ALK_LOG)
  OE <-
    model.predict.v4.1(bugcal.pa,
                       grps.final,
                       preds.final,
                       grpmns,
                       covpinv,
                       prednew,
                       bugnew,
                       Pc = 0.5)
  modelResults<-OE$OE.scores
}else if (def_models$modelId == 12) {# OR eastern region is a null model and no predictors are used
  modelResults <- OR_NBR_model(bugnew)

}else if (def_models$modelId == 1) {# CSCI has its own package and function
    report <- CSCI::CSCI(bugs = CSCIbugs, stations = prednew)
    modelResults = report$core
    rownames(modelResults)=modelResults$SampleID

     }else if (def_models$modelId == 8) {
  ## all MMIs will need their own function added here because there is a rf model for each metric
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
}  else if (def_models$modelId == 3) {# all MMIs will need their own function added here because there is a rf model for each metric
  # need to call conductivity model first before calling the NV model because predicted conductivity is a predictor for the NV model
  PrdCond = setNames(as.data.frame(
    predict(ranfor.mod, prednew, type = "response")
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
}else if (def_models$modelId %in% c(27, 28, 29, 30)) {#conductivity, tp, tn,temperature
  modelResults = as.data.frame(predict(ranfor.mod, prednew, type = "response"))# make sure prednew has sampleIds as the rows
}else{

}



# ---------------------------------------------------------------
# Always run model applicability test
# ---------------------------------------------------------------
# get site locations from database
def_sites = NAMCr::query(
  api_endpoint = "samples",
  include = c("sampleId","siteId", "siteName", "usState", "siteLocation","customerSiteCode","visitId"),
  boxId = boxId
)
points2process=geojsonsf::geojson_sf(def_sites$siteLocation)

# use nhdplusTools to get COMID for each site
for (p in 1:nrow(points2process)){
  tryCatch({start_comid <- nhdplusTools::discover_nhdplus_id(points2process[p,1])},error=function(e){
  cat(paste0("error_row_",p))
})
  if(p == 1){
    comids= start_comid
  }else{
    comids = rbind(comids, start_comid)
  }
}
def_sites$COMID=comids


# get elevation, watershed area, precip, and temp variable from stream cat SQLite database
# create function to loop over each site for input into SQL query
inLOOP<- function(inSTR,...) {
  inSTR=unlist(inSTR)
  if (inSTR[1]==''){loopSTR="''"} else{
    for (i in 1:length(inSTR)){
      comma=ifelse(i==length(inSTR),'',',')
      STRl=sprintf("'%s'%s",inSTR[i],comma)
      if(i==1){loopSTR=STRl} else{loopSTR=paste(loopSTR,STRl)}
    } }
  return(loopSTR)
}
#create database connection
conn<-DBI::dbConnect(RSQLite::SQLite(),SQLite_file_path)
preds = DBI::dbGetQuery(conn,sprintf("SELECT COMID, ElevCat,	Precip8110Ws,	Tmean8110Ws,	WsAreaSqKm FROM StreamCat_2022 WHERE COMID in (%s)",inLOOP(substr(def_sites$COMID,1,10))))
applicabilitypreds=merge(def_sites,preds, by="COMID")

# run model applicability function
ModelApplicability = ModelApplicability(CalPredsModelApplicability,
                                        modelId = modelId,
                                        applicabilitypreds) # add to config file or add an R object with calpreds

finalResults=merge(modelResults,ModelApplicability,by="row.names")




# ---------------------------------------------------------------
# Get additional bug metrics (fixed count and invasives)
# ---------------------------------------------------------------
##### get fixed count column #####
bugsOTU = NAMCr::query("sampleTaxaTranslationRarefied",
                       translationId = def_models$translationId,
                       fixedCount = def_models$fixedCount,
                       boxId = boxId
)
sumrarefiedOTUTaxa = bugsOTU  %>%
  dplyr::group_by(sampleId) %>%
  dplyr::summarize(fixedCount = sum(splitCount))


###### get invasives #####
# get raw bug data
bugRaw = NAMCr::query(
  "sampleTaxa",
  boxId = boxId
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




# ---------------------------------------------------------------
# join all together and write out final file
# ---------------------------------------------------------------

finalResults=dplyr::left_join(finalResults,additionalbugmetrics,by="sampleId")

# subset only the columns we need
if (def_models$modelTypeAbbreviation == "OE") {
  finalResults=finalResults[,c("sampleId","visitId","customerSiteCode","O","E","OoverE","ModelApplicability","fixedCount","InvasiveInvertSpecies")]
}else if (def_models$modelTypeAbbreviation == "Hybrid"){
  finalResults=finalResults[,c("sampleId","visitId","customerSiteCode","OoverE","MMI","CSCI","ModelApplicability","fixedCount","InvasiveInvertSpecies")]
}else if (def_models$modelTypeAbbreviation == "MMI"){
 finalResults=finalResults[,c("sampleId","visitId","customerSiteCode","MMI","ModelApplicability","fixedCount","InvasiveInvertSpecies")]
}else{
}

write.csv(finalResults,paste0("FinalResults_",def_models$abbreviation,"_",Sys.Date(),".csv"),row.names = FALSE)




  # ---------------------------------------------------------------
  # Save model results
  # ---------------------------------------------------------------
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
    }else if (def_models$modelTypeAbbreviation == "Hybrid") {
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
    cat(paste0("\n\tSAMPLE ERROR: ",finalResults$sampleId[i],"\n"))
    str(e,indent.str = "   "); cat("\n")
  })
}

