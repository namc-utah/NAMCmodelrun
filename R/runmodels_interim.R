library(randomForest)
library(NAMCr)
library(tidyverse)
library(DBI)
library(RSQLite)
#library(CSCI)
boxId=2150
modelId=25
#prednew=read.csv()
SQLite_file_path="C:/NAMC_S3/StreamCat/StreamCat.sqlite"
source("R/Bug_functions_box.R")
source("R/Model_functions.R")
source("R/model.predict.RanFor.4.2.r")
source("R/model.predict.RanFor.r")
source("R/model.predict.v4.1.r")

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
if (def_models$modelTypeAbbreviation == "OE") {
  bugnew = OE_bug_matrix_box(
    boxId = boxId,
    translationId = def_models$translationId,
    fixedCount = def_models$fixedCount
  )
} else if (def_models$modelId == 12) {#need a way to distinguish this model from others.. call NULL OE?
  bugnew = OR_NBR_bug_box(
    boxId = boxId,
    translationId = def_models$translationId,
    fixedCount = def_models$fixedCount
  )
} else if (def_models$modelId %in% c(4, 5, 6)) {# CO model must be written out as an excel file using a separate bank of code and function
  CObugs=CO_bug_export_box(boxId = boxId)
} else if (def_models$modelId %in% c(1)) {# CSCI requires just the raw taxa list translated for misspelling
  # add in model names in comments
  CSCIbugs = CSCI_bug_box(boxId = boxId)
}else if (def_models$modelTypeAbbreviation == "MMI") {# if modelType= bug MMI get
  bugnew = MMI_metrics_box(boxId = boxId, fixedCount = def_models$fixedCount)
} else {

}

# ---------------------------------------------------------------
# load model specific R objects which include reference bug data and predictors RF model objects
# ---------------------------------------------------------------
# every model has an R object that stores the random forest model and reference data
# the R objects are named with the model abbreviation
# instead of all these if statements the R file name could be stored in the database... and should be!!
#if CO or CSCI model no R data file needs loaded in
if (def_models$modelId %in% c(1,4,5,6)){

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
# models using latest version of van sickle function include: AREMP, UTDEQ15, Westwide
if (def_models$modelId %in% c(7, 2, 25, 26)) {
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
}else if (def_models$modelId %in% (13:23)) {# WY also uses version 4.1 of van sickle code but requires alkalinity model as a dependency
  ALK_LOG = setNames(as.data.frame(
    predict(ranfor.mod, prednew, type = "response")
  ), c("ALK_LOG"))# need to log value
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
}else if (def_models$modelId == 9) {# PIBO model was one of the earliest models built and used a version of van sickle's function that didnt have a version number
  OE <-
    PIBO_model(bug.otu,
               bugall,
               grps.final,
               preds.final,
               ranfor.mod,
               prednew)
}else if (def_models$modelId == 12) {# OR eastern region is a null model and no predictors are used
  OE <- OR_NBR_model(bugnew)

}else if (def_models$modelId == 1) {# CSCI has its own package and function
    report <- CSCI::CSCI(bugs = CSCIbugs, stations = prednew)
    OE = report$core

     }else if (def_models$modelId == 8) {
  ## all MMIs will need their own function added here because there is a rf model for each metric
  MMI <-
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
  MMI <-
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
  WQ = as.data.frame(predict(ranfor.mod, prednew, type = "response"))# make sure prednew has sampleIds as the rows
}else{

}



# ---------------------------------------------------------------
# Always run model applicability test
# ---------------------------------------------------------------
# get all predictor values needed for a box or project # note this either needs a loop written over it or a different API end point
conn<-DBI::dbConnect(RSQLite::SQLite(),SQLite_file_path)
applicabilitypreds = DBI::dbGetQuery(conn,sprintf("SELECT ElevCat,	Precip8110Ws,	Tmean8110Ws,	WsAreaSqKm FROM StreamCat_2016 WHERE COMID in (%s)",inLOOP(substr(COMIDs,1,10))))


# run model applicability function
ModelApplicability = ModelApplicability(CalPredsModelApplicability,
                                        modelId = def_models$modelId,
                                        applicabilitypreds) # add to config file or add an R object with calpreds


