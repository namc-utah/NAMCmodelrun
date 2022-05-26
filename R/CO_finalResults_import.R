# ---------------------------------------------------------------
# Read in csv with results from COEDAS access database
# ---------------------------------------------------------------
modelResults=read.csv("/Users/triparmstrong/Library/CloudStorage/Box-Box/NAMC/OE_Modeling/NAMC_Supported_OEmodels/CO/InputAndResults_CO2017MMI/AIM 2021/May19results.csv")
modelResults = modelResults[,c("StationID","SiteClassification","MMI", "TotalInd")] # what about TotalInd_ why is this over 300 in some cases... should this be the fixedcount we use instead of down below?
modelResults$BioType=ifelse(modelResults$SiteClassification=="1",4,
                            ifelse(modelResults$SiteClassification=="2",5,
                                   ifelse(modelResults$SiteClassification=="3",6,NA
                                   )))
modelResults=setNames(modelResults,c("sampleId","Biotype","MMI","Count","modelId"))


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
                                  sampleIds = modelResults$sampleId
) #need list of samples in database with values
applicabilitypreds = subset(applicabilitypreds, abbreviation %in% c('ElevCat','Tmean8110Ws','WsAreaSqKm','Precip8110Ws'))
applicabilitypreds = tidyr::pivot_wider(applicabilitypreds,
                                        id_cols="sampleId",
                                        names_from = "abbreviation",
                                        values_from = "predictorValue")# add id_cols=sampleId once it gets added to end point
applicabilitypreds=as.data.frame(applicabilitypreds)
# run model applicability function
ModelApplicability = ModelApplicability(CalPredsModelApplicability,
                                        modelId = 4,
                                        applicabilitypreds) # add to config file or add an R object with calpreds

finalResults=merge(modelResults,ModelApplicability,by="row.names")

# ---------------------------------------------------------------
# Get additional bug metrics (fixed count and invasives)
# ---------------------------------------------------------------
##### get fixed count column #####
bugsOTU = NAMCr::query("sampleTaxaTranslationRarefied",
                       translationId = def_models$translationId,
                       fixedCount = def_models$fixedCount,
                       sampleIds=modelResults$sampleId
)
sumrarefiedOTUTaxa = bugsOTU  %>%
  dplyr::group_by(sampleId) %>%
  dplyr::summarize(fixedCount = sum(splitCount))


###### get invasives #####
# get raw bug data
bugRaw = NAMCr::query(
  "sampleTaxa",
  sampleIds=modelResults$sampleId
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

for (i in 1:nrow(finalResults) ){# need to add invasives and extra metrics to the notes field in some easy fashion???
  #has permission to save then spit out result to console
  # pass Nas for anything not used
  tryCatch({
    NAMCr::save(
      api_endpoint = "setModelResult",
      sampleId = finalResults$sampleId[i],
      modelId = modelResults$modelId[1],
      modelResult = finalResults$MMI[i],
      fixedCount = finalResults$fixedCount[i],
      modelApplicability = finalResults$ModelApplicability[i],
      notes=finalResults$InvasiveInvertSpecies[i]
    )

  }, error =function(e){
    cat(paste0("\n\tSAMPLE ERROR: ",finalResults$sampleId[i],"\n"))
    str(e,indent.str = "   "); cat("\n")
  })
}
