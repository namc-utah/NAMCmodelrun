# ---------------------------------------------------------------
# Read in csv with results from COEDAS access database
# ---------------------------------------------------------------
#modelResults=read.csv("/Users/namc/Library/CloudStorage/Box-Box/NAMC/OE_Modeling/NAMC_Supported_OEmodels/CO/InputAndResults_CO2017MMI/Results1311.csv")
modelId=565
modelResults=read.csv("C://Users//andrew.caudillo.BUGLAB-I9//Box//NAMC//OEModeling//NAMC_Supported_OEmodels//Colorado//Documentation//2025 EDAS Distributable Package//exports//Proj6933.csv")
#modelResults=read.csv("C://Users//andrew.caudillo.BUGLAB-I9//Box//NAMC//OEModeling//NAMC_Supported_OEmodels//Colorado//Documentation//EDAS2017//working new version//exports//11845.csv")
modelResults = modelResults[,c("StationID","SiteClassification","MMI", "TotalInd")]
#modelResults=modelResults[modelResults$StationID %in% c(220496,220852),]
if(modelId %in% c(565,566,567)){
modelResults$BioType=ifelse(modelResults$SiteClassification=="1",565,
                            ifelse(modelResults$SiteClassification=="2",566,
                                   ifelse(modelResults$SiteClassification=="3",567,NA
                                   )))
}else(
  modelResults$BioType=ifelse(modelResults$SiteClassification=="1",4,
                              ifelse(modelResults$SiteClassification=="2",5,
                                     ifelse(modelResults$SiteClassification=="3",6,NA
                                     )))
)
modelResults=setNames(modelResults,c("sampleId","Biotype","MMI","Count","modelId"))
modelResults = modelResults %>%  na.omit()

# ---------------------------------------------------------------
# Always run model applicability test
# ---------------------------------------------------------------
# get all predictor values needed for a box or project # note this either needs a loop written over it or a different API end point
if(modelId %in% c(565,566,567)){

applicabilitypreds = NAMCr::query("samplePredictorValues",
                                  include = c(
                                    "sampleId",
                                    "predictorId",
                                    "status",
                                    "abbreviation",
                                    "predictorValue"
                                  ),
                                  sampleIds = modelResults$sampleId,
                                  modelIds=565
) #need list of samples in database with values
applicabilitypreds = subset(applicabilitypreds, abbreviation %in% c('ElevCat','Tmean8110Ws','WsAreaSqKm','Precip8110Ws'))
applicabilitypreds$predictorValue=as.numeric(applicabilitypreds$predictorValue)
applicabilitypreds = tidyr::pivot_wider(applicabilitypreds,
                                        id_cols="sampleId",
                                        names_from = "abbreviation",
                                        values_from = "predictorValue")# add id_cols=sampleId once it gets added to end point
applicabilitypreds=as.data.frame(applicabilitypreds)
# run model applicability function

ModelApplicability1 = ModelApplicability(CalPredsModelApplicability,
                                        modelId = 565,
                                        applicabilitypreds)
ModelApplicability1$modelId=565
ModelApplicability2 = ModelApplicability(CalPredsModelApplicability,
                                         modelId = 566,
                                         applicabilitypreds)
ModelApplicability2$modelId=566
ModelApplicability3 = ModelApplicability(CalPredsModelApplicability,
                                         modelId = 567,
                                         applicabilitypreds)
ModelApplicability3$modelId=567
}else{
  applicabilitypreds = NAMCr::query("samplePredictorValues",
                                    include = c(
                                      "sampleId",
                                      "predictorId",
                                      "status",
                                      "abbreviation",
                                      "predictorValue"
                                    ),
                                    sampleIds = modelResults$sampleId,
                                    modelIds=4
  ) #need list of samples in database with values
  applicabilitypreds = subset(applicabilitypreds, abbreviation %in% c('ElevCat','Tmean8110Ws','WsAreaSqKm','Precip8110Ws'))
  applicabilitypreds$predictorValue=as.numeric(applicabilitypreds$predictorValue)
  applicabilitypreds = tidyr::pivot_wider(applicabilitypreds,
                                          id_cols="sampleId",
                                          names_from = "abbreviation",
                                          values_from = "predictorValue")# add id_cols=sampleId once it gets added to end point
  applicabilitypreds=as.data.frame(applicabilitypreds)
  # run model applicability function

  ModelApplicability1 = ModelApplicability(CalPredsModelApplicability,
                                           modelId = 4,
                                           applicabilitypreds)
  ModelApplicability1$modelId=4
  ModelApplicability2 = ModelApplicability(CalPredsModelApplicability,
                                           modelId = 5,
                                           applicabilitypreds)
  ModelApplicability2$modelId=5
  ModelApplicability3 = ModelApplicability(CalPredsModelApplicability,
                                           modelId = 6,
                                           applicabilitypreds)
  ModelApplicability3$modelId=6
}
ModelAppy=rbind(ModelApplicability1,ModelApplicability2,ModelApplicability3)
finalResults=plyr::join(ModelAppy,modelResults,by='sampleId','left')
names(finalResults)[12]='MID'
finalResults2<-finalResults[!duplicated(finalResults$sampleId),]
#finalResults2=finalResults2[,names(finalResults2) != 'modelId.1']

# ---------------------------------------------------------------
# Get additional bug metrics (fixed count and invasives)
# ---------------------------------------------------------------
if(modelId %in% c(565,566,567)){
def_models = NAMCr::query(
  api_endpoint = "modelInfo",
  include = c("modelId",
              "modelTypeAbbreviation",
              "abbreviation",
              "translationId",
              "fixedCount"),
  modelId=565
)
}else{
  def_models = NAMCr::query(
    api_endpoint = "modelInfo",
    include = c("modelId",
                "modelTypeAbbreviation",
                "abbreviation",
                "translationId",
                "fixedCount"),
    modelId=4)
}
##### get fixed count column #####
bugsOTU = NAMCr::query("sampleTaxaTranslation",
                       translationId = def_models$translationId,
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
bugraw = subset(bugRaw,taxonomyId %in% c(1330,1331,2633, 2671,4933,4934,
                                         4935,4936,4937,4938,4939,4940,
                                         4941,4942,1019,1994,5096,
                                         1604,2000,4074,1369,2013,
                                         1579,9584,9585,9586,9587,601,632,611,
                                         650,936,1369,1514,1515,1518,
                                         1579,1604,1637,1658,1817,1953,
                                         1970,1971,1990,2000,2013,2053,8878,5510,
                                         607,820,3943,4194,3833,1030))
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

finalResults3=dplyr::left_join(finalResults2,additionalbugmetrics,by="sampleId")
finalResults3[finalResults3=='Orconectes']<-'Faxonius'
finalResults3$InvasiveInvertSpecies[is.na(finalResults3$InvasiveInvertSpecies)]<-'Absent'
#finalResults3$fixedCount=ifelse(finalResults3$Count>=300,300,finalResults3$Count)
b=NAMCr::query('samples',projectId=6933)
ress=NAMCr::query('modelResults',sampleIds=b$sampleId)
curr=ress[ress$modelId %in% c(565,566,567),]
for (i in 1:nrow(finalResults3) ){# need to add invasives and extra metrics to the notes field in some easy fashion???
  #has permission to save then spit out result to console
  # pass Nas for anything not used
  tryCatch({
    NAMCr::save(
      api_endpoint = "setModelResult",
      sampleId = finalResults3$sampleId[i],
      modelId = finalResults3$modelId[i],
      modelResult = finalResults3$MMI[i],
      fixedCount = finalResults3$Count[i],
      modelApplicability = finalResults3$ModelApplicability[i])
      notes=ifelse(finalResults3$InvasiveInvertSpecies[i]=='Absent','',finalResults3$InvasiveInvertSpecies[i])
    message(paste('saved result ',i,' of ',nrow(finalResults3)))

  }, error =function(e){
    cat(paste0("\n\tSAMPLE ERROR: ",finalResults3$sampleId[i],"\n"))
    str(e,indent.str = "   "); cat("\n")
  })
}
finalResults3

## Arizona
AZ_results<-data.frame(
  sampleId=
  c(218497,
    219183,
    219180,
    219185,
    219186,
    219182,
    218503,
    218504,
    219184,
    218500,
    218502,
    218507,
    218511,
    218505,
    218508,
    218512,
    218506,
    218510,
    218509,
    218501,
    219181,
    218498,
    218499),

IBI=
  c(32.61945687,
    30.11862601,
    20.7597006,
    33.26706496,
    49.84984515,
    13.32992151,
    48.2995001,
    38.63629338,
    53.60984936,
    36.31747273,
    29.94226164,
    48.27275826,
    44.26733863,
    26.08065393,
    64.87849214,
    40.14427199,
    66.13819974,
    39.93661668,
    33.41187011,
    42.4948816,
    16.70272591,
    22.31973582,
    15.17530466),

InvertReg=
  c('warm',
    'cold',
    'cold',
    'warm',
    'warm',
    'warm',
    'warm',
    'warm',
    'warm',
    'warm',
    'warm',
    'warm',
    'warm',
    'warm',
    'warm',
    'warm',
    'warm',
    'warm',
    'warm',
    'warm',
    'warm',
    'warm',
    'warm'),
fixedCount=
  c(500,
    500,
    500,
    71,
    500,
    500,
    500,
    499,
    372,
    500,
    500,
    498,
    500,
    500,
    500,
    500,
    500,
    500,
    500,
    499,
    175,
    172,
    500))

AZ_results$modelId=
  ifelse(AZ_results$InvertReg=='cold',236,169)


applicabilitypreds = NAMCr::query("samplePredictorValues",
                                  include = c(
                                    "sampleId",
                                    "predictorId",
                                    "status",
                                    "abbreviation",
                                    "predictorValue"
                                  ),
                                  sampleIds = AZ_results$sampleId,
                                  modelIds=169)
#need list of samples in database with values
applicabilitypreds = subset(applicabilitypreds, abbreviation %in% c('ElevCat','Tmean8110Ws','WsAreaSqKm','Precip8110Ws'))
applicabilitypreds$predictorValue=as.numeric(applicabilitypreds$predictorValue)
applicabilitypreds = tidyr::pivot_wider(applicabilitypreds,
                                        id_cols="sampleId",
                                        names_from = "abbreviation",
                                        values_from = "predictorValue")# add id_cols=sampleId once it gets added to end point
applicabilitypreds=as.data.frame(applicabilitypreds)
# run model applicability function
ModelApplicability = ModelApplicability(CalPredsModelApplicability,
                                        modelId = 169,
                                        applicabilitypreds) # add to config file or add an R object with calpreds

AZ_results=merge(AZ_results,ModelApplicability,by="sampleId")

for (i in 1:nrow(AZ_results) ){# need to add invasives and extra metrics to the notes field in some easy fashion???
  #has permission to save then spit out result to console
  # pass Nas for anything not used
  tryCatch({
    NAMCr::save(
      api_endpoint = "setModelResult",
      sampleId = AZ_results$sampleId[i],
      modelId = AZ_results$modelId[i],
      modelResult = AZ_results$IBI[i],
      fixedCount = AZ_results$fixedCount[i],
      modelApplicability = AZ_results$ModelApplicability[i])
    message(paste('saved result ',i,' of ',nrow(AZ_results)))

  }, error =function(e){
    cat(paste0("\n\tSAMPLE ERROR: ",AZ_results$sampleId[i],"\n"))
    str(e,indent.str = "   "); cat("\n")
  })
}

