model_id=2 #change as needed
empty_samps<-210561

fake_bugs<-data.frame(sampleId=rep(empty_samps,2),
                      taxonomyId=c(135,59),
                      scientificName=c('Optioservus','Hygrobatidae'),
                      levelId=c(23,19),
                      levelName=c('Genus','Family'),
                      otuName=c('Optioservus','Acarina'),
                      splitCount=c(1,0))

sumrarefiedOTUTaxa = fake_bugs  %>%
  dplyr::group_by(sampleId, otuName) %>%
  dplyr::summarize(sumSplitCount = sum(splitCount)) # why are multiple records exported here per OTU???

sumrarefiedOTUTaxa$presence = ifelse(sumrarefiedOTUTaxa$sumSplitCount >=1, 1, 0)
bug_fake = tidyr::pivot_wider(sumrarefiedOTUTaxa,id_cols = "sampleId", names_from = "otuName",values_from = "presence")
bug_fake[is.na(bug_fake)]<-0
bug_fake=as.data.frame(bug_fake)
rownames(bug_fake)<-bug_fake$sampleId
bug_fake<-bug_fake[,-1]



fake_OE <-model.predict.RanFor.4.2(
  bugcal.pa,
  grps.final,
  preds.final,
  ranfor.mod,
  prednew,
  bugnew=bug_fake,
  Pc = 0.5,
  Cal.OOB = FALSE)


applicabilitypreds = NAMCr::query("samplePredictorValues",
                                  include = c(
                                    "sampleId",
                                    "predictorId",
                                    "status",
                                    "abbreviation",
                                    "predictorValue"),sampleIds=empty_samps) #need list of samples in database with values
applicabilitypreds = subset(applicabilitypreds, abbreviation %in% c('ElevCat','Tmean8110Ws','WsAreaSqKm','Precip8110Ws'))
applicabilitypreds$predictorValue=as.numeric(applicabilitypreds$predictorValue)
applicabilitypreds = tidyr::pivot_wider(applicabilitypreds,
                                        id_cols="sampleId",
                                        names_from = "abbreviation",
                                        values_from = "predictorValue")# add id_cols=sampleId once it gets added to end point
applicabilitypreds=as.data.frame(applicabilitypreds)


ModelApplicability_obj = ModelApplicability(CalPredsModelApplicability,
                                            modelId = model_id,
                                            applicabilitypreds)
Empty_OE<-merge(fake_OE,ModelApplicability_obj,by='row.names')

dat_to_pass<-list(sampleId = 210561,
                  modelId = 2,
                  oResult = Empty_OE$OE.scores.O,
                  eResult = Empty_OE$OE.scores.E,
                  modelResult = Empty_OE$OE.scores.OoverE ,
                  fixedCount = 0,
                  modelApplicability = ModelApplicability_obj$ModelApplicability,
                  #notes=General_OE_results[[i]]$InvasiveInvertSpecies[j])
                  notes='National')

NAMCr::save(
  api_endpoint = "setModelResult",
  args=dat_to_pass)




