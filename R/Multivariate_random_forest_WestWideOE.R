
  #Query to pull down predictors and taxa
  #but will just use the raw P/A matrix that we read into the
  #file for taxa
  query<- ("SELECT
    r.sample_id AS sampleId,
    r.translation_id AS translationId,
    z.taxonomy_id AS taxonomyId,
    t.scientific_name AS scientificName,
    z.split_count,
    p.abbreviation,
    sp.predictor_value,
    (
    SELECT DISTINCT COALESCE((
        WITH RECURSIVE TaxonomyHierarchy AS (/* create the CTE, essentially a for loop*/
            SELECT taxonomy_id, scientific_name, parent_id /* this is the base query of the CTE. Checking for matches {*/
            FROM taxonomy
            WHERE taxonomy_id = z.taxonomy_id  /* }*/
            UNION ALL /* combines recursive and base queries*/
            SELECT t.taxonomy_id, t.scientific_name, t.parent_id /* this is the looping part {*/
            FROM taxonomy t
            JOIN TaxonomyHierarchy th ON t.taxonomy_id = th.parent_id
        ) /* }*/ /*leave the recursive section*/
        SELECT g.scientific_name /* now let's take the scientific name*/
        FROM TaxonomyHierarchy th /* from the big nasty CTE we just made*/
        LEFT JOIN translation_taxa f ON th.taxonomy_id = f.taxonomy_id /* join f and th, by taxId, to start getting the OTU name*/
        LEFT JOIN translation_groups g ON f.translation_group_id = g.translation_group_id /*join g and f. Now we can get the OTU name*/
        WHERE g.translation_id = 21 /* filter*/
        LIMIT 1 /* make sure only one row is returned*/
    ), NULL)
) As OTU, /*name the column*/
  /*the rest is all pretty standard SQL and joins*/
  z.split_count AS splitCount
FROM
rarefactions r
JOIN
rarefied_organisms z ON r.rarefaction_id = z.rarefaction_id
JOIN
taxonomy t ON z.taxonomy_id = t.taxonomy_id
join samples s on s.sample_id = r.sample_id
join model_samples ms on ms.sample_id = r.sample_id
join site_predictors sp on sp.site_id = s.site_id
join predictors p on sp.predictor_id = p.predictor_id
join model_predictors mp on mp.predictor_id = sp.predictor_id
WHERE
r.sample_id in (select sample_id from samples where sample_id
in (select sample_id from model_samples where model_id = 25 ))
AND r.translation_id = 21
and mp.model_id = 25;")
  #connect to the database
  mydb=DBI::dbConnect(RSQLite::SQLite(),'C://NAMC_S3//LegacyDatabases//instar.sqlite')
  #now execute the query and save as a data frame
  ref_dat<-DBI::dbGetQuery(mydb,query)
  library(tidyr)
  library(dplyr)
  #read in ref data
  ben_dat=read.csv("C://Users//andrew.caudillo.BUGLAB-I9//Box//NAMC//Research Projects//AIM//IncreaserDecreaser_OE//Observed_PA.csv")
  #reference data has an extra Tanypodinae for some reason due to the model run code
  #just combine them together.
  ben_dat<-ben_dat%>%
    mutate(Tanypodinae = pmax(Tanypodinae, Tanypodinae.1, na.rm = TRUE)) %>%
    select(-Tanypodinae.1)
  refs<-read.csv('C://Users//andrew.caudillo.BUGLAB-I9//Box//NAMC//Research Projects//AIM//IncreaserDecreaser_OE//WestWide_development_samples.csv')
  #remove 3 sites for which predictors cannot be calculated
  ben_dat<-ben_dat[ben_dat$X %in% c(116830,132243,132807, 185022,116435,116436,117875,116849)==F,] #last 3 are new
  names(ben_dat)[1]<-'sampleId'


  # Wide-format predictor values — ensures only one column per predictor_id
  predictor_wide <- ref_dat %>%
     select(sampleId, abbreviation, predictor_value) %>%
     distinct() %>%
     pivot_wider(names_from = abbreviation, values_from = predictor_value)

 #join the two datasets together
  ben_dat=plyr::join(ben_dat,predictor_wide)
#omit Carabidae
  ben_dat<-ben_dat[,names(ben_dat) %in% 'Carabidae'==F]
  #ben_dat<-ben_dat[,names(ben_dat) %in% taxa_notraits_rare$taxon==F]


#load in MRF
  library(randomForestSRC)  # for multivariate random forest

  # ---------------------------
  # set the random seed and isolate predictor names
  set.seed(42)
  pred_index<-which(names(ben_dat)==names(predictor_wide)[2])

 #force the predictors to be numbers, not strings, as SQL shows
  ben_dat[,pred_index:ncol(ben_dat)]<-lapply(ben_dat[,pred_index:ncol(ben_dat)],as.numeric)

  #just renaming ben_dat for a cleanliness here forward
  train <- ben_dat
  #assign the sampleId to row names, then drop it
  row.names(train)<-train$sampleId;train=train[,-1]
  #the names of the taxa (responses)
   taxa_names<-names(train[,1:(pred_index-2)])
   #just_tax=train[,taxa_names]
   #tax_filter=just_tax[, apply(just_tax, 2, function(x) length(unique(x)) < 2)]

#variables that the model finds as the top variables via the
#"trees" function.
topvars=c("CaOWs","NCat","PctSalLake","Precip8110","MAST_mean08091314")
  formula_string <- paste0(
    "cbind(", paste(names(train[,1:(pred_index-2)]), collapse = ", "), ") ~ ",
    paste(names(train[,topvars]),collapse = '+'))
    #paste(names(train[,(ncol(train)-13):ncol(train)]),collapse = '+'))
  # Convert string to formula using parse + eval
  #Else MRF will not read it correctly
  rf_formula <- eval(parse(text = formula_string))
  cases=1 / (rowSums(train[, taxa_names]) + 1)
  case_wts=cases / sum(cases)

  #case weights
  if(0){compute_hybrid_weights <- function(train, resp_cols, power = 0.5, max_scale = 5) {

    # Step 1: count presences per site
    presences <- rowSums(train[, resp_cols])

    # Step 2: invert to favor low-presence sites
    raw_w <- 1 / (presences + 1)

    # Step 3: apply a smoothing power to reduce extreme weights
    smooth_w <- raw_w^power

    # Step 4: optionally cap extreme values
    smooth_w[smooth_w > max_scale] <- max_scale

    # Step 5: normalize mean = 1
    smooth_w <- smooth_w / mean(smooth_w)

    return(smooth_w)
  }}


  #case_weights <- compute_hybrid_weights(train, taxa_names, power = 0.5, max_scale = 3)


  #this model has already been fit with the best hyperparameters
  #as found by the tune() function
  rf_model <- rfsrc(
    formula=rf_formula,
    data=train,
    ntree = 500,
    mtry=4, #changed from 4 on Oct 13
    importance = 'permute',
    nodedepth = 5, #changed from 5
    nodesize=20,
    nsplit=1,
    case.wt = case_wts#changed from 10
    ) #added on 13 Oct
#view OOB errors for each response
OOB_pred_each=get.mv.error(rf_model)
mean(OOB_pred_each)
#pseuod RMSE
mean(1- colSums((rf_model$yvar - OOB_pred_each)^2, na.rm = TRUE) /
  colSums((rf_model$yvar - colMeans(rf_model$yvar, na.rm = TRUE))^2, na.rm = TRUE))
  #get the predicted Pcs for each taxon here
  results<-data.frame(sapply(rf_model$regrOutput, function(x) x$predicted.oob))
  results_inbag=data.frame(sapply(rf_model$regrOutput, function(x) x$predicted))
  #set row names
  row.names(results)=row.names(train)
  #adjust hyperparamters for best model
  #if set to if(0), it has already been run!
  if(0){
  tune(rf_formula,
             train,
             doBest = T,
             improve=1e-1, #define what improvement means for OOB error
       sampsize = 500, #number of samples for the subsamples
       mtryStart =round(length(pred_index:ncol(train))/3), #default from rfsrc doc
       trace=T, #show the progress in console
       nodesizeTry = c(5, 10,20,30),
       ntreeTry = 500,
       nsplit=5) #a small number of trees for memory usage
  }


quick_results=cbind.data.frame(Pred=t(round(results[1,],2)),
             Obs=t(train[1,names(train) %in% names(results)]))
names(quick_results)=c('Pred','Obs')

  #calculating Fo and Fe, then subsetting out those which are uncommon
  #and/or do not have traits
  Fos=colSums(train[,1:(pred_index-2)])
  Fes=colSums(results)
  Fos_small=Fos[names(Fos) %in% taxa_notraits_rare$taxon==F]
  Fes_small=Fes[names(Fes) %in% taxa_notraits_rare$taxon==F]
  ratio=Fos/Fes
  ratio_dat=data.frame(Fo=Fos_small,Fe=Fes_small)
  ratio_small=ratio[names(ratio) %in% taxa_notraits_rare$taxon==F]
  boxplot(ratio_small,ylab='Fo/Fe ratio')
  savp(10,8,'C://Users//andrew.caudillo.BUGLAB-I9//Box//NAMC//Research Projects//AIM//IncreaserDecreaser_OE//MRF_OE//Fo_Fe_boxplot_251023.png')
  plot(Fes,Fos,ylab='Fo',xlab='Fe',main='Reference Fo vs Fe')
  abline(0,1,col='red')
  savp(10,8,'C://Users//andrew.caudillo.BUGLAB-I9//Box//NAMC//Research Projects//AIM//IncreaserDecreaser_OE//MRF_OE//Ref_FovsFe_251023.png')


  # Predict on prob. sites


  #read in prob. sites and their predictors

  query2<- ("SELECT
    r.sample_id AS sampleId,
    r.translation_id AS translationId,
    z.taxonomy_id AS taxonomyId,
    t.scientific_name AS scientificName,
    z.split_count,
    p.abbreviation,
    sp.predictor_value,
    (
    SELECT DISTINCT COALESCE((
        WITH RECURSIVE TaxonomyHierarchy AS (/* create the CTE, essentially a for loop*/
            SELECT taxonomy_id, scientific_name, parent_id /* this is the base query of the CTE. Checking for matches {*/
            FROM taxonomy
            WHERE taxonomy_id = z.taxonomy_id  /* }*/
            UNION ALL /* combines recursive and base queries*/
            SELECT t.taxonomy_id, t.scientific_name, t.parent_id /* this is the looping part {*/
            FROM taxonomy t
            JOIN TaxonomyHierarchy th ON t.taxonomy_id = th.parent_id
        ) /* }*/ /*leave the recursive section*/
        SELECT g.scientific_name /* now let's take the scientific name*/
        FROM TaxonomyHierarchy th /* from the big nasty CTE we just made*/
        LEFT JOIN translation_taxa f ON th.taxonomy_id = f.taxonomy_id /* join f and th, by taxId, to start getting the OTU name*/
        LEFT JOIN translation_groups g ON f.translation_group_id = g.translation_group_id /*join g and f. Now we can get the OTU name*/
        WHERE g.translation_id = 21 /* filter*/
        LIMIT 1 /* make sure only one row is returned*/
    ), NULL)
) As OTU, /*name the column*/
  /*the rest is all pretty standard SQL and joins*/
  z.split_count AS splitCount
FROM
rarefactions r
JOIN
rarefied_organisms z ON r.rarefaction_id = z.rarefaction_id
JOIN
taxonomy t ON z.taxonomy_id = t.taxonomy_id
join samples s on s.sample_id = r.sample_id
join project_samples ps on ps.sample_id = r.sample_id
join site_predictors sp on sp.site_id = s.site_id
join predictors p on sp.predictor_id = p.predictor_id
join model_predictors mp on mp.predictor_id = sp.predictor_id
WHERE
r.sample_id in (select sample_id from samples where sample_id
in (select sample_id from project_samples where project_id = 5964 ))
AND r.translation_id = 21
and mp.model_id = 25")

  #this is identical to the reference steps, except
  #we only give MRF the predictors
  test_pred_dat<-DBI::dbGetQuery(mydb,query2)
  test_pred_dat<-test_pred_dat[!duplicated(test_pred_dat),]

  test_dat=read.csv('C://Users//andrew.caudillo.BUGLAB-I9//Box//NAMC//Research Projects//AIM//IncreaserDecreaser_OE//Prob_O.csv')
  test_dat$Tanypodinae = test_dat$Tanypodinae + test_dat$Tanypodinae.1;test_dat<-test_dat[,names(test_dat) %in% 'Tanypodinae.1'==F,]
  names(test_dat)[1]<-'sampleId'
  failed_sites<-read.csv('C://Users//andrew.caudillo.BUGLAB-I9//Box//NAMC//Research Projects//AIM//IncreaserDecreaser_OE//failed_sites.csv')
  #test_dat<-test_dat[test_dat$sampleId %in% failed_sites$sampleId==F,]


  predictor_test_wide <- test_pred_dat %>%
    select(sampleId, abbreviation, predictor_value) %>%
    distinct() %>%
    pivot_wider(names_from = abbreviation, values_from = predictor_value)

  # 3. Merge: join both by sample_id (one row per sample)
  test_dat<-test_dat[test_dat$sampleId %in% failed_sites$sampleId==F,]
  test_dat=plyr::join(test_dat,predictor_test_wide)

  test_dat<-test_dat[,names(test_dat) %in% 'Carabidae'==F]
  pred_test_index<-which(names(test_dat)==names(predictor_test_wide)[2])
  rownames(test_dat) <- test_dat$sampleId;test_dat<-test_dat[,names(test_dat) %in% 'sampleId'==F]

  new_sites <- test_dat[,(pred_test_index-1):ncol(test_dat)]
  row.names(new_sites)=row.names(test_dat)
  test_site_names<-row.names(new_sites)
  new_sites<-as.data.frame(lapply(new_sites,as.numeric))
  row.names(new_sites)<-test_site_names
#predict the model onto the prob. sites
  pred <- predict(rf_model, newdata = new_sites)


  # Extract predicted probabilities for each taxon at prob sites
  pred_probs <- as.data.frame(sapply(pred$regrOutput, function(x) x$predicted))
  row.names(pred_probs)<-row.names(new_sites)
  #write.csv(pred_probs,'C://Users//andrew.caudillo.BUGLAB-I9//Box//NAMC//Research Projects//AIM//IncreaserDecreaser_OE//MRF_Prob_Pcs.csv')
#calculating Fo and Fe
  Pr_PA=test_dat[,names(test_dat) %in% names(pred_probs)]
  Pr_PA=as.data.frame(ifelse(Pr_PA >0,1,0))
  #write.csv(Pr_PA, 'C://Users//andrew.caudillo.BUGLAB-I9//Box//NAMC//Research Projects//AIM//IncreaserDecreaser_OE//MRF_Prob_Os.csv')
  PFos=colSums(Pr_PA)
  PFes=colSums(pred_probs)

  Pr_PA_oth=Pr_PA[row.names(Pr_PA) %in% row.names(Pesgrp2)==F,]
  pred_probs_oth=pred_probs[row.names(pred_probs) %in% row.names(Pesgrp2)==F,]

  Pr_PA_x=Pr_PA[row.names(Pr_PA) %in% row.names(Pesgrp2),]
  pred_probs_x=pred_probs[row.names(pred_probs) %in% row.names(Pesgrp2),]

  Pratio=PFos/PFes
  Pratio_small=Pratio[names(Pratio) %in% taxa_notraits_rare$taxon==F]
  PFos_small=PFos[names(PFos) %in% taxa_notraits_rare$taxon==F]
  PFes_small=PFes[names(PFes) %in% taxa_notraits_rare$taxon==F]
  Poth_Fo=colSums(Pr_PA_oth);Poth_Fe=colSums(pred_probs_oth)
  Poth_Fosmall=Poth_Fo[names(Poth_Fo) %in% taxa_notraits_rare$taxon==F]
  Poth_Fesmall=Poth_Fe[names(Poth_Fe) %in% taxa_notraits_rare$taxon==F]
  Px_Fo=colSums(Pr_PA_x);Px_Fe=colSums(pred_probs_x)
  Px_Fosmall=Px_Fo[names(Px_Fo) %in% taxa_notraits_rare$taxon==F]
  Px_Fesmall=Px_Fe[names(Px_Fe) %in% taxa_notraits_rare$taxon==F]
  #Poth_ratio=Poth_Fo/Poth_Fe
  #Px_ratio=Px_Fo/Px_Fe
  Poth_ratio_small=Poth_Fosmall/Poth_Fesmall
  Px_ratio_small=Px_Fosmall/Poth_Fesmall
  plot(PFes_small,PFos_small,ylab='Fo',xlab='Fe',main='Prob. Fo vs Fe')
  abline(0,1,col='red')

  Poth_plotdat=data.frame(Fo=Poth_Fosmall,Fe=Poth_Fesmall)
  Px_plotdat=data.frame(Fo=Px_Fosmall,Fe=Px_Fesmall)
  boxplot(Pratio_small,ylab='Fo/Fe ratio')
  savp(10,8,'C://Users//andrew.caudillo.BUGLAB-I9//Box//NAMC//Research Projects//AIM//IncreaserDecreaser_OE//MRF_OE//Prob_FoFe_ratio_251023.png')
  savp(10,8,'C://Users//andrew.caudillo.BUGLAB-I9//Box//NAMC//Research Projects//AIM//IncreaserDecreaser_OE//MRF_OE//Prob_FovsFe_box_251023.png')


  boxplot(Pratio_small,at=2,xlim=c(1,4),ylab='Fo/Fe')
  boxplot(ratio_small,at=3,xlim=c(1,4),add=T)

  axis(side=1,at=c(2,3),labels = c('Prob','Train'))


  savp(10,8,'C://Users//andrew.caudillo.BUGLAB-I9//Box//NAMC//Research Projects//AIM//IncreaserDecreaser_OE//MRF_OE//Box_compare_FoFes_251023.png')

boxplot(ratio_small,at=1,xlim=c(0,3),ylim=c(0,max(na.omit(Pratio_small))),col='purple3',ylab='Fo/Fe ratio')
boxplot(Pratio_small,at=2,add=T,col='yellow3')
legend('topright',
       leg=c('Reference',
             'Probabilistic'),
       pch=rep(22,2),
       pt.bg=c('purple4',
               'yellow3'),
       bty='n',
       cex=0.8)
savp(10,8,'C://Users//andrew.caudillo.BUGLAB-I9//Box//NAMC//Research Projects//AIM//IncreaserDecreaser_OE//MRF_OE//regional_boxplot_compare.png')

boxplot(Oth_ratio_small,at=1,xlim=c(0,5),ylim=c(0,max(na.omit(Poth_ratio_small))),col='purple3',ylab='Fo/Fe ratio')
boxplot(Poth_ratio_small,at=2,add=T,col='yellow3')
boxplot(Xer_ratio_small,at=3,add=T,col='purple3')
boxplot(Px_ratio_small,at=4,add=T,col='yellow3')
mtext(text = c("Other", "Xeric"), side = 1, line = 1, at = c(1.5, 3.5), cex = 1)
legend('topright',
       leg=c('Reference',
             'Probabilistic'),
       pch=rep(22,2),
       pt.bg=c('purple4',
               'yellow3'),
       bty='n',
       cex=0.8)
savp(10,8,'C://Users//andrew.caudillo.BUGLAB-I9//Box//NAMC//Research Projects//AIM//IncreaserDecreaser_OE//MRF_OE//ecoregion_boxplot_compare.png')
# ---------------------------
#creating O/E scores, quasi-RIVPACS style
#set the Pc threshold, in this case, 0.5
threshold <- 0.5
  #define fxn that will calculate the O/E scores
OE_calc=function(results_data, PA,threshold){
  Os=list()
  Es=list()
  BC=list()
  for(i in 1:nrow(results_data)) {

    cur.prd<-results_data[i,]; #vector of predicted probs for current sample;
    spdyn<-names(cur.prd)[cur.prd>=threshold];  #subset of taxa with Pi>=Pcutoff;
    cur.prd<-cur.prd[spdyn]
    cur.obs<-PA[i,spdyn]; #vector of observed P/A for those taxa;
    Os[[i]]<-sum(cur.obs); #observed richness (O);
    Es[[i]]<-sum(cur.prd); #Expected richness (E);
    BC[[i]]<-sum(abs(cur.obs-cur.prd))/ (sum(cur.obs)+sum(cur.prd)); #BC value;
  }
  OE.dat=as.data.frame(cbind(unlist(Os),unlist(Es),unlist(BC)))

  names(OE.dat)<-c('O','E','BC')
  OE.dat$OtoE=OE.dat$O/OE.dat$E
  #OE.dat$E.adj = OE.dat$E*(mean(OE.dat$OtoE))
  #OE.dat$OtoE.adj=OE.dat$O / OE.dat$E.adj
  return(OE.dat)
}
#run it on ref and prob. sites
EastXer=results[row.names(results) %in% groups$`Eastern Xeric Plateaus`$sampleId,]
OthEco=results[row.names(results) %in% groups$`Eastern Xeric Plateaus`$sampleId==F,]
XerTrain=train[row.names(train) %in% groups$`Eastern Xeric Plateaus`$sampleId,]
XerTrain=XerTrain[,1:(ncol(XerTrain)-15)]
OthTrain=train[row.names(train) %in% groups$`Eastern Xeric Plateaus`$sampleId==F,]
OthTrain=OthTrain[,1:(ncol(OthTrain)-15)]
XerFo=colSums(XerTrain)
XerFe=colSums(EastXer)
Xer_ratio=XerFo/XerFe

OthFo=colSums(OthTrain)
OthFe=colSums(OthEco)
Oth_ratio=OthFo/OthFe
Oth_ratio_small=Oth_ratio[names(Oth_ratio) %in% names(taxa_notraits_rare)==F]
Xer_ratio_small=Xer_ratio[names(Xer_ratio) %in% names(taxa_notraits_rare)==F]


par(mfrow=c(1,2))
plot(OthFe,OthFo)
abline(0,1,col='red')
plot(XerFe,XerFo)
abline(0,1,col='red')

OthEco_plotdat=data.frame(Fo=OthFo,Fe=OthFe)
EastXer_plotdat=data.frame(Fo=XerFo,Fe=XerFe)

# WW_Xer=Os[Os$X %in% groups$`Eastern Xeric Plateaus`$sampleId==T,]
# WW_E_Xer=Es[Es$X %in% groups$`Eastern Xeric Plateaus`$sampleId==T,]
# WW_XerFo=colSums(WW_Xer[,2:ncol(WW_Xer)])
# WW_XerFe=colSums(WW_E_Xer[,2:ncol(WW_E_Xer)])
# WW_XerFe_ratio=WW_XerFo/WW_XerFe
# WW_XerFe_plotdat=data.frame(Fo=WW_XerFo,Fe=WW_XerFe)
#
# WW=Os[Os$X %in% groups$`Eastern Xeric Plateaus`$sampleId==F,]
# WW_E=Es[Es$X %in% groups$`Eastern Xeric Plateaus`$sampleId==F,]
# WW_Fo=colSums(WW[,2:ncol(WW)])
# WW_Fe=colSums(WW_E[,2:ncol(WW_E)])
# WW_Fe_ratio=WW_Fo/WW_Fe
# WW_plotdat=data.frame(Fo=WW_Fo,Fe=WW_Fe)

combined_dat=data.frame(Fo=Fos_small,Fe=Fes_small)
A<-ggplot(data=OthEco_plotdat,aes(y=Fo,x=Fe))+geom_point()+geom_abline(intercept = 0,slope = 1,col='red')+ggtitle('Other Ecoregions')
B<-ggplot(data=EastXer_plotdat,aes(y=Fo,x=Fe))+geom_point()+geom_abline(intercept = 0,slope = 1,col='red')+ggtitle('Eastern Xeric')
ggplot(data=combined_dat,aes(x=Fe,y=Fo))+geom_point()+geom_abline(intercept = 0,slope = 1,col='red')
savp(10,8,'C://Users//andrew.caudillo.BUGLAB-I9//Box//NAMC//Research Projects//AIM//IncreaserDecreaser_OE//MRF_OE//MRF_allsites_FoFe.png')
C<-ggplot(data=Poth_plotdat,aes(y=Fo,x=Fe))+geom_point()+geom_abline(intercept = 0,slope = 1,col='red')
D<-ggplot(data=Px_plotdat,aes(y=Fo,x=Fe))+geom_point()+geom_abline(intercept = 0,slope = 1,col='red')

gridExtra::grid.arrange(A,B,C,D,ncol=2)
savp(10,8,'C://Users//andrew.caudillo.BUGLAB-I9//Box//NAMC//Research Projects//AIM//IncreaserDecreaser_OE//MRF_OE//MRF_ecoregions_FoFe_sitecompare.png')
MRF_ref_OEs=OE_calc(results_data = results,
                    PA=train,
                    threshold=threshold)
row.names(MRF_ref_OEs)=row.names(results)
MRF_test_OEs=OE_calc(results_data=pred_probs,
                     PA=Pr_PA,
                     threshold=threshold)
XerOE=OE_calc(results_data = EastXer,
              PA=XerTrain,
              threshold=threshold)
OthEco_OE=OE_calc(results_data = OthEco,
              PA=OthTrain,
              threshold=threshold)

POthOE=OE_calc(results_data = pred_probs_oth,
               PA=Pr_PA_oth,
               threshold = threshold)
PXOE=OE_calc(results_data = pred_probs_x,
             PA=Pr_PA_x,
             threshold = threshold)
# WW_calc=OE_calc(results_data = WW_E[,-1],
#                 PA=WW[,-1],
#                 threshold = threshold)
# WW_X_calc=OE_calc(results_data = WW_E_Xer[,-1],
#                 PA=WW_Xer[,-1],
#                 threshold = threshold)
sd(OthEco_OE$OtoE)
boxplot(OthEco_OE$OtoE,at=1,xlim=c(0,5),ylab='O/E score',col='purple4',ylim=c(0,1.4),
        names=c('A'))
boxplot(POthOE$OtoE,at=2,add=T,col='yellow3')
boxplot(XerOE$OtoE,at=3,add=T,col='purple4')
boxplot(PXOE$OtoE,at=4,add=T,col='yellow3')
#axis(1,at=c(1:4),labels = c('WW','MRF WW','WW EXP','MRF EXP'),cex.axis=0.8 )
mtext(text = c("Other", "Xeric"), side = 1, line = 1, at = c(1.5, 3.5), cex = 1)
legend('topright',
       leg=c('Reference',
             'Probabilistic'),
       pch=rep(22,2),
       pt.bg=c('purple4',
               'yellow3'),
       bty='n',
       cex=0.8)
#
savp(10,8,'C://Users//andrew.caudillo.BUGLAB-I9//Box//NAMC//Research Projects//AIM//IncreaserDecreaser_OE//MRF_OE//OtoE_boxes_ecoregions_251107.png')
#savp(10,8,'C://Users//andrew.caudillo.BUGLAB-I9//Box//NAMC//Research Projects//AIM//IncreaserDecreaser_OE//MRF_OE//OtoE_boxes_Pc05_251106.png')
### This is just looking at O/E performance and metrics surrounding it
hist(MRF_ref_OEs$OtoE)
hist(MRF_test_OEs$OtoE)

sd(MRF_ref_OEs$OtoE)
sd(MRF_test_OEs$OtoE)

rmse_train <- sqrt(mean((MRF_ref_OEs$OtoE - 1)^2)); rmse_test <- sqrt(mean((MRF_test_OEs$OtoE - 1)^2))
rmse_xer=sqrt(mean((XerOE$OtoE - 1)^2))
rmse_oth=sqrt(mean((OthEco_OE$OtoE - 1)^2))
rmse_train
rmse_test

quantile(MRF_ref_OEs$OtoE, probs = c(0.01,0.05,0.1,0.25,0.5,0.75,0.9,0.95,0.99))

boxplot(MRF_ref_OEs$OtoE,at=1,xlim=c(0,3),ylab='O/E',ylim=c(min(MRF_test_OEs$OtoE),max(MRF_ref_OEs$OtoE)),col='purple4')
boxplot(MRF_test_OEs$OtoE,at=2,add=T,col='yellow3')
legend('topright',
       leg=c('Reference',
             'Probabilistic'),
       pch=rep(22,2),
       pt.bg=c('purple4',
               'yellow3'),
       bty='n',
       cex=0.8)

#axis(side=1,at=c(1,2),labels=c('Ref','Prob.'))
abline(h=median(MRF_ref_OEs$OtoE),lty=2)
abline(h=median(MRF_test_OEs$OtoE),lty=2)

savp(10,8,'C://Users//andrew.caudillo.BUGLAB-I9//Box//NAMC//Research Projects//AIM//IncreaserDecreaser_OE//MRF_OE//OtoE_boxes_Pc05_251107.png')



MRF_WW_ratios=data.frame(Median_ratio=median(na.omit(ratio_small)),
           Median_ratio_p=median(na.omit(Pratio_small)),
           mean_ratio=mean(na.omit(ratio_small)),
           mean_ratio_p=mean(na.omit(Pratio_small))
)
MRF_WW_ratios

graphics.off()
plot(MRF_ref_OEs$E,MRF_ref_OEs$O,ylab='O',xlab='E',main='Ref. O vs E',ylim=c(5,max(MRF_ref_OEs$O)),
     xlim=c(5,max(MRF_ref_OEs$O)))
abline(0,1,col='red')
OE_lm<-lm(MRF_ref_OEs$O~MRF_ref_OEs$E)
abline(OE_lm,col='blue')
legend('topleft',
       leg=c('1:1',
             'Regression'),
       lty=rep(1,2),
       col=c('red','blue'),
       bty='n')
savp(10,8,'C://Users//andrew.caudillo.BUGLAB-I9//Box//NAMC//Research Projects//AIM//IncreaserDecreaser_OE//MRF_OE//OvsE_Pc0_251023.png')

coeffs=coef(OE_lm)
coeffs
MRF_WW_OEs=data.frame(Median_OE=median(na.omit(MRF_ref_OEs$OtoE)),
                      Median_OE_p=median(na.omit(MRF_test_OEs$OtoE)),
                      mean_OE=mean(na.omit(MRF_ref_OEs$OtoE)),
                      mean_OE_p=mean(na.omit(MRF_test_OEs$OtoE)),
                      SD=sd(na.omit(MRF_ref_OEs$OtoE)),
                      SD_p=sd(na.omit(MRF_test_OEs$OtoE)),
                      Slope=coeffs[2],
                      R2=summary(OE_lm)$adj.r.squared)
MRF_WW_OEs
summary(OE_lm)$adj.r.squared



#clipr::write_clip(MRF_WW_OEs)


boxplot(get.mv.error(rf_model),ylab='taxon OOB Error')
savp(10,8,'C://Users//andrew.caudillo.BUGLAB-I9//Box//NAMC//Research Projects//AIM//IncreaserDecreaser_OE//MRF_OE//MRF_OOB_per_response_box.png')


#testing model on known low OE score sites.

deg=read.csv('C://Users//andrew.caudillo.BUGLAB-I9//Box//NAMC//Research Projects//AIM//IncreaserDecreaser_OE//MRF_OE//low_OE_sites.csv')
predictor_deg_wide <- as.data.frame(deg %>%
  select(sample_id, abbreviation, predictor_value) %>%
  distinct() %>%
  pivot_wider(names_from = abbreviation, values_from = predictor_value))

row.names(predictor_deg_wide)=predictor_deg_wide$sample_id
predictor_deg_wide<-predictor_deg_wide[,names(predictor_deg_wide) %in% 'sample_id'==F]
predictor_deg_wide_sub=predictor_deg_wide[sample(nrow(predictor_deg_wide)*0.1),]
deg_samps=row.names(predictor_deg_wide_sub)
predictor_deg_wide_sub<-as.data.frame(lapply(predictor_deg_wide_sub,as.numeric))
row.names(predictor_deg_wide_sub)=deg_samps
deg_pred=predict(rf_model,predictor_deg_wide_sub)
deg_results<-as.data.frame(sapply(deg_pred$regrOutput, function(x) x$predicted))
row.names(deg_results)=row.names(predictor_deg_wide_sub)
deg_sites=read.csv('C://Users//andrew.caudillo.BUGLAB-I9//Box//NAMC//Research Projects//AIM//IncreaserDecreaser_OE//MRF_OE//Low_score_sites.csv')
row.names(deg_sites)<-deg_sites$X;deg_sites=deg_sites[,-1]
deg_sites<-deg_sites[row.names(deg_sites) %in% row.names(deg_results),]

deg_sites<-deg_sites[match(row.names(deg_results),row.names(deg_sites)),]
deg_sites=deg_sites[,names(deg_sites) %in% names(deg_results)]
miss_deg=names(deg_results)[names(deg_results) %in% names(deg_sites)==F]
for (taxon in miss_deg){
  deg_sites[[taxon]]<-0
}
deg_sites=data.table::setcolorder(deg_sites,neworder = names(deg_results))
deg_OEs=OE_calc(results_data = deg_results,PA=deg_sites, threshold)

boxplot(MRF_ref_OEs$OtoE ,at=1,xlim=c(0,4),ylim=c(0,1.8))
boxplot(MRF_test_OEs$OtoE,add=T,at=2)
boxplot(deg_OEs$OtoE,add=T,at=3)
abline(h=median(MRF_ref_OEs$OtoE),lty=2)
abline(h=median(MRF_test_OEs$OtoE),lty=2)
abline(h=median(deg_OEs$OtoE),lty=2)

axis(side=1,at=1:3,labels=c('Reference','Probabilistic',"'Impaired'"))


savp(10,8,'C://Users//andrew.caudillo.BUGLAB-I9//Box//NAMC//Research Projects//AIM//IncreaserDecreaser_OE//MRF_OE//OE_types_boxes_Pc0_251023.png')

plot(density(MRF_ref_OEs$OtoE),ty='l',main='O/E densities',xlim=c(0,1.6),xlab='O/E score',ylim=c(0,2.2))
lines(density(MRF_test_OEs$OtoE),col='red')
lines(density(deg_OEs$OtoE),col='blue')
legend('topleft',
       leg=c('Reference',
             'Prob.',
             "'Impaired'"),
       lty=rep(1,3),
       col=c('black','red','blue'),
       bty='n')
savp(10,8,'C://Users//andrew.caudillo.BUGLAB-I9//Box//NAMC//Research Projects//AIM//IncreaserDecreaser_OE//MRF_OE//OE_types_densities_Pc0_251023.png')

#ANOVA
OE_dat_comp=data.frame(Group=as.factor(c(rep('Ref.',nrow(MRF_ref_OEs)),
                               rep('Prob.',nrow(MRF_test_OEs)),
                               rep('Deg.',nrow(deg_OEs)))),
                       OE=c(MRF_ref_OEs$OtoE,MRF_test_OEs$OtoE,deg_OEs$OtoE))
anova_res=aov(OE~Group,data=OE_dat_comp)
summary(anova_res)

median(MRF_test_OEs$OtoE)
median(deg_OEs$OtoE)
median(MRF_ref_OEs$OtoE)

t.test(MRF_ref_OEs$OtoE, MRF_test_OEs$OtoE)




#pois binomial distrib



CI_dat=data.frame(X=rep(NA,ncol(OthEco)))#UL=rep(NA,ncol(OthEco)),LL=rep(NA,ncol(OthEco)),
row.names(CI_dat)=names(OthEco)
for(i in 1:ncol(OthEco)){
  x=poibin::qpoibin(qq=c(0.025, 0.975),pp=OthEco[,i])/nrow(OthEco)
  #CI_dat$LL[i]=x[1]
  #CI_dat$UL[i]=x[2]
  CI_dat$X[i]=x*100
}


if(0){
  # #k-folds
  library(caret)

  k <- 5  # choose 5- or 10-fold depending on sample size
  folds <- createFolds(1:nrow(train), k = k, list = TRUE, returnTrain = TRUE)
  perform_metrics=numeric(length(folds))
  #
  cv_results <- data.frame(fold =foldscv_results <- data.frame(fold = integer(), R2 = numeric(), RMSE = numeric()))
  #
  for (i in seq_along(folds)) {
    cat("Running fold", i, "of", k, "\n")
    train_data<-train[-folds[[i]],]
    val_dat<-train[folds[[i]],]
    train_idx <- folds[[i]]
    train_fold <- train[train_idx, ]
    test_fold  <- train[-train_idx, ]
    message('fitting model')
    #   # fit model
    fit_fold <- rfsrc(
      rf_formula,
      data = train_data,
      ntree = 1000,
      nodesize = 5,
      mtry = 4,
      nsplit = 10,
      nodedepth = 5
    )
    #
    #   # predict on hold-out
    pred_fold <- predict(fit_fold, newdata = val_dat)
    pred_fold_preds=as.data.frame(sapply(pred_fold$regrOutput, function(x) x$predicted))
    row.names(pred_fold_preds)=row.names(val_dat)
    num_cols <- sapply(val_dat[, 1:(pred_index-2)], is.numeric)
    val_numeric <- val_dat[, 1:(pred_index-2)][, num_cols]



    obs_numeric  <- val_numeric[, sapply(val_numeric, is.numeric)]
    pred_numeric <- pred_fold_preds[, sapply(pred_fold_preds, is.numeric)]

    # Step 2: Check dimensions match
    if(!all(dim(obs_numeric) == dim(pred_numeric))) {
      stop("Observed and predicted data frames have different dimensions")
    }

    # Step 3: Compute RMSE
    rmse <-  sqrt(mean((as.matrix(obs_numeric) - as.matrix(pred_numeric))^2, na.rm = TRUE))
    #
    #   # compute observed and expected richness

    E <- rowSums(pred_fold_preds)

    O <- rowSums(val_dat[,1:(pred_index-2)])
    #E_test=rowSums(pred_fold_preds)




    #
    R2 <- 1 - sum((O - E)^2) / sum((O - mean(O))^2)
    RMSE <- sqrt(mean((O - E)^2))
    #
    cv_results <- rbind(cv_results,
                        data.frame(fold = i,
                                   taxa_RMSE = rmse,
                                   richness_RMSE = RMSE,
                                   R2 = R2))
  }
  #
  # # summarize
  cv_summary <- cv_results %>%
    summarise(mean_R2 = mean(R2, na.rm=TRUE),
              sd_R2   = sd(R2, na.rm=TRUE),
              mean_RMSE = mean(richness_RMSE, na.rm=TRUE),
              sd_RMSE   = sd(richness_RMSE, na.rm=TRUE))
  print(cv_summary)

  cv_results
  mean(perform_metrics)

  #individual RF to tell overall variable importance
  #across responses.
  response_vars=names(train)[1:(pred_index-2)]
  #pred_vars=names(train)[pred_index:ncol(train)]
  pred_vars=names(train[,topvars])
  vimplist<-list()

  for (resp in response_vars) {
    form <- as.formula(paste(resp, "~", paste(pred_vars, collapse = "+")))
    m <- rfsrc(form, data = train, importance = TRUE, ntree = 500)
    vimplist[[resp]] <- m$importance
    print(resp)
  }

  vimp_df <- do.call(cbind, vimplist)
  vimp_mean <- rowMeans(vimp_df, na.rm = TRUE)
  sort(vimp_mean, decreasing = TRUE)[1:10]
  vimp_sorted <- sort(vimp_mean, decreasing = TRUE)
  plot(vimp_sorted, type = "b", main = "Variable Importance (Aggregated)")
  abline(h = quantile(vimp_sorted, 0.25), col = "red", lty = 2)
  top_vars <- names(vimp_sorted)[vimp_sorted > quantile(vimp_sorted, 0.5)]


  vimp_summary <- data.frame(
    Predictor = row.names(vimp_df),
    Mean = apply(vimp_df, 1, mean, na.rm = TRUE),
    Median = apply(vimp_df, 1, median, na.rm = TRUE),
    Max = apply(vimp_df, 1, max, na.rm = TRUE),
    Consistency = apply(vimp_df, 1, function(x) mean(x > quantile(x, 0.9), na.rm = TRUE))
  )

  library(ggplot2)
  ggplot(vimp_summary, aes(x = reorder(Predictor, Mean), y = Mean)) +
    geom_col(fill = "forestgreen", alpha = 0.7) +
    coord_flip() +
    labs(x = "Predictor", y = "Mean Variable Importance",
         title = "Mean Variable Importance Across Taxa") +
    theme_minimal(base_size = 13)

  savp(10,8,'C://Users//andrew.caudillo.BUGLAB-I9//Box//NAMC//Research Projects//AIM//IncreaserDecreaser_OE//MRF_OE//Varimp_MRF_251023.png')

  ggplot(vimp_summary, aes(x = Mean, y = Max, label = Predictor)) +
    geom_point(color = "dodgerblue3", size = 3) +
    ggrepel::geom_text_repel(size=3.25) +
    labs(x = "Mean Importance (Consistency Across Taxa)",
         y = "Max Importance (Peak Effect for a Taxon)",
         title = "Predictor Roles in Multivariate RF Model") +
    theme_minimal(base_size = 13)
  savp(10,8,'C://Users//andrew.caudillo.BUGLAB-I9//Box//NAMC//Research Projects//AIM//IncreaserDecreaser_OE//MRF_OE//Predictor_roles_MRF_251023.png')

}

if(0){
  #fxn to get overal OOB error for a MRF.
  #inherently hard to do from the rfsrc function
  #because of so many response variables, so rfsrc does not compute oob.
  summarize_mrf_oob <- function(rf_model, threshold = 0.5, weight_by = c("none", "prevalence", "variance")) {
    weight_by <- match.arg(weight_by)

    # --- Extract observed and predicted OOB ---
    obs_mat <- rf_model$yvar
    pred_oob <- data.frame(sapply(rf_model$regrOutput, function(x) x$predicted.oob))

    if (is.null(pred_oob)) stop("rf_model$predicted.oob is NULL — make sure you fit with OOB predictions enabled.")
    if (!all(dim(obs_mat) == dim(pred_oob))) stop("Dimensions of observed and OOB-predicted data do not match.")

    n_taxa <- ncol(obs_mat)
    taxa_names <- colnames(obs_mat)

    # --- Detect binary responses ---
    is_binary <- apply(obs_mat, 2, function(x) {
      ux <- unique(na.omit(x))
      length(ux) == 2 && all(sort(ux) == c(0, 1))
    })

    # --- Initialize containers ---
    taxon_misclass <- rep(NA_real_, n_taxa)
    taxon_rmse <- rep(NA_real_, n_taxa)
    taxon_r2 <- rep(NA_real_, n_taxa)

    for (j in seq_len(n_taxa)) {
      y <- obs_mat[, j]
      yhat <- pred_oob[, j]

      if (var(y, na.rm = TRUE) == 0 || all(is.na(y))) next

      if (is_binary[j]) {
        preds_class <- ifelse(yhat > threshold, 1, 0)
        taxon_misclass[j] <- mean(preds_class != y, na.rm = TRUE)
        taxon_rmse[j] <- sqrt(mean((y - yhat)^2, na.rm = TRUE))
      } else {
        taxon_rmse[j] <- sqrt(mean((y - yhat)^2, na.rm = TRUE))
      }

      ss_res <- sum((y - yhat)^2, na.rm = TRUE)
      ss_tot <- sum((y - mean(y, na.rm = TRUE))^2, na.rm = TRUE)
      taxon_r2[j] <- 1 - (ss_res / ss_tot)
    }

    per_taxon <- data.frame(
      Taxon = taxa_names,
      is_binary = is_binary,
      Misclass = taxon_misclass,
      RMSE = taxon_rmse,
      R2 = taxon_r2,
      stringsAsFactors = FALSE
    )

    # --- Define weights if requested ---
    weights <- rep(1, n_taxa)
    if (weight_by == "prevalence") {
      weights <- colMeans(obs_mat == 1, na.rm = TRUE)
    } else if (weight_by == "variance") {
      weights <- apply(obs_mat, 2, var, na.rm = TRUE)
    }
    weights[is.na(weights)] <- 0
    weights <- weights / sum(weights, na.rm = TRUE)

    # --- Aggregate per-taxon metrics ---
    overall_mean_R2 <- mean(per_taxon$R2, na.rm = TRUE)
    overall_median_R2 <- median(per_taxon$R2, na.rm = TRUE)
    weighted_mean_R2 <- sum(per_taxon$R2 * weights, na.rm = TRUE)
    overall_mean_RMSE <- mean(per_taxon$RMSE, na.rm = TRUE)

    # --- Community-level metrics (richness) ---
    Fe_oob <- rowSums(pred_oob, na.rm = TRUE)          # expected richness (probability sum)
    Fo_obs <- rowSums(obs_mat, na.rm = TRUE)           # observed richness

    comm_rmse <- sqrt(mean((Fo_obs - Fe_oob)^2, na.rm = TRUE))
    ss_res_comm <- sum((Fo_obs - Fe_oob)^2, na.rm = TRUE)
    ss_tot_comm <- sum((Fo_obs - mean(Fo_obs, na.rm = TRUE))^2, na.rm = TRUE)
    comm_r2 <- 1 - (ss_res_comm / ss_tot_comm)

    # --- Return summary list ---
    return(list(
      per_taxon = per_taxon,
      summary = data.frame(
        Metric = c("Mean_R2", "Median_R2", "Weighted_R2", "Mean_RMSE", "Community_R2", "Community_RMSE"),
        Value = c(overall_mean_R2, overall_median_R2, weighted_mean_R2, overall_mean_RMSE, comm_r2, comm_rmse)
      )
    ))
  }
  oob=summarize_mrf_oob(rf_model = rf_model,threshold = 0.5,weight_by = 'none')
  oob$summary
  #savp(10,8,'C://Users//andrew.caudillo.BUGLAB-I9//Box//NAMC//Research Projects//AIM//IncreaserDecreaser_OE//MRF_OE//Predictor_roles_MRF.png')

  #This is variable selection for the
  trees=max.subtree.rfsrc(rf_model,conservative = F)

}
