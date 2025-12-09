
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
  #remove sites for which predictors cannot be calculated
  ben_dat<-ben_dat[ben_dat$X %in% c(116830,132243,132807, 185022,116435,116436,117875,116849)==F,] #last 3 are new
  names(ben_dat)[1]<-'sampleId'


  # Wide-format predictor values — ensures only one column per predictor_id
  predictor_wide <- ref_dat %>%
     select(sampleId, abbreviation, predictor_value) %>%
     distinct() %>%
     pivot_wider(names_from = abbreviation, values_from = predictor_value)

 #join the two datasets together
  ben_dat=plyr::join(ben_dat,predictor_wide)
#omit Carabidae and Curuclionidae
  ben_dat<-ben_dat[,names(ben_dat) %in% c('Carabidae','Curculionidae')==F]
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
  #create simple case weights for sites with fewer bugs to be removed more often
  cases=1 / (rowSums(train[, taxa_names]) + 1)
  case_wts=cases / sum(cases)


  #this model has already been fit with the best hyperparameters
  #as found by the tune() function
  rf_model <- rfsrc(
    formula=rf_formula,
    data=train,
    ntree = 500,
    mtry=4,
    importance = 'permute',
    nodedepth = 5,
    nodesize=20,
    nsplit=1,
    case.wt = case_wts
    )
#view OOB errors for each response
OOB_pred_each=get.mv.error(rf_model)
mean(OOB_pred_each)
#pseudo RMSE
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

All_ratios=read.csv('C://Users//andrew.caudillo.BUGLAB-I9//Box//NAMC//Research Projects//AIM//IncreaserDecreaser_OE//regional_responses_MRF_251114.csv')
write.csv(All_ratios,'C://Users//andrew.caudillo.BUGLAB-I9//Box//NAMC//Research Projects//AIM//IncreaserDecreaser_OE//regional_responses_MRF_251117.csv')
#calculating Fo and Fe, then subsetting out those which are uncommon
  #and/or do not have traits
  Fos=colSums(train[,1:(pred_index-2)])
  Fes=colSums(results)
  #Fos_small=Fos[names(Fos) %in% taxa_notraits_rare$taxon==F]
  #Fes_small=Fes[names(Fes) %in% taxa_notraits_rare$taxon==F]
  ratio=Fos/Fes
  ratio_dat=data.frame(Fo=Fos,Fe=Fes)
  #ratio_small=ratio[names(ratio) %in% taxa_notraits_rare$taxon==F]
  boxplot(ratio,ylab='Fo/Fe ratio')
  savp(10,8,'C://Users//andrew.caudillo.BUGLAB-I9//Box//NAMC//Research Projects//AIM//IncreaserDecreaser_OE//MRF_OE//Fo_Fe_boxplot_251114.png')
  plot(Fes,Fos,ylab='Fo',xlab='Fe',main='Reference Fo vs Fe')
  abline(0,1,col='red')
  savp(10,8,'C://Users//andrew.caudillo.BUGLAB-I9//Box//NAMC//Research Projects//AIM//IncreaserDecreaser_OE//MRF_OE//Ref_FovsFe_251114.png')


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

  test_dat=read.csv('C://Users//andrew.caudillo.BUGLAB-I9//Box//NAMC//Research Projects//AIM//IncreaserDecreaser_OE//MRF_Prob_Os.csv')
  #test_dat$Tanypodinae = test_dat$Tanypodinae + test_dat$Tanypodinae.1;test_dat<-test_dat[,names(test_dat) %in% 'Tanypodinae.1'==F,]
  test_dat<-test_dat[,-1]
  names(test_dat)[1]<-'sampleId'
  failed_sites<-read.csv('C://Users//andrew.caudillo.BUGLAB-I9//Box//NAMC//Research Projects//AIM//IncreaserDecreaser_OE//failed_sites.csv')
  failed_dat=test_dat[test_dat$sampleId %in% failed_sites$sampleId,]
  test_dat<-test_dat[test_dat$sampleId %in% failed_sites$sampleId==F,]

  predictor_test_wide <- test_pred_dat %>%
    select(sampleId, abbreviation, predictor_value) %>%
    distinct() %>%
    pivot_wider(names_from = abbreviation, values_from = predictor_value)

  # 3. Merge: join both by sample_id (one row per sample)
  #test_dat<-test_dat[test_dat$sampleId %in% failed_sites$sampleId==F,]
  test_dat=plyr::join(test_dat,predictor_test_wide)
#drop Curcilionid and Carabid
  test_dat<-test_dat[,names(test_dat) %in% c('Carabidae','Curculionidae')==F]
  pred_test_index<-which(names(test_dat)==names(predictor_test_wide)[2])
  rownames(test_dat) <- test_dat$sampleId;test_dat<-test_dat[,names(test_dat) %in% 'sampleId'==F]

  new_sites <- test_dat[,(pred_test_index-1):ncol(test_dat)]
  row.names(new_sites)=row.names(test_dat)
  test_site_names<-row.names(new_sites)
  new_sites<-as.data.frame(lapply(new_sites,as.numeric))
  row.names(new_sites)<-test_site_names
#predict the model onto the prob. sites
  pred <- predict.rfsrc(rf_model, newdata = new_sites)
  # Extract predicted probabilities for each taxon at prob sites
  #no out of bag for test sites because the model uses all trees to make a choice.
  pred_probs <- as.data.frame(sapply(pred$regrOutput, function(x) x$predicted))
  row.names(pred_probs)<-row.names(new_sites)
  #write.csv(pred_probs,'C://Users//andrew.caudillo.BUGLAB-I9//Box//NAMC//Research Projects//AIM//IncreaserDecreaser_OE//MRF_Prob_Pcs.csv')
#calculating Fo and Fe
  Pr_PA=test_dat[,names(test_dat) %in% names(pred_probs)]
  Pr_PA=as.data.frame(ifelse(Pr_PA >0,1,0))
  #write.csv(Pr_PA, 'C://Users//andrew.caudillo.BUGLAB-I9//Box//NAMC//Research Projects//AIM//IncreaserDecreaser_OE//MRF_Prob_Os.csv')
  PFos=colSums(Pr_PA)
  PFes=colSums(pred_probs)

  Xeric_P=c(171583, 171664, 171668, 171669, 171777, 173398, 173404,
            173407, 173408, 173409, 178667, 178671, 178676, 184764,
            184768, 185060, 185063, 185699, 185702, 187162, 187182,
            187607,187662, 187779, 187780, 187849, 190728, 190762,
            190916, 191507)
  Pr_PA_oth=Pr_PA[row.names(Pr_PA) %in% Xeric_P==F,]
  pred_probs_oth=pred_probs[row.names(pred_probs) %in% Xeric_P==F,]

  Pr_PA_x=Pr_PA[row.names(Pr_PA) %in% Xeric_P,]
  pred_probs_x=pred_probs[row.names(pred_probs) %in% Xeric_P,]

  Pratio=PFos/PFes
  Poth_Fo=colSums(Pr_PA_oth);Poth_Fe=colSums(pred_probs_oth)
  Poth_ratio=Poth_Fo/Poth_Fe
  Prob_combined_FoFes=as.data.frame(cbind(PFos,PFes))
  #Pratio_small=Pratio[names(Pratio) %in% taxa_notraits_rare$taxon==F]
  #PFos_small=PFos[names(PFos) %in% taxa_notraits_rare$taxon==F]
  #PFes_small=PFes[names(PFes) %in% taxa_notraits_rare$taxon==F]

  #Poth_Fosmall=Poth_Fo[names(Poth_Fo) %in% taxa_notraits_rare$taxon==F]
  #Poth_Fesmall=Poth_Fe[names(Poth_Fe) %in% taxa_notraits_rare$taxon==F]
  Px_Fo=colSums(Pr_PA_x);Px_Fe=colSums(pred_probs_x)
  #Px_Fosmall=Px_Fo[names(Px_Fo) %in% taxa_notraits_rare$taxon==F]
  #Px_Fesmall=Px_Fe[names(Px_Fe) %in% taxa_notraits_rare$taxon==F]
  #Poth_ratio=Poth_Fo/Poth_Fe
  Px_ratio=Px_Fo/Px_Fe
  #Poth_ratio_small=Poth_Fosmall/Poth_Fesmall
  #Px_ratio_small=Px_Fosmall/Poth_Fesmall
  plot(PFes,PFos,ylab='Fo',xlab='Fe',main='Prob. Fo vs Fe')
  abline(0,1,col='red')

  Poth_plotdat=data.frame(Fo=Poth_Fo,Fe=Poth_Fe)
  Px_plotdat=data.frame(Fo=Px_Fo,Fe=Px_Fe)
  boxplot(Pratio,ylab='Fo/Fe ratio')
  savp(10,8,'C://Users//andrew.caudillo.BUGLAB-I9//Box//NAMC//Research Projects//AIM//IncreaserDecreaser_OE//MRF_OE//Prob_FoFe_ratio_251117.png')
  savp(10,8,'C://Users//andrew.caudillo.BUGLAB-I9//Box//NAMC//Research Projects//AIM//IncreaserDecreaser_OE//MRF_OE//Prob_FovsFe_box_25117.png')


  boxplot(Pratio,at=2,xlim=c(1,4),ylab='Fo/Fe')
  boxplot(ratio,at=3,xlim=c(1,4),add=T)

  axis(side=1,at=c(2,3),labels = c('Prob','Train'))


  savp(10,8,'C://Users//andrew.caudillo.BUGLAB-I9//Box//NAMC//Research Projects//AIM//IncreaserDecreaser_OE//MRF_OE//Box_compare_FoFes_251023.png')

  boxplot(ratio,at=1,col='purple3',ylim=c(0,max(Pratio[is.finite(Pratio)])),xlim=c(0,3),
          ylab='Fo/Fe ratio')
  boxplot(Pratio,at=2,col='yellow3',add=T)
points(x=rep(2, length(Pratio[Pratio>=8])),y=(Pratio[Pratio>=8]),bg=c('blue','red'),pch=21)

legend('topright',
       leg=c('Reference',
             'Probabilistic'),
       pch=rep(22,2),
       pt.bg=c('purple4',
               'yellow3'),
       bty='n',
       cex=0.8)
legend('topleft',
       leg=c('Mayatrichia','Cambaridae'),
       pch=c(21,21),
       pt.bg=c('red','blue'),
       bty='n',
       cex=0.8)
savp(10,8,'C://Users//andrew.caudillo.BUGLAB-I9//Box//NAMC//Research Projects//AIM//IncreaserDecreaser_OE//MRF_OE//regional_boxplot_compare_251125_colored.png')
graphics.off()
boxplot(Oth_ratio,at=1,xlim=c(0,5),ylim=c(0,max(Px_ratio[is.finite(Px_ratio)])),col='purple3',ylab='Fo/Fe ratio')
boxplot(Poth_ratio,at=2,add=T,col='yellow3')
boxplot(Xer_ratio,at=3,add=T,col='purple3')
boxplot(Px_ratio,at=4,add=T,col='yellow3')
points(x=rep(4,6),y=Px_ratio[Px_ratio >= 20 & is.finite(Px_ratio)],bg=c('red','blue','purple',
                                                                        'black','dodgerblue3','orange'),pch=21)
mtext(text = c("Other", "Xeric"), side = 1, line = 1, at = c(1.5, 3.5), cex = 1)
legend('topleft',
       leg=c('Cambaridae',
             'Paracloedes',
             'Gammarus',
             'Caecidotea',
             'Lestes',
             'Mayatrichia'),
       pt.bg=c('blue','orange','purple',
               'red','black','dodgerblue3'),
       pch=rep(21,6),
       bty='n',
       cex=0.8,
       ncol=2)
legend('topright',
       leg=c('Reference',
             'Probabilistic'),
       pch=rep(22,2),
       pt.bg=c('purple4',
               'yellow3'),
       bty='n',
       cex=0.8)
savp(10,8,'C://Users//andrew.caudillo.BUGLAB-I9//Box//NAMC//Research Projects//AIM//IncreaserDecreaser_OE//MRF_OE//ecoregion_boxplot_compare_251125_colored.png')
# ---------------------------
#creating O/E scores, quasi-RIVPACS style
#set the Pc threshold, in this case, 0.5
threshold <- 0
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

All_Xer=c(186872,187130,187162,187182,
          187183,187222,187250,187275,
          187290,187314,187348,187392,
          187448,187455,187456,187468,
          187542,187587,187588,187591,
          187604,187607,187638,187661,
          187662,187709,187721,187722,
          187779,187780,187831,187849,
          190552,190553,190717,190728,
          190736,190748,190762,190846,
          190880,190916,190933,190934,
          191067,191068,191069,191224,
          191231,191401,191446,191464,
          191507,191508,171582,171583,
          171663,171664,171668,171669,
          171777,173398,173404,173405,
          173406,173407,173408,173409,
          173578,173579,176421,176422,
          178188,178189,178418,178667,
          178671,178676,184757,184764,
          184768,185060,185062,185063,
          185699,185702,190591,190619,
          190620,190641,190642,190647,
          190650,190654,190667,190682,
          190709,116386,116389,116849,
          117866,117867,117868,117869,
          117882,131904,132140,132431,
          132660,132665,132667,133381,
          133382,133743,133816,133820)
#run it on ref and prob. sites
EastXer=results[row.names(results) %in% All_Xer,]
OthEco=results[row.names(results) %in% All_Xer==F,]
XerTrain=train[row.names(train) %in% All_Xer,]
XerTrain=XerTrain[,1:(ncol(XerTrain)-15)]
OthTrain=train[row.names(train) %in%All_Xer==F,]
OthTrain=OthTrain[,1:(ncol(OthTrain)-15)]
XerFo=colSums(XerTrain)
XerFe=colSums(EastXer)
Xer_ratio=XerFo/XerFe


OthFo=colSums(OthTrain)
OthFe=colSums(OthEco)
Oth_ratio=OthFo/OthFe
#Oth_ratio_small=Oth_ratio[names(Oth_ratio) %in% names(taxa_notraits_rare)==F]
#Xer_ratio_small=Xer_ratio[names(Xer_ratio) %in% names(taxa_notraits_rare)==F]


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
graphics.off()
#combined_dat=data.frame(Fo=Fos_small,Fe=Fes_small)
#4 panel plot showing Ref vs Prob
A<-ggplot(data=OthEco_plotdat,aes(y=Fo,x=Fe))+geom_point()+geom_abline(intercept = 0,slope = 1,col='red')+ggtitle('Other Ecoregions')+
  xlim(0,max(Poth_Fe))+ylim(0,max(Poth_Fo))+
  annotate("text",
           x = -Inf, y = Inf, # Position at top-left corner
           label = 'Reference',
           hjust = 0, vjust = 1, # Justify text relative to corner
           size = 5, color = "black")
B<-ggplot(data=EastXer_plotdat,aes(y=Fo,x=Fe))+geom_point()+geom_abline(intercept = 0,slope = 1,col='red')+ggtitle('Eastern Xeric')+
  annotate("text",
           x = -Inf, y = Inf, # Position at top-left corner
           label = 'Reference',
           hjust = 0, vjust = 1, # Justify text relative to corner
           size = 5, color = "black")
#ggplot(data=combined_dat,aes(x=Fe,y=Fo))+geom_point()+geom_abline(intercept = 0,slope = 1,col='red')+
 # xlim(0,max(Px_plotdat$Fe))
#savp(10,8,'C://Users//andrew.caudillo.BUGLAB-I9//Box//NAMC//Research Projects//AIM//IncreaserDecreaser_OE//MRF_OE//MRF_allsites_FoFe.png')

#getting highlight taxa for showing some inc/decs
Poth_plotdat$col=ifelse(row.names(Poth_plotdat) %in% c('Serratella','Ephemerella','Antocha','DRUNELLA_DODDSI'),
                     'red',NA)
Px_plotdat$col=ifelse(row.names(Px_plotdat) %in% c('Serratella','Ephemerella','Antocha','DRUNELLA_DODDSI'),
                      'red',NA)
Poth_highlight=Poth_plotdat[which(Poth_plotdat$col=='red'),]
Poth_highlight$taxon=row.names(Poth_highlight)
Poth_highlight$taxon[Poth_highlight$taxon=='DRUNELLA_DODDSI']<-'D. dodsii'
Px_highlight=Px_plotdat[which(Px_plotdat$col=='red'),]
Px_highlight$taxon=row.names(Px_highlight)
Px_highlight$taxon[Px_highlight$taxon=='DRUNELLA_DODDSI']<-'D. dodsii'
C<-ggplot(data=Poth_plotdat,aes(y=Fo,x=Fe))+geom_point()+geom_abline(intercept = 0,slope = 1,col='red')+
  geom_point(data=Poth_highlight,aes(x=Fe,y=Fo,color=taxon))+
  scale_color_manual(values=c('Antocha'='red',
                              'Ephemerella'='purple',
                              'Serratella'='orange',
                              'D. dodsii' = 'dodgerblue'))+
  annotate("text",
           x = -Inf, y = Inf, # Position at top-left corner
           label = 'Probabilistic',
           hjust = 0, vjust = 1, # Justify text relative to corner
           size = 5, color = "black")+
  theme(legend.position = "bottom")
D<-ggplot(data=Px_plotdat,aes(y=Fo,x=Fe))+geom_point()+geom_abline(intercept = 0,slope = 1,col='red')+
  ylim(0,max(EastXer_plotdat$Fo))+xlim(0,max(EastXer_plotdat$Fe))+
  geom_point(data=Px_highlight,aes(x=Fe,y=Fo,color=taxon))+
  scale_color_manual(values=c('Antocha'='red',
                              'Ephemerella'='purple',
                              'Serratella'='orange',
                              'D. dodsii' = 'dodgerblue'))+
  annotate("text",
           x = -Inf, y = Inf, # Position at top-left corner
           label = 'Probabilistic',
           hjust = 0, vjust = 1, # Justify text relative to corner
           size = 5, color = "black")+
  theme(legend.position = "none")
#getting a shared legend for this large panel is hard.
#use a function
get_only_legend <- function(plot) {

  # get tabular interpretation of plot
  plot_table <- ggplot_gtable(ggplot_build(plot))

  #  Mark only legend in plot
  legend_plot <- which(sapply(plot_table$grobs, function(x) x$name) == "guide-box")

  # extract legend
  legend_tiles <- plot_table$grobs[[legend_plot]]

  # return legend
  return(legend_tiles)
}
tax_legend=get_only_legend(C)
#redefine C with no legend
C<-ggplot(data=Poth_plotdat,aes(y=Fo,x=Fe))+geom_point()+geom_abline(intercept = 0,slope = 1,col='red')+
  geom_point(data=Poth_highlight,aes(x=Fe,y=Fo,color=taxon))+
  scale_color_manual(values=c('Antocha'='red',
                              'Ephemerella'='purple',
                              'Serratella'='orange',
                              'D. dodsii' = 'dodgerblue'))+
  annotate("text",
           x = -Inf, y = Inf, # Position at top-left corner
           label = 'Probabilistic',
           hjust = 0, vjust = 1, # Justify text relative to corner
           size = 5, color = "black")+
  theme(legend.position = "none")
#define the plot, but don't plot it
cmbin_plot=gridExtra::grid.arrange(A,B,C,D,ncol=2)
#plot the final graph with shared legend
gridExtra::grid.arrange(cmbin_plot,tax_legend,heights=c(10,1))

savp(10,8,'C://Users//andrew.caudillo.BUGLAB-I9//Box//NAMC//Research Projects//AIM//IncreaserDecreaser_OE//MRF_OE//MRF_ecoregions_FoFe_sitecompare_251125_Colored.png')
#get O/E scores for all sites / ecoregion subsets
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
boxplot(OthEco_OE$OtoE,at=1,xlim=c(0,5),ylab='O/E score',col='purple4',ylim=c(0,2),
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
savp(10,8,'C://Users//andrew.caudillo.BUGLAB-I9//Box//NAMC//Research Projects//AIM//IncreaserDecreaser_OE//MRF_OE//OtoE_boxes_ecoregions_251117.png')
#savp(10,8,'C://Users//andrew.caudillo.BUGLAB-I9//Box//NAMC//Research Projects//AIM//IncreaserDecreaser_OE//MRF_OE//OtoE_boxes_Pc05_251106.png')
### This is just looking at O/E performance and metrics surrounding it

boxplot(Pratio,at=2,col='yellow3',ylim=c(0,70),xlim=c(0,5))
boxplot(Px_ratio,at=3,add=T,col='purple4')
boxplot(Poth_ratio,at=4,add=T,col='yellow3')
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



#get RMSE from model
rmse_train <- sqrt(mean((MRF_ref_OEs$OtoE - 1)^2)); rmse_test <- sqrt(mean((MRF_test_OEs$OtoE - 1)^2))
rmse_xer=sqrt(mean((XerOE$OtoE - 1)^2))
rmse_oth=sqrt(mean((OthEco_OE$OtoE - 1)^2))
rmse_train
rmse_test


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

savp(10,8,'C://Users//andrew.caudillo.BUGLAB-I9//Box//NAMC//Research Projects//AIM//IncreaserDecreaser_OE//MRF_OE//OtoE_boxes_Pc0_251117.png')



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
#optional, just testing for overfit.

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
#this gets the confidence intervals for the predictions
#quick and simple vectorized calculations
#copy/paste results into a spreadsheet with overall responses and
#should have enough info for a descriptive table.

CI_dat <- data.frame(
  Fe_mean = rep(NA, ncol(results)),
  LL      = rep(NA, ncol(results)),
  UL      = rep(NA, ncol(results))
)
row.names(CI_dat)=names(results)
n=nrow(results)
for(i in 1:ncol(results)){
  x=poibin::qpoibin(qq=c(0.025, 0.975),pp=results[,i])
  Fo
  Fe_mean = sum(results[,i])/n
  ci_Fe=x/n
  CI_dat$LL[i]=x[1]/n
  CI_dat$UL[i]=x[2]/n
  CI_dat$Fe_mean[i]<-Fe_mean
  #CI_dat$X[i]=x
}
CI_dat
clipr::write_clip(CI_dat)

#this just shows how many times that taxa appeared
#so you can compare it to the CIs
pr_sums=colSums(Pr_PA[,taxa_names])
pr_sums_df=data.frame(OTU=names(pr_sums),
                      N=pr_sums)



