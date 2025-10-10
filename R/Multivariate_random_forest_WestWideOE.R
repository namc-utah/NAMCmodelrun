
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
  ben_dat<-ben_dat[ben_dat$X %in% c(116830,132243,132807, 185022)==F,]
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

  # Fit Multivariate Random Forest
  # 'regr.multivariate' formula, response = multiple taxa columns
  #writ ehte formula
  formula_string <- paste0(
    "cbind(", paste(taxa_names, collapse = ", "), ") ~ ",
    paste(names(train[pred_index:ncol(train)]), collapse = " + ")
  )

  # Convert string to formula using parse + eval
  #Else MRF will not read it correctly
  rf_formula <- eval(parse(text = formula_string))
  #this model has already been fit with the best hyperparameters
  #as found by the tune() function
  rf_model <- rfsrc(
    formula=rf_formula,
    data=train,
    ntree = 500,
    mtry=14,
    importance = T
  )
#view OOB errors for each response
get.mv.error(rf_model)
  #get the predicted Pcs for each taxon here
  results<-data.frame(sapply(rf_model$regrOutput, function(x) x$predicted))
  #set row names
  row.names(results)=row.names(train)
  #adjust hyperparamters for best model
  #if set to if(0), it has already been run!
  if(0){
  tune(rf_formula,
             train,
             doBest = T,
             improve=1e-2, #define what improvement means for OOB error
       sampsize = 400, #number of samples for the subsamples
       mtryStart = max(2, round(ncol(train) / 4)), #default from rfsrc doc
       trace=T, #show the progress in console
       nodesizeTry = c(5,10,15),
       ntreeTry = 100) #a small number of trees for memory usage
  }
  #calculating Fo and Fe, then subsetting out those which are uncommon
  #and/or do not have traits
  Fos=colSums(train[,taxa_names])
  Fes=colSums(results)
  Fos_small=Fos[names(Fos) %in% taxa_notraits_rare$taxon==F]
  Fes_small=Fes[names(Fes) %in% taxa_notraits_rare$taxon==F]
  ratio=Fos/Fes
  ratio_small=ratio[names(ratio) %in% taxa_notraits_rare$taxon==F]
  boxplot(ratio_small)
  plot(Fes_small,Fos_small,ylab='Fo',xlab='Fe',main='Reference Fo vs Fe')
  abline(0,1,col='red')
  savp(10,8,'C://Users//andrew.caudillo.BUGLAB-I9//Box//NAMC//Research Projects//AIM//IncreaserDecreaser_OE//MRF_OE//Ref_FovsFe.png')


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
  #failed_sites<-read.csv('C://Users//andrew.caudillo.BUGLAB-I9//Box//NAMC//Research Projects//AIM//IncreaserDecreaser_OE//failed_sites.csv')
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

  new_sites<-as.data.frame(lapply(new_sites,as.numeric))
#predict the model onto the prob. sites
  pred <- predict(rf_model, newdata = new_sites)


  # Extract predicted probabilities for each taxon at prob sites
  pred_probs <- as.data.frame(sapply(pred$regrOutput, function(x) x$predicted))
  row.names(pred_probs)<-row.names(new_sites)

#calculating Fo and Fe
  Pr_PA=test_dat[,1:(pred_test_index-2)]
  Pr_PA=as.data.frame(ifelse(Pr_PA >0,1,0))
  PFos=colSums(Pr_PA)
  PFes=colSums(pred_probs)

  Pratio=PFos/PFes
  Pratio_small=Pratio[names(Pratio) %in% taxa_notraits_rare$taxon==F]
  PFos_small=PFos[names(PFos) %in% taxa_notraits_rare$taxon==F]
  PFes_small=PFes[names(PFes) %in% taxa_notraits_rare$taxon==F]
  plot(PFes_small,PFos_small,ylab='Fo',xlab='Fe',main='Prob. Fo vs Fe')
  abline(0,1,col='red')
  savp(10,8,'C://Users//andrew.caudillo.BUGLAB-I9//Box//NAMC//Research Projects//AIM//IncreaserDecreaser_OE//MRF_OE//Prob_FovsFe.png')


  boxplot(log(Pratio_small),at=2,xlim=c(1,4),ylab='log Fo/Fe')
  boxplot(log(ratio_small),at=3,xlim=c(1,4),add=T)

  axis(side=1,at=c(2,3),labels = c('Prob','Train'))


  savp(10,8,'C://Users//andrew.caudillo.BUGLAB-I9//Box//NAMC//Research Projects//AIM//IncreaserDecreaser_OE//MRF_OE//Box_compare_FoFes.png')



  # ---------------------------
#creating O/E scores, quasi-RIVPACS style
#set the Pc threshold, in this case, 0.5
threshold <- 0.5
  #define fxn that will calculate the O/E scores
OE_calc=function(results_data, threshold){
  Os=list()
  Es=list()
  BC=list()
  for(i in 1:nrow(results_data)) {

    cur.prd<-results[i,]; #vector of predicted probs for current sample;
    spdyn<-names(cur.prd)[cur.prd>=threshold];  #subset of taxa with Pi>=Pcutoff;
    cur.prd<-cur.prd[spdyn]
    cur.obs<-train[i,spdyn]; #vector of observed P/A for those taxa;
    Os[[i]]<-sum(cur.obs); #observed richness (O);
    Es[[i]]<-sum(cur.prd); #Expected richness (E);
    BC[[i]]<-sum(abs(cur.obs-cur.prd))/ (sum(cur.obs)+sum(cur.prd)); #BC value;
  }
  OE.dat=as.data.frame(cbind(unlist(Os),unlist(Es),unlist(BC)))

  names(OE.dat)<-c('O','E','BC')
  OE.dat$OtoE=OE.dat$O/OE.dat$E
  OE.dat$E.adj = OE.dat$E*(mean(OE.dat$OtoE))
  OE.dat$OtoE.adj=OE.dat$O / OE.dat$E.adj
  return(OE.dat)
}
#run it on ref and prob. sites
MRF_ref_OEs=OE_calc(results,threshold)
row.names(MRF_ref_OEs)=row.names(results)
MRF_test_OEs=OE_calc(pred_probs,threshold)

### This is just looking at O/E performance and metrics surrounding it
hist(MRF_ref_OEs$OtoE)
hist(MRF_test_OEs$OtoE)

sd(MRF_ref_OEs$OtoE)
sd(MRF_test_OEs$OtoE)

rmse_train <- sqrt(mean((MRF_ref_OEs$OtoE - 1)^2)); rmse_test <- sqrt(mean((MRF_test_OEs$OtoE - 1)^2))

rmse_train
rmse_test

quantile(MRF_ref_OEs$OtoE, probs = c(0.01,0.05,0.1,0.25,0.5,0.75,0.9,0.95,0.99))

boxplot(MRF_ref_OEs$OtoE,at=1,xlim=c(0,3),ylab='O/E')
boxplot(MRF_test_OEs$OtoE,at=2,add=T)
axis(side=1,at=c(1,2),labels=c('Train','Prob.'))
abline(h=median(MRF_ref_OEs$OtoE),lty=2)
abline(h=median(MRF_test_OEs$OtoE),lty=2)

savp(10,8,'C://Users//andrew.caudillo.BUGLAB-I9//Box//NAMC//Research Projects//AIM//IncreaserDecreaser_OE//MRF_OE//OtoE_boxes.png')



MRF_WW_ratios=data.frame(Median_ratio=median(na.omit(ratio_small)),
           Median_ratio_p=median(na.omit(Pratio_small)),
           mean_ratio=mean(na.omit(ratio_small)),
           mean_ratio_p=mean(na.omit(Pratio_small))
)
MRF_WW_ratios

graphics.off()
plot(MRF_ref_OEs$E,MRF_ref_OEs$O,ylab='O',xlab='E')
abline(0,1,col='red')
savp(10,8,'C://Users//andrew.caudillo.BUGLAB-I9//Box//NAMC//Research Projects//AIM//IncreaserDecreaser_OE//MRF_OE//OvsE.png')
OE_lm<-lm(MRF_ref_OEs$O~MRF_ref_OEs$E)
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


plot(density(MRF_ref_OEs$OtoE.adj),ty='l',xlim=c(min(MRF_ref_OEs$OtoE.adj),max(MRF_ref_OEs$OtoE)),xlab='O/E Score',main='')
lines(density(MRF_ref_OEs$OtoE),col='red')
legend('topleft',
       leg=c('O/E',
             'O/E adjusted'),
       lty=rep(1,2),
       col=c('red','black'),
       bty='n')
savp(10,8,'C://Users//andrew.caudillo.BUGLAB-I9//Box//NAMC//Research Projects//AIM//IncreaserDecreaser_OE//MRF_OE//Density_compare.png')


boxplot(MRF_ref_OEs$O,at=1,xlim=c(0,4),ylab='Number of OTUs')
boxplot(MRF_ref_OEs$E,at=2,add=T)
boxplot(MRF_ref_OEs$E.adj,at=3,add=T)
axis(side=1,at=c(1,2,3),labels=c('O','E','E.adj'))
savp(10,8,'C://Users//andrew.caudillo.BUGLAB-I9//Box//NAMC//Research Projects//AIM//IncreaserDecreaser_OE//MRF_OE//O_E_adj_boxes.png')


plot(density(MRF_ref_OEs$E),ty='l',xlim=c(min(MRF_ref_OEs$E),max(MRF_ref_OEs$E.adj)),
     xlab='E value',main='')
lines(density(MRF_ref_OEs$E.adj),col='red')
legend('topright',
       leg=c('E',
             'E adjusted'),
       lty=rep(1,2),
       col=c('black','red'),
       bty='n')
