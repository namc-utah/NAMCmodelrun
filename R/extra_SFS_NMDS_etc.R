#addtl I/D analyses
#traits and %I and %D

#this came from the FoFe_script... file
trait_table=read.csv('C://Users//andrew.caudillo.BUGLAB-I9//Box//NAMC//Research Projects//AIM//IncreaserDecreaser_OE//Updated_w_modelObj//New_for_SFS//Maximized_OTUmatrix_trimmed_more_withvariance.csv')
names(trait_table)[1]<-'taxon'
rbind(O_long,ProbOs_long)
O_for_traits=O_long

first_cols <- names(O_for_traits)[1:656]

# 2. Look at columns 659+ and keep ONLY the ones found in trait_table
target_cols <- names(O_for_traits)[657:ncol(O_for_traits)]
valid_target_cols <- target_cols[target_cols %in% names(trait_table)]

# 3. Combine both sets of columns to subset the data frame
O_long2 <- O_for_traits[, c(first_cols, valid_target_cols)]
#subset just the traits we want for both Prob and Ref
O_long2$status='Ref'
O_long2=O_long2[,names(O_long2) %in% c('HAWK-33','WE-131','WE-209')]
P_for_traits=ProbOs_long
first_cols <- names(P_for_traits)[1:349]
target_cols <- names(P_for_traits)[350:ncol(P_for_traits)]
P_for_traits=P_for_traits[,names(P_for_traits) %in% sites_w_no_bugs==F]
valid_target_cols <- target_cols[target_cols %in% names(trait_table)]
P_long2=P_for_traits[, c(first_cols, valid_target_cols)]
#P_long2=P_for_traits[,which(names(P_for_traits)[352:(ncol(P_for_traits))] %in% names(trait_table))]
P_long2$status='Prob'
P_long2=P_long2[,names(P_long2) %in% odd_sites_samps==F]
#O_long2=O_for_traits[,names(O_for_traits)[659:(ncol(O_for_traits))] %in% names(trait_table)]


#OTU21<-read.csv('C://Users//andrew.caudillo.BUGLAB-I9//Box//NAMC//Research Projects//AIM//IncreaserDecreaser_OE//OTU21.csv')
#trait_table<-trait_table[trait_table$OTU %in% OTU21$taxonOTU,]
names(trait_table)[1]='taxon'

#subset the OTUs and their mapping taxa

P_Results2=P_Results
P_w_traits_full=plyr::join(P_long2,P_Results2[,c('taxon','Regional_Response')],by='taxon')
P_w_traits_full=P_w_traits_full[P_w_traits_full$taxon %in% trait_table$taxon,]
P_w_traits_full=P_w_traits_full[P_w_traits_full$taxon %in% c('Farula','Stactobiella')==F,]
P_w_traits_full=P_w_traits_full[P_w_traits_full$Regional_Response!='Neutral' & !is.na(P_w_traits_full$Regional_Response),]
P_Inc=P_w_traits_full[P_w_traits_full$Regional_Response=='Increaser',]
P_Dec=P_w_traits_full[P_w_traits_full$Regional_Response!='Increaser',]

P_Inc=P_Inc[,c(350,352:491)]
P_Dec=P_Dec[,c(351,353:491)]

P_everything=rbind(P_Inc,P_Dec)
P_everything[P_everything==T]<-'1'
P_everything[P_everything==F]<-'0'
names(P_everything)[names(P_everything)=='taxon']<-'OTU'

final_traits=names(trait_table[,2:ncol(trait_table)])
# identify categorical columns
categorical_cols <- names(P_everything)[names(P_everything) %in% names(trait_table)[2:ncol(trait_table)]]

summary_df <- as.data.frame(P_everything %>%
  group_by(
    OTU,
    Regional_Response
  ) %>%
  dplyr::summarise(

    # categorical: sum presences
   dplyr::across(
      all_of(categorical_cols),
      ~sum(. == 1, na.rm = TRUE)
    ),

    .groups = "drop"
  ))

summary_df

final_traits=names(trait_table[,2:ncol(trait_table)])
names(P_everything)
#-----------------------------------
# CONVERT TO LONG FORMAT
#-----------------------------------

long_df <- P_everything %>%
  pivot_longer(
    cols = all_of(c(final_traits)),
    names_to = "Trait",
    values_to = "Value"
  )



# find traits we cannot assess:
# weird_traits1=long_df %>%
#   dplyr::group_by(Trait, Regional_Response) %>%
#   dplyr::summarise(
#     n = sum(!is.na(Value)),
#     sd = sd(Value, na.rm = TRUE),
#     .groups = "drop"
#   )
# good_traits=weird_traits1[!(weird_traits1$n <2 | is.na(weird_traits1$sd)),]
# good_traits=unique(good_traits$Trait)


#Fisher

fisher_results <- long_df %>%

  dplyr::filter(Trait %in% categorical_cols) %>%
  dplyr::mutate(Present = Value == 1) %>%
  dplyr::group_by(Trait) %>%

  dplyr::summarise(

    PercentInc = {
      inc <- Present[Regional_Response == "Increaser"]
      inc <- inc[!is.na(inc)]
      if(length(inc) == 0) NA_real_ else mean(inc) * 100
    },

    PercentDec = {
      dec <- Present[Regional_Response == "Decreaser"]
      dec <- dec[!is.na(dec)]
      if(length(dec) == 0) NA_real_ else mean(dec) * 100
    },

    odds_ratio = {
      inc <- Present[Regional_Response == "Increaser"]
      dec <- Present[Regional_Response == "Decreaser"]

      inc <- inc[!is.na(inc)]
      dec <- dec[!is.na(dec)]

      if(length(inc) == 0 || length(dec) == 0) {
        NA_real_
      } else {

        tbl <- matrix(
          c(
            sum(inc), length(inc) - sum(inc),
            sum(dec), length(dec) - sum(dec)
          ),
          nrow = 2,
          byrow = TRUE
        )

        tryCatch(
          unname(fisher.test(tbl)$estimate),
          error = function(e) NA_real_
        )
      }
    },

    p_value = {
      inc <- Present[Regional_Response == "Increaser"]
      dec <- Present[Regional_Response == "Decreaser"]

      inc <- inc[!is.na(inc)]
      dec <- dec[!is.na(dec)]

      if(length(inc) == 0 || length(dec) == 0) {
        NA_real_
      } else {

        tbl <- matrix(
          c(
            sum(inc), length(inc) - sum(inc),
            sum(dec), length(dec) - sum(dec)
          ),
          nrow = 2,
          byrow = TRUE
        )

        tryCatch(
          fisher.test(tbl)$p.value,
          error = function(e) NA_real_
        )
      }
    },
    n_Increaser = {
      inc <- Present[Regional_Response == "Increaser"]
      sum(!is.na(inc))
    },

    n_Decreaser = {
      dec <- Present[Regional_Response == "Decreaser"]
      sum(!is.na(dec))
    },
    .groups = "drop"
  )


fisher_results2 <- long_df %>%

  dplyr::filter(Trait %in% categorical_cols) %>%
  dplyr::mutate(Present = Value == 1) %>%
  dplyr::group_by(Trait) %>%

  dplyr::summarise(

    PercentInc = {
      inc <- Present[Regional_Response == "Increaser"]
      inc <- inc[!is.na(inc)]

      if(length(inc) == 0) {
        NA_real_
      } else {
        mean(inc) * 100
      }
    },

    PercentDec = {
      dec <- Present[Regional_Response == "Decreaser"]
      dec <- dec[!is.na(dec)]

      if(length(dec) == 0) {
        NA_real_
      } else {
        mean(dec) * 100
      }
    },

    # Counts used in Fisher test
    inc_present = {
      inc <- Present[Regional_Response == "Increaser"]
      sum(inc, na.rm = TRUE)
    },

    inc_absent = {
      inc <- Present[Regional_Response == "Increaser"]
      sum(!inc, na.rm = TRUE)
    },

    dec_present = {
      dec <- Present[Regional_Response == "Decreaser"]
      sum(dec, na.rm = TRUE)
    },

    dec_absent = {
      dec <- Present[Regional_Response == "Decreaser"]
      sum(!dec, na.rm = TRUE)
    },

    # Sample sizes
    n_Increaser = {
      inc <- Present[Regional_Response == "Increaser"]
      sum(!is.na(inc))
    },

    n_Decreaser = {
      dec <- Present[Regional_Response == "Decreaser"]
      sum(!is.na(dec))
    },

    n_total = {
      inc <- Present[Regional_Response == "Increaser"]
      dec <- Present[Regional_Response == "Decreaser"]

      sum(!is.na(inc)) + sum(!is.na(dec))
    },

    odds_ratio = {
      inc <- Present[Regional_Response == "Increaser"]
      dec <- Present[Regional_Response == "Decreaser"]

      inc <- inc[!is.na(inc)]
      dec <- dec[!is.na(dec)]

      if(length(inc) == 0 || length(dec) == 0) {

        NA_real_

      } else {

        tbl <- matrix(
          c(
            sum(inc), length(inc) - sum(inc),
            sum(dec), length(dec) - sum(dec)
          ),
          nrow = 2,
          byrow = TRUE
        )

        tryCatch(
          unname(fisher.test(tbl)$estimate),
          error = function(e) NA_real_
        )
      }
    },

    p_value = {
      inc <- Present[Regional_Response == "Increaser"]
      dec <- Present[Regional_Response == "Decreaser"]

      inc <- inc[!is.na(inc)]
      dec <- dec[!is.na(dec)]

      if(length(inc) == 0 || length(dec) == 0) {

        NA_real_

      } else {

        tbl <- matrix(
          c(
            sum(inc), length(inc) - sum(inc),
            sum(dec), length(dec) - sum(dec)
          ),
          nrow = 2,
          byrow = TRUE
        )

        tryCatch(
          fisher.test(tbl)$p.value,
          error = function(e) NA_real_
        )
      }
    },

    .groups = "drop"
  )
clipr::write_clip(as.data.frame(fisher_results2))


site_traits= P_long2 %>% dplyr::left_join(long_df %>% dplyr::select(taxon,all_of(final_traits)),
                                                  by='taxon')


Prob_site_cols=names(P_everything)[1:338]
Prob_Trait_cols=names(P_everything)[340:375]#ncol(P_everything)]
Prob_site_cols=Prob_site_cols[Prob_site_cols]
long_sites <- P_everything %>%

  tidyr::pivot_longer(
    cols = all_of(Prob_site_cols),
    names_to = "Site",
    values_to = "Presence"
  )
long_sites <- long_sites %>%
  dplyr::filter(Presence == 1)

names(long_sites)

final_traits=names(trait_table)[2:ncol(trait_table)]

P_for_RF=P_everything[,names(P_everything) %in% c(names(P_everything)[1:349], final_traits,'OTU','Regional_Response')]

P_for_RF=P_for_RF[,names(P_for_RF) %in% c(names(P_everything)[1:349],final_traits,'OTU','Regional_Response')]
# P_for_RF=P_for_RF[,!(names(P_for_RF) %in% c('Develop_fast_season',
#                                           'Resp_abbrev_Gills',
#                                           'Occurence_drift_common',
#                                           'Swim_strong',
#                                           'Attach_free_range',
#                                           'Attach_both',
#                                           'Female_disp_abbrev_High','Swim_none',
#                                           'Attach_free_range_both',
#                                           'Emerge_season_1_Spring',
#                                           'Emerge_season_2_Summer'))]

#P_for_RF=P_for_RF[!duplicated(P_for_RF$OTU),]
O_w_traits=O_for_traits
# O_for_RF=O_w_Traits[,names(O_w_Traits) %in% c(names(O_w_Traits)[1:656], final_traits,'OTU','Regional_Response')]
# O_for_RF=O_for_RF[,names(O_for_RF) %in% c(names(O_for_RF)[1:656], final_traits,'OTU','Regional_Response')]

maxxed_OTUs2=trait_table
names(maxxed_OTUs2)[1]='taxon'
P_for_RF=plyr::join(maxxed_OTUs2,P_Results2[,c('taxon','Regional_Response')],by='taxon','left')

names(P_for_RF)
P_for_RF=P_for_RF[!duplicated(P_for_RF$taxon),]
P_for_RF=P_for_RF[which(!is.na(P_for_RF$Regional_Response) & P_for_RF$Regional_Response!='Neutral'),]
row.names(P_for_RF)=P_for_RF$taxon;P_for_RF=P_for_RF[,-1]
swapP_for_RF=P_for_RF

Psitecols=names(P_for_RF)[1:350]
Ptraitcols=names(P_for_RF)[352:436]

Psite_taxa=ProbOs
Osite_taxa=Os
maxxed_OTUs3=maxxed_OTUs2
maxxed_OTUs3=maxxed_OTUs3[!duplicated(maxxed_OTUs3$taxon),]
#row.names(maxxed_OTUs3)=maxxed_OTUs3$taxon;maxxed_OTUs3=maxxed_OTUs3[,-1]
common_taxa = intersect(colnames(Psite_taxa),
                        row.names(maxxed_OTUs3))


Psite_taxa  <- Psite_taxa[, common_taxa]
taxa_traits <- maxxed_OTUs3[common_taxa, ]


Psite_trait_counts <- as.matrix(Psite_taxa) %*%
  as.matrix(taxa_traits)
Psite_trait_pa=as.data.frame((Psite_trait_counts >0)*1)

Psite_rich=rowSums(Psite_taxa)

Psite_trait_pct=as.data.frame(sweep(Psite_trait_counts,1,Psite_rich,"/")*100)
Psite_trait_pct$sampleId=row.names(Psite_trait_pct)
Psite_trait_pa$sampleId=row.names(Psite_trait_pa)




#O_pres_RF=plyr::join(Osite_trait_matrix,Ref_preds,by='Site')
P_Pres_RF=plyr::join(Psite_trait_matrix,Prob_preds,by='Site')


#Same but for Ref sites.
maxxed_OTUs3=maxxed_OTUs2
maxxed_OTUs3=maxxed_OTUs3[!duplicated(maxxed_OTUs3$taxon),]
row.names(maxxed_OTUs3)=maxxed_OTUs3$taxon;maxxed_OTUs3=maxxed_OTUs3[,-1]
common_taxa = intersect(colnames(Osite_taxa),
                        row.names(maxxed_OTUs3))


Osite_taxa  <- Osite_taxa[, common_taxa]
taxa_traits <- maxxed_OTUs3[common_taxa, ]


Osite_trait_counts <- as.matrix(Osite_taxa) %*%
  as.matrix(taxa_traits)
Osite_trait_pa=as.data.frame((Osite_trait_counts >0)*1)

Osite_rich=rowSums(Osite_taxa)

Osite_trait_pct=as.data.frame(sweep(Osite_trait_counts,1,Osite_rich,"/")*100)
Osite_trait_pct$sampleId=row.names(Osite_trait_pct)
Osite_trait_pa$sampleId=row.names(Osite_trait_pa)



env_dat2=rbind(NMDS_cal_preds[,names(NMDS_prob_preds)],NMDS_prob_preds)
AgUrb_dat=read.csv('C://Users//andrew.caudillo.BUGLAB-I9//Box//NAMC//Research Projects//AIM//IncreaserDecreaser_OE//all_site_UrbAg.csv')

Psite_trait_pa=plyr::join(Psite_trait_pa,env_dat2,by='sampleId')
Psite_trait_pa=plyr::join(Psite_trait_pa,AgUrb_dat,by='COMID')

Psite_trait_pct=plyr::join(Psite_trait_pct,env_dat2,by='sampleId')
Psite_trait_pct=plyr::join(Psite_trait_pct,AgUrb_dat,by='COMID')

Osite_trait_pa=plyr::join(Osite_trait_pa,env_dat2,by='sampleId')
Osite_trait_pa=plyr::join(Osite_trait_pa,AgUrb_dat,by='COMID')

Osite_trait_pct=plyr::join(Osite_trait_pct,env_dat2,by='sampleId')
Osite_trait_pct=plyr::join(Osite_trait_pct,AgUrb_dat,by='COMID')

row.names(Psite_trait_pct)=Psite_trait_pct$sampleId

library(randomForestSRC)
library(randomForestSRC)  # for multivariate random forest

Psite_trait_PCtIndivs=read.csv('C://Users//andrew.caudillo.BUGLAB-I9//Box//NAMC//Research Projects//AIM//IncreaserDecreaser_OE//Updated_w_modelObj//New_for_SFS//SiteStatus_models//Prob_Pct_indivs_wtrait_and_envs.csv')
psitetrat
# ---------------------------
# set the random seed and isolate predictor names
set.seed(42)
#setting the index for where predictors start
pred_index=names(Psite_trait_PCtIndivs)[c(26:39,42,46)]
#pred_index<-names(Psite_trait_pct)[c(27:40,43,47)]

#force the predictors to be numbers, not strings, as SQL shows
library(randomForestSRC)

Psite_trait_pa=read.csv('C://Users//andrew.caudillo.BUGLAB-I9//Box//NAMC//Research Projects//AIM//IncreaserDecreaser_OE//Updated_w_modelObj//New_for_SFS//SiteStatus_models//PA_all_traits_probabilistic.csv')
Psite_trait_pa=Psite_trait_pa[Psite_trait_pa$sampleId %in% sites_w_no_bugs==F,]

train <- Psite_trait_pa
train=Psite_trait_PCtIndivs
trait_names=names(train)[1:25]
#only for pct MRF
#variables that the model finds as the top variables via the
#max.subtree function.
# vars=c(  "Precip8110",               "CaOWs",             "NCat",
#           "RunoffWs",          "SWs",               "PctSalLake",               "MAST_mean08091314",
#            "MSST_mean08091314", "MWST_mean08091314",'DOY')
topvarsPA=c("WsAreaSqKm",        "Precip8110",        "MSST_mean08091314", "ElevCat",
            "Tmean8110Ws"  )
topvarsPCT=c("PctSalLake" ,       "MSST_mean08091314", "Precip8110",        "WsAreaSqKm",
             "SWs"   )
topvarsIndv=c("WsAreaSqKm",        "Precip8110",        "MSST_mean08091314", "ElevCat",
              "Tmean8110Ws" )
vars=pred_index
formula_string <- paste0(
  "cbind(", paste(trait_names, collapse = ", "), ") ~ ",
  #paste(pred_names,collapse='+'))
  #paste(names(train[,(ncol(train)-13):ncol(train)]),collapse = '+'))
  paste(names(train[,vars]),collapse = '+'))

# Convert string to formula using parse + eval
#Else MRF will not read it correctly
rf_formula <- eval(parse(text = formula_string))
rf_formula <- as.formula(formula_string)

#only if using P/A of trait states
train[trait_names]<-lapply(train[trait_names],as.factor)

#this model has already been fit with the best hyperparameters
#as found by the tune() function

tune(rf_formula,
     train,
     doBest = T,
     improve=1e-1, #define what improvement means for OOB error, larger values
     #= faster run-times, but less improvement.
     #sampsize = 30, #number of samples for the subsamples, #default from rfsrc doc
     trace=T, #show the progress in console
     nodesizeTry = c(5, 10,20,30), #various terinal node sizes. Bigger node
     #means better fit, but higher chance of overfit
     ntreeTry = 500,
     nsplit=5)

#% traits
rf_model <- rfsrc(
  formula=rf_formula,
  data=train,
  ntree = 500,
  mtry=14,
  importance = 'permute',
  nodedepth = 5,
  nodesize=10,
  nsplit=5)
#For P/A
rf_model <- rfsrc(
  formula=rf_formula,
  data=train,
  ntree = 500,
  mtry=13,
  importance = 'permute',
  nodedepth = 20,
  nodesize=20,
  nsplit=5
  #case.wt = case_wts,
  #na.action='na.impute' #replace any NAs in predictors with median for that column
)
#indivs
rf_model <- rfsrc(
  formula=rf_formula,
  data=train,
  ntree = 500,
  mtry=4,
  importance = 'permute',
  nodedepth = 20,
  nodesize=20,
  nsplit=5
  #case.wt = case_wts,
  #na.action='na.impute' #replace any NAs in predictors with median for that column
)
#view OOB errors for each response



OOB_pred_each=get.mv.error(rf_model)
mean(OOB_pred_each)
mean(OOB_pred_each / 100^2)
#on average, how erroneous is each taxon's predictions?
sqrt(mean(OOB_pred_each))
#for Count
Y=train[,trait_names]
pseudo_r2 <- sapply(names(rf_model$regrOutput), function(v) {

  observed <- Y[[v]]

  predicted <- rf_model$regrOutput[[v]]$predicted

  mse <- mean((observed - predicted)^2,
              na.rm = TRUE)

  1 - mse / var(observed, na.rm = TRUE)
})

mean(pseudo_r2,na.rm=T)









response_vars=trait_names
pred_vars=pred_index

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


top_vars

vimp_sorted


vimp_summary <- data.frame(
  Predictor = row.names(vimp_df),
  Mean = apply(vimp_df, 1, mean, na.rm = TRUE),
  Median = apply(vimp_df, 1, median, na.rm = TRUE),
  Max = apply(vimp_df, 1, max, na.rm = TRUE),
  Consistency = apply(vimp_df, 1, function(x) mean(x > quantile(x, 0.9), na.rm = TRUE))
)

library(ggplot2)
ggplot(vimp_summary, aes(x = reorder(Predictor, Mean), y = Mean)) +
  geom_col(fill = "dodgerblue4", alpha = 0.7) +
  coord_flip() +
  labs(x = "Predictor", y = "Mean Variable Importance\nAcross Traits")+
       #title = "Mean Variable Importance Across Taxa") +
  theme_minimal(base_size = 11)+
  theme_classic()+
    # scale_x_discrete(labels=c('Catchment Elevation',
    #                            '30-year Mean Annual Precip',
    #                            'Mean Summer Stream Temp',
    #                             '30-year Mean Annual Temp',
    #                            'Basin Area (Sq. km)'))+
  theme(axis.text = element_text(size=20),
        axis.title = element_text(size=20))


savp(10,8,'C://Users//andrew.caudillo.BUGLAB-I9//Box//NAMC//Research Projects//AIM//IncreaserDecreaser_OE//Updated_w_modelObj//New_for_SFS//I_D_PA_Mod_vimp_full_260518.png')

#VAR IMP

for(i in 1:nrow(vimp_summary)){
  partial_data <- ggRandomForests::gg_partial(plot.variable(rf_model,xvar.names=vimp_summary$Predictor[i], partial = TRUE,smooth.lines=T))
  if(vimp_summary$Predictor[i]=='PctSalLake'){
  # Plot and adjust fonts
  partial_data$continuous=partial_data$categorical
  partial_data$continuous$x=as.numeric(as.character(partial_data$continuous$x))
  #if PctSalLake
  df <- partial_data$continuous
  p <- ggplot(df, aes(x = x, y = yhat)) +
    labs(y='Partial Effect')+
    geom_line(linewidth = 2) +
    theme_classic(base_size = 18) +
    labs(title='PctSalLake')+
    theme(
      axis.title = element_text(size = 20),
      axis.text = element_text(size = 20),
      strip.text = element_text(size = 24),
      strip.background = element_blank(),
      plot.title = element_text(size=20,hjust=0.5),
      legend.text = element_text(size=20)
    )

  print(p)
  savp(10,8,paste0('C://Users//andrew.caudillo.BUGLAB-I9//Box//NAMC//Research Projects//AIM//IncreaserDecreaser_OE//Updated_w_modelObj//New_for_SFS//I_D_PAMod_PDP_full_2605018',
              vimp_summary$Predictor[i],'.png'))
  }else{
  p=plot(partial_data)
 p=p+geom_line(linewidth=2)+
    theme_classic(base_size=18)+
    theme(text = element_text(size = 20),
          axis.title = element_text(size = 20),
          axis.text=element_text(size=20),
          strip.text = element_text(size = 24),
          strip.background = element_blank()
    )
 print(p)
 savp(10,8,paste0('C://Users//andrew.caudillo.BUGLAB-I9//Box//NAMC//Research Projects//AIM//IncreaserDecreaser_OE//Updated_w_modelObj//New_for_SFS//I_D_PAMod_PDP_full_260518',
             vimp_summary$Predictor[i],'.png'))
  }

}

#PA
windows(10,8)
ggplot(vimp_summary, aes(x = reorder(Predictor, Mean), y = Mean)) +
  geom_col(fill = "dodgerblue4", alpha = 0.7) +
  coord_flip() +
  labs(x = "Predictor", y = "Mean Variable Importance\n Across Traits")+
       #title = "Mean Variable Importance Across Taxa") +
  #theme_minimal(base_size = 11)+
  theme_classic()+
    scale_x_discrete(labels=c(
                              'Catchment Elevation',
                              '30-year Annual Temp',
                              'Mean Summer Stream Temp',
                              '30-year Mean Annual Precip',
                             'Basin Area (Sq. km)'))+
   theme(axis.text = element_text(size=12),
        axis.title = element_text(size=15))

savp(10,8,'C://Users//andrew.caudillo.BUGLAB-I9//Box//NAMC//Research Projects//AIM//IncreaserDecreaser_OE//Updated_w_modelObj//New_for_SFS//I_D_PAMod_vimp_Trimmed.png')


#VAR IMP

for(i in 1:nrow(vimp_summary)){
  partial_data <- ggRandomForests::gg_partial(plot.variable(rf_model,xvar.names=vimp_summary$Predictor[i], partial = TRUE,smooth.lines=T))
  if(vimp_summary$Predictor[i]=='PctSalLake'){
    # Plot and adjust fonts
    partial_data$continuous=partial_data$categorical
    partial_data$continuous$x=as.numeric(as.character(partial_data$continuous$x))
    #if PctSalLake
    df <- partial_data$continuous
    p <- ggplot(df, aes(x = x, y = yhat)) +
      labs(y='Partial Effect')+
      geom_line(linewidth = 2) +
      theme_classic(base_size = 18) +
      labs(title='PctSalLake')+
      theme(
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 20),
        strip.text = element_text(size = 24),
        strip.background = element_blank(),
        plot.title = element_text(size=20,hjust=0.5),
        legend.text = element_text(size=20)
      )

    print(p)
    savp(10,8,paste0('C://Users//andrew.caudillo.BUGLAB-I9//Box//NAMC//Research Projects//AIM//IncreaserDecreaser_OE//Updated_w_modelObj//New_for_SFS//I_D_PAMod_PDP_trimmed',
                     vimp_summary$Predictor[i],'.png'))
  }else{
    p=plot(partial_data)
    p=p+geom_line(linewidth=2)+
      theme_classic(base_size=18)+
      theme(text = element_text(size = 20),
            axis.title = element_text(size = 20),
            axis.text=element_text(size=20),
            strip.text = element_text(size = 24),
            strip.background = element_blank()
      )
    print(p)
    savp(10,8,paste0('C://Users//andrew.caudillo.BUGLAB-I9//Box//NAMC//Research Projects//AIM//IncreaserDecreaser_OE//Updated_w_modelObj//New_for_SFS//I_D_PAMod_PDP_trimmed',
                     vimp_summary$Predictor[i],'.png'))
  }

}

#Pct Indivs

windows(10,8)
ggplot(vimp_summary, aes(x = reorder(Predictor, Mean), y = Mean)) +
  geom_col(fill = "dodgerblue4", alpha = 0.7) +
  coord_flip() +
  labs(x = "Predictor", y = "Mean Variable Importance\n Across Traits")+
  #title = "Mean Variable Importance Across Taxa") +
  #theme_minimal(base_size = 11)+
  theme_classic()+
   scale_x_discrete(labels=c(
    'Catchment Elevation',
     '30-year Annual Temp',
     'Mean Summer Stream Temp',
     '30-year Mean Annual Precip',
     'Basin Area (Sq. km)'))+
  theme(axis.text = element_text(size=12),
        axis.title = element_text(size=15))

savp(10,8,'C://Users//andrew.caudillo.BUGLAB-I9//Box//NAMC//Research Projects//AIM//IncreaserDecreaser_OE//Updated_w_modelObj//New_for_SFS//I_D_PctIndiv_vimp_Trimmed.png')


#VAR IMP

for(i in 1:nrow(vimp_summary)){
  partial_data <- ggRandomForests::gg_partial(plot.variable(rf_model,xvar.names=vimp_summary$Predictor[i], partial = TRUE,smooth.lines=T))
  if(vimp_summary$Predictor[i]=='PctSalLake'){
    # Plot and adjust fonts
    partial_data$continuous=partial_data$categorical
    partial_data$continuous$x=as.numeric(as.character(partial_data$continuous$x))
    #if PctSalLake
    df <- partial_data$continuous
    p <- ggplot(df, aes(x = x, y = yhat)) +
      labs(y='Partial Effect')+
      geom_line(linewidth = 2) +
      theme_classic(base_size = 18) +
      labs(title='PctSalLake')+
      theme(
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 20),
        strip.text = element_text(size = 24),
        strip.background = element_blank(),
        plot.title = element_text(size=20,hjust=0.5),
        legend.text = element_text(size=20)
      )

    print(p)
    savp(10,8,paste0('C://Users//andrew.caudillo.BUGLAB-I9//Box//NAMC//Research Projects//AIM//IncreaserDecreaser_OE//Updated_w_modelObj//New_for_SFS//I_D_PCtIndivsMod_PDP_trimmed',
                     vimp_summary$Predictor[i],'.png'))
  }else{
    p=plot(partial_data)
    p=p+geom_line(linewidth=2)+
      theme_classic(base_size=18)+
      theme(text = element_text(size = 20),
            axis.title = element_text(size = 20),
            axis.text=element_text(size=20),
            strip.text = element_text(size = 24),
            strip.background = element_blank()
      )
    print(p)
    savp(10,8,paste0('C://Users//andrew.caudillo.BUGLAB-I9//Box//NAMC//Research Projects//AIM//IncreaserDecreaser_OE//Updated_w_modelObj//New_for_SFS//I_D_PctIndivs_PDP_trimmed',
                     vimp_summary$Predictor[i],'.png'))
  }

}




## same as above, but for reference sites###############

library(randomForestSRC)
library(randomForestSRC)  # for multivariate random forest

# ---------------------------
# set the random seed and isolate predictor names
set.seed(42)
#setting the index for where predictors start
pred_index2<-names(Osite_trait_pct)[c(27:40,43,47)]

#force the predictors to be numbers, not strings, as SQL shows
train <- Osite_trait_pa
trait_names=names(train)[1:25]
#only for pct MRF
train=train[complete.cases(train[,trait_names]),]

#assign the sampleId to row names, then drop it
#the names of the taxa
constant_cols <- names(train)[
  sapply(train, function(x)
    length(unique(x[!is.na(x)])) <= 1)
]

constant_cols
#only for MRF
trait_names=trait_names[trait_names %in% constant_cols==F]


response_var <- sapply(
  train[trait_names],
  function(x) var(x, na.rm = TRUE)
)

sort(response_var)
#variables that the model finds as the top variables via the
#max.subtree function.
# vars=c(  "Precip8110",               "CaOWs",             "NCat",
#           "RunoffWs",          "SWs",               "PctSalLake",               "MAST_mean08091314",
#            "MSST_mean08091314", "MWST_mean08091314",'DOY')
topvarsPA_REF=c("WsAreaSqKm",        "Precip8110",        "MSST_mean08091314", "LAT",
                "ElevCat"   )
topvarsPCT_REF=c("MSST_mean08091314", "LAT" ,              "Tmean8110Ws",       "LONG",
                 "MAST_mean08091314")
vars=topvarsPA_REF

formula_string <- paste0(
  "cbind(", paste(trait_names, collapse = ", "), ") ~ ",
  #paste(pred_names,collapse='+'))
  #paste(names(train[,(ncol(train)-13):ncol(train)]),collapse = '+'))
  paste(names(train[,vars]),collapse = '+'))

# Convert string to formula using parse + eval
#Else MRF will not read it correctly
rf_formula <- eval(parse(text = formula_string))
rf_formula <- as.formula(formula_string)

#only if using P/A of trait states
train[trait_names]<-lapply(train[trait_names],as.factor)

#this model has already been fit with the best hyperparameters
#as found by the tune() function
row.names(train)=train$sampleId
#only for pct:
train=train[row.names(train) %in% c('HAWK-33','WE-131','WE-209')==F,]
tune(rf_formula,
     train,
     doBest = T,
     improve=1e-3, #define what improvement means for OOB error, larger values
     #= faster run-times, but less improvement.
     #sampsize = 30, #number of samples for the subsamples, #default from rfsrc doc
     trace=T, #show the progress in console
     nodesizeTry = c(5, 10,20,30), #various terinal node sizes. Bigger node
     #means better fit, but higher chance of overfit
     ntreeTry = 500,
     nsplit=5)

#% traits
rf_model <- rfsrc(
  formula=rf_formula,
  data=train,
  ntree = 500,
  mtry=4,
  importance = 'permute',
  nodedepth = 5,
  nodesize=10,
  nsplit=5)
#For P/A
rf_model <- rfsrc(
  formula=rf_formula,
  data=train,
  ntree = 500,
  mtry=4,
  importance = 'permute',
  nodedepth = 5,
  nodesize=20,
  nsplit=5
  #case.wt = case_wts,
  #na.action='na.impute' #replace any NAs in predictors with median for that column
)
#view OOB errors for each response
OOB_pred_each=get.mv.error(rf_model)
#on average, how erroneous is each taxon's predictions?
mean(OOB_pred_each)
sqrt(mean(OOB_pred_each))
#for Count
Y=train[,trait_names]
pseudo_r2 <- sapply(names(rf_model$regrOutput), function(v) {

  observed <- Y[[v]]

  predicted <- rf_model$regrOutput[[v]]$predicted

  mse <- mean((observed - predicted)^2,
              na.rm = TRUE)

  1 - mse / var(observed, na.rm = TRUE)
})

mean(pseudo_r2,na.rm=T)









response_vars=trait_names
pred_vars=topvarsPCT_REF

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


top_vars

vimp_sorted


vimp_summary <- data.frame(
  Predictor = row.names(vimp_df),
  Mean = apply(vimp_df, 1, mean, na.rm = TRUE),
  Median = apply(vimp_df, 1, median, na.rm = TRUE),
  Max = apply(vimp_df, 1, max, na.rm = TRUE),
  Consistency = apply(vimp_df, 1, function(x) mean(x > quantile(x, 0.9), na.rm = TRUE))
)

library(ggplot2)
ggplot(vimp_summary, aes(x = reorder(Predictor, Mean), y = Mean)) +
  geom_col(fill = "dodgerblue4", alpha = 0.7) +
  coord_flip() +
  labs(x = "Predictor", y = "Mean Variable Importance\nAcross Traits")+
  #title = "Mean Variable Importance Across Taxa") +
  theme_minimal(base_size = 11)+
  theme_classic()+
      scale_x_discrete(labels=c('Mean Summer Stream Temp',
                                  'Latitude',
                                     'Longitude',
                              '30-year Mean Temp',
                               'Mean Annual Stream Temp'))+
  theme(axis.text = element_text(size=14),
        axis.title = element_text(size=15))


savp(10,8,'C://Users//andrew.caudillo.BUGLAB-I9//Box//NAMC//Research Projects//AIM//IncreaserDecreaser_OE//Updated_w_modelObj//New_for_SFS//Reference_I_D_Pct_Mod_vimp_trimmed.png')
#PA
windows(10,8)
ggplot(vimp_summary, aes(x = reorder(Predictor, Mean), y = Mean)) +
  geom_col(fill = "dodgerblue4", alpha = 0.7) +
  coord_flip() +
  labs(x = "Predictor", y = "Mean Variable Importance\n Across Traits")+
  #title = "Mean Variable Importance Across Taxa") +
  #theme_minimal(base_size = 11)+
  theme_classic()+
  scale_x_discrete(labels=c(
    'Mean Summer Stream Temp',
    'Basin Area (Sq. km)',
    'Longitude',
    'Mean Winter Stream Temp',
    '30-year Mean Annual Precip'))+
  theme(axis.text = element_text(size=12),
        axis.title = element_text(size=15))

savp(10,8,'C://Users//andrew.caudillo.BUGLAB-I9//Box//NAMC//Research Projects//AIM//IncreaserDecreaser_OE//Updated_w_modelObj//New_for_SFS//I_D_PAMod_vimp_trimmed.png')

#PDPs
for(i in 1:nrow(vimp_summary)){
  partial_data <- ggRandomForests::gg_partial(plot.variable(rf_model,xvar.names=vimp_summary$Predictor[i], partial = TRUE,smooth.lines=T))
  if(vimp_summary$Predictor[i]=='PctSalLake'){
    # Plot and adjust fonts
    partial_data$continuous=partial_data$categorical
    partial_data$continuous$x=as.numeric(as.character(partial_data$continuous$x))
    #if PctSalLake
    df <- partial_data$continuous
    p <- ggplot(df, aes(x = x, y = yhat)) +
      labs(y='Partial Effect')+
      geom_line(linewidth = 2) +
      theme_classic(base_size = 18) +
      labs(title='PctSalLake')+
      theme(
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 20),
        strip.text = element_text(size = 24),
        strip.background = element_blank(),
        plot.title = element_text(size=20,hjust=0.5),
        legend.text = element_text(size=20)
      )

    print(p)
    savp(10,8,paste0('C://Users//andrew.caudillo.BUGLAB-I9//Box//NAMC//Research Projects//AIM//IncreaserDecreaser_OE//Updated_w_modelObj//New_for_SFS//Ref_PctMod_PDP_trimmed',
                     vimp_summary$Predictor[i],'.png'))
  }else{
    p=plot(partial_data)
    p=p+geom_line(linewidth=2)+
      theme_classic(base_size=18)+
      theme(text = element_text(size = 20),
            axis.title = element_text(size = 20),
            axis.text=element_text(size=20),
            strip.text = element_text(size = 24),
            strip.background = element_blank()
      )
    print(p)
    savp(10,8,paste0('C://Users//andrew.caudillo.BUGLAB-I9//Box//NAMC//Research Projects//AIM//IncreaserDecreaser_OE//Updated_w_modelObj//New_for_SFS//REF_PctMod_PDP_trimmed',
                     vimp_summary$Predictor[i],'.png'))
  }

}










if(0){
  library(randomForestSRC)  # for multivariate random forest

  df_work=O_w_Traits[,658:ncol(O_w_Traits)]

  repeat {

    row_na <- rowMeans(is.na(df_work))
    col_na <- colMeans(is.na(df_work))

    if(max(row_na) == 0 && max(col_na) == 0) break

    if(max(row_na) > max(col_na)) {
      df_work <- df_work[row_na < max(row_na), ]
    } else {
      df_work <- df_work[, col_na < max(col_na)]
    }
  }
  full=merge(df_work,Ptrait_sub[,'OTU'],by='row.names')
  clipr::write_clip(full)


  P_everything
}


PA_traits_ForNMDS1=Psite_trait_pa
PA_traits_ForNMDS1$status='Probabilistic'
row.names(PA_traits_ForNMDS1)=row.names(Psite_trait_counts)
PA_traits_ForNMDS2=Osite_trait_pa
row.names(PA_traits_ForNMDS2)=row.names(Osite_trait_counts)
PA_traits_ForNMDS2$status='Reference'



PA_traits_ForNMDS_all=rbind(PA_traits_ForNMDS1,PA_traits_ForNMDS2)
PA_traits_ForNMDS_all=PA_traits_ForNMDS_all[,c(1:25,48)]
PA_traits_ForNMDS_all =PA_traits_ForNMDS_all[rowSums(PA_traits_ForNMDS_all[,-ncol(PA_traits_ForNMDS_all)])>0,]

PA_traits_ForNMDS_all=PA_traits_ForNMDS_all[match(env_dat3$sampleId, row.names(PA_traits_ForNMDS_all)),]
NMDS_PA_sor=vegan::vegdist(PA_traits_ForNMDS_all[,1:25],method='bray',binary=T)
set.seed(99)
NMDS_PA_MDS=vegan::metaMDS(NMDS_PA_sor)

row.names(PA_traits_ForNMDS_all)==env_dat3$sampleId
set.seed(99)
NMDS_PA_fit=envfit(NMDS_PA_MDS,scale(env_dat3[,!names(env_dat3)%in% c('COMID','sampleId')]))
NMDS_PA_scores=as.data.frame(scores(NMDS_PA_MDS))
NMDS_PA_scores$status=PA_traits_ForNMDS_all$status
env_vect<-as.data.frame(NMDS_PA_fit$vectors$arrows)
env_vec=as.data.frame(scores(NMDS_PA_fit, display = 'vectors'))
env_vect$env <- rownames(env_vect)
colnames(env_vect)[1:2] <- c("NMDS1", "NMDS2")

#PA_scores2=PA_scores[PA_scores$NMDS1<35 & PA_scores$NMDS2 < 2,]

arrow_multiplier <- 0.5 * max(abs(NMDS_PA_scores[,c(1,2)])) / max(abs(env_vect[,1:2]))
env_vect[,1:2] <- env_vect[,1:2] * arrow_multiplier
env_vect$env=c('Catchment Elev',
               'Mean Precip',
               'Mean Temp','Basin Area',
               '% Crops',
               '% Urban')

env_vect <- env_vect %>%
  mutate(
    nudge_x = case_when(
      env == "Basin Area"    ~  0.1,
      env == "Catchment Elev"~ 0.1,
      env == "% Urban"       ~  0.1,
      env == "Mean Precip"   ~0,
      env== 'Mean Temp' ~ -0.4,
      TRUE                   ~  0
    ),

    nudge_y = case_when(
      env=='% Urban' ~ 0.2,
      env=='% Mean Temp' ~ 0.5,
      env == "Catchment Elev"~ 0.2,
      env == "% Crops" ~ 0.028,
      env=='Mean Precip'~-0.05,
      TRUE             ~ 0
    )
  )
ggplot(NMDS_PA_scores,aes(x=NMDS1,y=NMDS2,fill=status))+geom_point(pch=21,size=3)+
  scale_fill_manual(name='Status',
                    values=c('Reference' = 'blue',
                             'Probabilistic' = 'orange3'))+
  geom_segment(data=env_vect,
               aes(x=0,y=0,xend=NMDS1,yend=NMDS2),#                arrow=arrow(length=unit(0.25,'cm')),
               linewidth=1,
               alpha=0.9,inherit.aes = F, arrow=arrow(length=unit(.25,"centimeters")))+
   # geom_text_repel(
   #   data = env_vect,
   #   aes(
   #     x = NMDS1,
   #     y = NMDS2,
   #     label = env
   #   ),
   #   nudge_x = env_vect$nudge_x,
   #   nudge_y = env_vect$nudge_y,
   #   size = 6,
   #   inherit.aes = FALSE,
   #   box.padding = 0.2,
   #   point.padding = 0.1,
   #   min.segment.length = 0,
   #   max.overlaps = Inf
   # )+
  stat_ellipse(level=0.8,data = subset(NMDS_PA_scores, status == "Reference"), color = "blue", size = 1,type='t') +
  stat_ellipse(level=0.8,data = subset(NMDS_PA_scores, status == "Probabilistic"), color = "orange3", size = 1)+
  theme_classic()+
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.title = element_text(size = 25),
    axis.title.x = element_text(margin=margin(t=20)),
    axis.title.y  = element_text(margin=margin(r=20)),
    strip.text = element_text(size = 14),
    legend.text = element_text(size = 20),
    legend.title = element_text(size = 20),
    legend.key.size = unit(1, "lines"),
    axis.ticks=element_blank(),
    legend.position = 'inside',
    legend.position.inside = c(0.88,0.15)
  )#+stat_ellipse(level=0.8,type='t')
savp(10,8,'C://Users//andrew.caudillo.BUGLAB-I9//Box//NAMC//Research Projects//AIM//IncreaserDecreaser_OE//Updated_w_modelObj//New_for_SFS//PA_trait_NMDS_260515_nolabs.png')


####
set.seed(99)
NMDS_PA_fit2=envfit(NMDS_PA_MDS,PA_traits_ForNMDS_all[,1:25])
NMDS_PA_scores2=as.data.frame(scores(NMDS_PA_MDS))
NMDS_PA_scores2$status=PA_traits_ForNMDS_all$status
env_vect<-as.data.frame(NMDS_PA_fit2$vectors$arrows)
env_vec=as.data.frame(scores(NMDS_PA_fit2, display = 'vectors'))
env_vect$env <- rownames(env_vect)
colnames(env_vect)[1:2] <- c("NMDS1", "NMDS2")

#PA_scores2=PA_scores[PA_scores$NMDS1<35 & PA_scores$NMDS2 < 2,]

arrow_multiplier <- 0.5 * max(abs(NMDS_PA_scores[,c(1,2)])) / max(abs(env_vect[,1:2]))
env_vect[,1:2] <- env_vect[,1:2] * arrow_multiplier
env_vect$env=c('Catchment Elev',
               'Mean Precip',
               'Mean Temp','Basin Area',
               '% Crops',
               '% Urban')

env_vect <- env_vect %>%
  mutate(
    nudge_x = case_when(
      env == "Basin Area"    ~  0.1,
      env == "Catchment Elev"~ 0.1,
      env == "% Urban"       ~  0.1,
      env == "Mean Precip"   ~0,
      env== 'Mean Temp' ~ -0.4,
      TRUE                   ~  0
    ),

    nudge_y = case_when(
      env=='% Urban' ~ 0.2,
      env=='% Mean Temp' ~ 0.5,
      env == "Catchment Elev"~ 0.2,
      env == "% Crops" ~ 0.028,
      env=='Mean Precip'~-0.05,
      TRUE             ~ 0
    )
  )
ggplot(NMDS_PA_scores,aes(x=NMDS1,y=NMDS2,fill=status))+geom_point(pch=21,size=3)+
  scale_fill_manual(name='Status',
                    values=c('Reference' = 'blue',
                             'Probabilistic' = 'orange3'))+
  geom_segment(data=env_vect,
               aes(x=0,y=0,xend=NMDS1,yend=NMDS2),#                arrow=arrow(length=unit(0.25,'cm')),
               linewidth=1,
               alpha=0.9,inherit.aes = F, arrow=arrow(length=unit(.25,"centimeters")))+
  # geom_text_repel(
  #   data = env_vect,
  #   aes(
  #     x = NMDS1,
  #     y = NMDS2,
  #     label = env
  #   ),
  #   nudge_x = env_vect$nudge_x,
  #   nudge_y = env_vect$nudge_y,
  #   size = 6,
  #   inherit.aes = FALSE,
  #   box.padding = 0.2,
  #   point.padding = 0.1,
  #   min.segment.length = 0,
  #   max.overlaps = Inf
  # )+
  stat_ellipse(level=0.8,data = subset(NMDS_PA_scores, status == "Reference"), color = "blue", size = 1,type='t') +
  stat_ellipse(level=0.8,data = subset(NMDS_PA_scores, status == "Probabilistic"), color = "orange3", size = 1)+
  theme_classic()+
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.title = element_text(size = 25),
    axis.title.x = element_text(margin=margin(t=20)),
    axis.title.y  = element_text(margin=margin(r=20)),
    strip.text = element_text(size = 14),
    legend.text = element_text(size = 20),
    legend.title = element_text(size = 20),
    legend.key.size = unit(1, "lines"),
    axis.ticks=element_blank(),
    legend.position = 'inside',
    legend.position.inside = c(0.88,0.15)
  )#+stat_ellipse(level=0.8,type='t')
savp(10,8,'C://Users//andrew.caudillo.BUGLAB-I9//Box//NAMC//Research Projects//AIM//IncreaserDecreaser_OE//Updated_w_modelObj//New_for_SFS//PA_trait_NMDS_260515_nolabs.png')




PCT_traits_ForNMDS1=Psite_trait_pct
PCT_traits_ForNMDS1$status='Probabilistic'
PCT_traits_ForNMDS1=PCT_traits_ForNMDS1[,names(PCT_traits_ForNMDS1)!='sampleId']
#row.names(PCT_traits_ForNMDS1)=#row.names(Psite_trait_counts)
PCT_traits_ForNMDS2=Osite_trait_pct
row.names(PCT_traits_ForNMDS2)=row.names(Osite_trait_counts)
PCT_traits_ForNMDS2$status='Reference'



PCT_traits_ForNMDS_all=rbind(PCT_traits_ForNMDS1,PCT_traits_ForNMDS2)
#PCT_traits_ForNMDS_all=PCT_traits_ForNMDS_all[,c(1:37,60)]
PCT_traits_ForNMDS_all=PCT_traits_ForNMDS_all[,names(PCT_traits_ForNMDS_all)%in% c('status',final_traits)]
PCT_traits_ForNMDS_all =PCT_traits_ForNMDS_all[rowSums(PCT_traits_ForNMDS_all[,-ncol(PCT_traits_ForNMDS_all)])>0,]

PCT_traits_ForNMDS_all=PCT_traits_ForNMDS_all[match(env_dat3$sampleId, row.names(PCT_traits_ForNMDS_all)),]
NMDS_PCT_sor=vegan::vegdist(PCT_traits_ForNMDS_all[,1:25],method='bray',binary=F)
set.seed(99)
NMDS_PCT_MDS=vegan::metaMDS(NMDS_PCT_sor)
set.seed(99)
NMDS_PCTfit=envfit(NMDS_PCT_MDS,scale(env_dat3[,!names(env_dat3)%in% c('COMID','sampleId')]))
NMDS_PCTscores=as.data.frame(scores(NMDS_PCT_MDS))
NMDS_PCTscores$status=PCT_traits_ForNMDS_all$status
env_vect<-as.data.frame(NMDS_PCTfit$vectors$arrows)
env_vec=as.data.frame(scores(NMDS_PCTfit, display = 'vectors'))
env_vect$env <- rownames(env_vect)
colnames(env_vect)[1:2] <- c("NMDS1", "NMDS2")

#PCTscores2=PCTscores[PCTscores$NMDS1<35 & PCTscores$NMDS2 < 2,]

arrow_multiplier <- 0.5 * max(abs(NMDS_PCTscores[,c(1,2)])) / max(abs(env_vect[,1:2]))
env_vect[,1:2] <- env_vect[,1:2] * arrow_multiplier
env_vect$env=c('Catchment Elev',
               'Mean Precip',
               'Mean Temp','Basin Area',
               '% Crops',
               '% Urban')

env_vect <- env_vect %>%
  mutate(
    nudge_x = case_when(
      env == "Basin Area"    ~  0.1,
      env == "Catchment Elev"~ 0.1,
      env == "% Urban"       ~  0.1,
      env == "Mean Precip"   ~0,
      env== 'Mean Temp' ~ -0.4,
      TRUE                   ~  0
    ),

    nudge_y = case_when(
      env=='% Urban' ~ 0.2,
      env=='% Mean Temp' ~ 0.5,
      env == "Catchment Elev"~ 0.2,
      env == "% Crops" ~ 0.028,
      env=='Mean Precip'~-0.05,
      TRUE             ~ 0
    )
  )
ggplot(NMDS_PCTscores,aes(x=NMDS1,y=NMDS2,fill=status))+geom_point(pch=21,size=3)+
  scale_fill_manual(name='Status',
                    values=c('Reference' = 'blue',
                             'Probabilistic' = 'orange3'))+
  geom_segment(data=env_vect,
               aes(x=0,y=0,xend=NMDS1,yend=NMDS2),#                arrow=arrow(length=unit(0.25,'cm')),
               linewidth=1,
               alpha=0.9,inherit.aes = F, arrow=arrow(length=unit(.25,"centimeters")))+
   geom_text_repel(
     data = env_vect,
     aes(
       x = NMDS1,
       y = NMDS2,
       label = env
     ),
     nudge_x = env_vect$nudge_x,
     nudge_y = env_vect$nudge_y,
     size = 6,
     inherit.aes = FALSE,
     box.padding = 0.2,
     point.padding = 0.1,
     min.segment.length = 0,
     max.overlaps = Inf
   )+
  stat_ellipse(level=0.8,data = subset(NMDS_PCTscores, status == "Reference"), color = "blue", size = 1,type='t') +
  stat_ellipse(level=0.8,data = subset(NMDS_PCTscores, status == "Probabilistic"), color = "orange3", size = 1)+
  theme_classic()+
  theme(
    axis.text.x = element_blank(),
    axis.text.y=element_blank(),
    axis.title = element_text(size = 25),
    axis.title.x = element_text(margin=margin(t=20)),
    axis.title.y  = element_text(margin=margin(r=20)),
    strip.text = element_text(size = 14),
    legend.text = element_text(size = 20),
    legend.title = element_text(size = 20),
    legend.key.size = unit(1, "lines"),
    axis.ticks=element_blank(),
    legend.position = 'inside',
    legend.position.inside = c(.15,0.15)
  )#+stat_ellipse(level=0.8,type='t')
savp(10,8,'C://Users//andrew.caudillo.BUGLAB-I9//Box//NAMC//Research Projects//AIM//IncreaserDecreaser_OE//Updated_w_modelObj//New_for_SFS//PCT_trait_NMDS_nolabs_260514.png')
Count_traits_ForNMDS1=as.data.frame(Psite_trait_counts)
Count_traits_ForNMDS1$status='Probabilistic'
row.names(Count_traits_ForNMDS1)=row.names(Psite_trait_counts)
Count_traits_ForNMDS2=as.data.frame(Osite_trait_counts)
row.names(Count_traits_ForNMDS2)=row.names(Osite_trait_counts)
Count_traits_ForNMDS2$status='Reference'



Count_traits_ForNMDS_all=rbind(Count_traits_ForNMDS1,Count_traits_ForNMDS2)
Count_traits_ForNMDS_all =Count_traits_ForNMDS_all[rowSums(Count_traits_ForNMDS_all[,-ncol(Count_traits_ForNMDS_all)])>0,]

Count_traits_ForNMDS_all=Count_traits_ForNMDS_all[match(env_dat3$sampleId, row.names(Count_traits_ForNMDS_all)),]
NMDS_Count_sor=vegan::vegdist(Count_traits_ForNMDS_all[,1:25],method='bray',binary=F)
set.seed(99)
NMDS_Count_MDS=vegan::metaMDS(NMDS_Count_sor)
set.seed(99)
NMDS_Countfit=envfit(NMDS_Count_MDS,scale(env_dat3[,!names(env_dat3)%in% c('COMID','sampleId')]))
NMDS_Countscores=as.data.frame(scores(NMDS_Count_MDS))

NMDS_Countscores$status=Count_traits_ForNMDS_all$status
env_vect<-as.data.frame(NMDS_Countfit$vectors$arrows)
env_vec=as.data.frame(scores(NMDS_Countfit, display = 'vectors'))
env_vect$env <- rownames(env_vect)
colnames(env_vect)[1:2] <- c("NMDS1", "NMDS2")

#Countscores2=Countscores[Countscores$NMDS1<35 & Countscores$NMDS2 < 2,]

arrow_multiplier <- 0.5 * max(abs(NMDS_Countscores[,c(1,2)])) / max(abs(env_vect[,1:2]))
env_vect[,1:2] <- env_vect[,1:2] * arrow_multiplier
env_vect$env=c('Catchment Elev',
               'Mean Precip',
               'Mean Temp','Basin Area',
               '% Crops',
               '% Urban')

env_vect <- env_vect %>%
  mutate(
    nudge_x = case_when(
      env == "Basin Area"    ~  0.5,
      env == "Catchment Elev"~ 0.3,
      env == "% Urban"       ~  0.2,
      env == "Mean Precip"   ~0,
      env== 'Mean Temp' ~ 0.4,
     env == '% Crops'~ 0.4,
     TRUE ~  0
    ),

    nudge_y = case_when(
      env=='% Urban' ~ 0.2,
      env=='% Mean Temp' ~ 0.5,
      env == "Catchment Elev"~ 0.2,
      env == "% Crops" ~ 0.028,
      env=='Mean Precip'~-0.1,
      TRUE             ~ 0
    )
  )
ggplot(NMDS_Countscores,aes(x=NMDS1,y=NMDS2,fill=status))+geom_point(pch=21,size=3)+
  scale_fill_manual(name='Status',
                    values=c('Reference' = 'blue',
                             'Probabilistic' = 'orange3'))+
  geom_segment(data=env_vect,
               aes(x=0,y=0,xend=NMDS1,yend=NMDS2),#                arrow=arrow(length=unit(0.25,'cm')),
               linewidth=1,
               alpha=0.9,inherit.aes = F, arrow=arrow(length=unit(.25,"centimeters")))+
   # geom_text_repel(
   #   data = env_vect,
   #   aes(
   #     x = NMDS1,
   #     y = NMDS2,
   #     label = env
   #   ),
   #   nudge_x = env_vect$nudge_x,
   #   nudge_y = env_vect$nudge_y,
   #   size = 6,
   #   inherit.aes = FALSE,
   #   box.padding = 0.1,
   #   point.padding = 0.1,
   #   min.segment.length = 0,
   #   max.overlaps = Inf
   # )+
  stat_ellipse(level=0.8,data = subset(NMDS_Countscores, status == "Reference"), color = "blue", size = 1,type='t') +
  stat_ellipse(level=0.8,data = subset(NMDS_Countscores, status == "Probabilistic"), color = "orange3", size = 1)+
  theme_classic()+
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.title = element_text(size = 25),
    axis.title.x = element_text(margin=margin(t=20)),
    axis.title.y  = element_text(margin=margin(r=20)),
    strip.text = element_text(size = 14),
    legend.text = element_text(size = 20),
    legend.title = element_text(size = 20),
    legend.key.size = unit(1, "lines"),
    axis.ticks=element_blank(),
    legend.position = 'inside',
    legend.position.inside = c(0.9,0.13)
  )#+stat_ellipse(level=0.8,type='t')
savp(10,8,'C://Users//andrew.caudillo.BUGLAB-I9//Box//NAMC//Research Projects//AIM//IncreaserDecreaser_OE//Updated_w_modelObj//New_for_SFS//Count_trait_NMDS260515.png')


summary(Osite_trait_pct[,1:37])

summ_mat=apply(Psite_trait_pct[,1:25], 2, function(x) {
  c(Min = min(x, na.rm = TRUE),
    Median = median(x, na.rm = TRUE),
    Mean = mean(x,na.rm=T),
    Max = max(x, na.rm = TRUE))
})
clipr::write_clip(summ_mat)
# Convert to data frame to keep the row labels intact
summary_df <- as.data.frame(summary_matrix)

Osite_trait_pct=read.csv('C://Users//andrew.caudillo.BUGLAB-I9//Box//NAMC//Research Projects//AIM//IncreaserDecreaser_OE//Updated_w_modelObj//New_for_SFS//SiteStatus_models//Pct_all_traits_ref.csv')
Osite_trait_pct=Osite_trait_pct[Osite_trait_pct$sampleId %in% c('HAWK-33','WE-131','WE-209')==F,]
Psite_trait_pct=read.csv('C://Users//andrew.caudillo.BUGLAB-I9//Box//NAMC//Research Projects//AIM//IncreaserDecreaser_OE//Updated_w_modelObj//New_for_SFS//SiteStatus_models//Pct_all_traits_probabilistic.csv')
Psite_trait_pct=Psite_trait_pct[Psite_trait_pct$sampleId %in% sites_w_no_bugs==F,]
summ_mat=apply(Psite_trait_pct[,1:25], 2, function(x) {
  c(Min = min(x, na.rm = TRUE),
    Median = median(x, na.rm = TRUE),
    Mean = mean(x,na.rm=T),
    Max = max(x, na.rm = TRUE))
})
clipr::write_clip(summ_mat)
