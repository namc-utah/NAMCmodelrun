#PcoA for sites using traits as variables
#this chunk is readind in the data,
#subseting, etc
top5s=c("Emerge_synch_abbrev_Poorly",  "Habit_prim_Swimmer",          "Crawl_rate_high",
        "Voltinism_abbrev_Univoltine", "Feed_prim_abbrev_PR")
if(0){
abun<-read.csv("C://Users//andrew.caudillo.BUGLAB-I9//Box//NAMC//Research Projects//AIM//IncreaserDecreaser_OE//abundances.csv")
abun$sample_id<-as.character(abun$sample_id)
abun<-abun[,-1]
}
refs=Os
Probs=ProbOs
refs$status='Reference'
Probs$status='Probabilistic'
PA=rbind(refs,Probs)#PA<-ifelse(abun>0,1,0)
failed_sites<-read.csv('C://Users//andrew.caudillo.BUGLAB-I9//Box//NAMC//Research Projects//AIM//IncreaserDecreaser_OE//failed_sites.csv')
PA=PA[row.names(PA) %in% failed_sites$sampleId==F,]
#ref_a<-abun[abun$Status=='Reference',]
#ref_a<-ref_a[ref_a$sample_id %in% c(116830,132243,132807, 185022)==F,]
#prob_a<-abun[abun$Status=='Probabilistic',]
#prob_a<-prob_a[prob_a$sample_id %in% failed_sites$sampleId==F,]
# traits<-read.csv('C://Users//andrew.caudillo.BUGLAB-I9//Box//NAMC//Research Projects//AIM//IncreaserDecreaser_OE//Updated_w_modelObj//OTU_21_traits.csv')
# OTU21<-read.csv('C://Users//andrew.caudillo.BUGLAB-I9//Box//NAMC//Research Projects//AIM//IncreaserDecreaser_OE//OTU21.csv')
traits=maxxed_OTUs
#traits<-traits[traits$OTU %in% OTU21$taxonOTU,]
traits=traits[!duplicated(traits$OTU),]
PA=PA[,names(PA) %in% c('status',traits$OTU)]

ref_a<-PA[PA$status=='Reference',]

library(tidyverse)
#drop the OTU column from the trait table
#trait_table<-trait_table[,-2]
#make the trait table wide format (since the O/E table is wide)
#
# wide_traits<-as.data.frame(traits |>
#                              pivot_wider(names_from = attribute_name, values_from = attribute_value,
#                                          values_fn = list) |>
#                              unnest(cols = everything() ))
# wide_traits=wide_traits[!duplicated(wide_traits$OTU),]
#add the OTU names to the dataset

#set the top5 traits identified from the Gower matrix
#
#subset out traits we won't be using,
#such as affinities, secondary feeding, etc.
# wide_traits<- wide_traits %>%
#   select(-matches(c("_af",'HBI','_sec_','_PH')))
# #convert the T/F codes to 1 or 0 for numeric operations
# taxa_traits_num<- lapply(wide_traits[,2:ncol(wide_traits)], function(x) {
#   ifelse(x == "TRUE", 1, ifelse(x == "FALSE", 0, NA))
# })

#now just combine the lapplied object into one data frame
taxa_traits_num=traits
#set OTU so we can join
#taxa_traits_num$OTU=wide_traits$OTU

#optional#
#get just the top 5 traits
taxa_traits_num<-taxa_traits_num[,names(taxa_traits_num) %in% c('OTU', top5s)]
#look at just the taxa who appear at > 10 sites and have traits.
#taxa_traits_num<-taxa_traits_num[taxa_traits_num$OTU %in% taxa_notraits_rare$taxon==F,]
#subset sampleId for later, if needed
#ref_a_samps<-ref_a$sample_id
#ref_a<-ref_a[,-2]
#join the traits to the respsective community data
PA_long=PA %>%
  tibble::rownames_to_column(var='sample_id') %>%
  pivot_longer(cols=-c(sample_id,status),
                    names_to = 'OTU',values_to = 'PA')
joined<-plyr::join(PA_long,taxa_traits_num,type='left')

# ref_P_joined<-plyr::join(ref_PA,taxa_traits_num,type='left')
# p_joined<-plyr::join(prob_a,taxa_traits_num,type='left')
# p_P_joined<-plyr::join(prob_PA,taxa_traits_num,type='left')

#set the trait columns
#for dplyr to calculate the weighted traits
trait_cols=5:ncol(joined)
trait_cols=names(joined)[trait_cols]

library(dplyr)

#row.names(taxa_traits_num)<-taxa_traits_num$OTU;taxa_traits_num=taxa_traits_num[,names(taxa_traits_num) %in% 'OTU'==F]

#weighting the traits by abundnace or PA
#multiply the trait relative abundnace for each taxon
#then sum it all for a community weighted mean
site_traits_weighted <- joined %>%
  # mutate(across(all_of(trait_cols), ~ .x * split_count)) %>%
  group_by(sample_id) %>%
  mutate(
    rel_abund = split_count / sum(split_count, na.rm = TRUE)) %>%
  dplyr::summarise(
    # abundance-weighted trait per site
    across(all_of(trait_cols), ~ sum(.x * rel_abund, na.rm = TRUE), .names = "abund_{.col}"),
    .groups='drop')

  # dplyr::summarise(across(all_of(trait_cols),
  #                         ~ sum(.x, na.rm = TRUE) / sum(split_count, na.rm = TRUE)),
  #                  .groups = "drop")
# joined %>%
#   mutate(across(where(is.numeric), ~.x*split_count))


# trait_cols=5:ncol(ref_P_joined)
# trait_cols=names(ref_P_joined)[trait_cols]
 site_PA_traits_weighted<-joined %>%
  #mutate(across(all_of(trait_cols), ~ .x * (split_count))) %>%  # presence = 1/0
  group_by(sample_id,status) %>%
   dplyr::summarise(
     # PA-weighted trait per site
     across(all_of(trait_cols), ~ sum(.x * PA, na.rm = TRUE) / sum(PA, na.rm = TRUE), .names = "PA_{.col}"),
     .groups = "drop"
   )
 #   dplyr::summarise(across(all_of(trait_cols),
 #                           ~ sum(.x, na.rm = TRUE) / sum(split_count, na.rm = TRUE)),
 #                    .groups = "drop")
 # # joined %>%
  # ref_P_joined%>%
  # mutate(across(where(is.numeric), ~.x*split_count))



# site_trait_agg <- site_traits_weighted %>%
#   group_by(sample_id) %>%
#   dplyr::summarise(across(where(is.numeric), sum, na.rm = TRUE))
#
# site_PA_trait_agg <- site_PA_traits_weighted %>%
#   group_by(sample_id) %>%
#   dplyr::summarise(across(where(is.numeric), sum, na.rm = TRUE))
# nrow(site_trait_agg)
 trait_cols=5:ncol(p_joined)
 trait_cols=names(p_joined)[trait_cols]
psite_traits_weighted <- p_joined %>%
  #mutate(across(all_of(trait_cols), ~ .x * (split_count))) %>%  # presence = 1/0
  group_by(sample_id) %>%
  mutate(
    rel_abund = split_count / sum(split_count, na.rm = TRUE)) %>%
  dplyr::summarise(
    # abundance-weighted trait per site
    across(all_of(trait_cols), ~ sum(.x * rel_abund, na.rm = TRUE), .names = "abund_{.col}"),
    .groups='drop')
  # dplyr::summarise(across(all_of(trait_cols), ~ sum(.x, na.rm = TRUE) / sum(split_count > 0, na.rm = TRUE)),
  #                  .groups = "drop")
# ref_P_joined%>%
# p_joined %>%
#   mutate(across(where(is.numeric), ~.x*split_count))
trait_cols=5:ncol(p_P_joined)
trait_cols=names(p_P_joined)[trait_cols]
psite_PA_traits_weighted<- p_P_joined %>%
  #mutate(across(all_of(trait_cols), ~ .x * (split_count))) %>%  # presence = 1/0
  group_by(sample_id) %>%
  dplyr::summarise(
    # PA-weighted trait per site
    across(all_of(trait_cols), ~ sum(.x * split_count, na.rm = TRUE) / sum(split_count, na.rm = TRUE), .names = "PA_{.col}"),
    .groups = "drop"
  )
#   dplyr::summarise(across(all_of(trait_cols), ~ sum(.x, na.rm = TRUE) / sum(split_count > 0, na.rm = TRUE)),
#                    .groups = "drop")
# # ref_P_joined%>%
# p_P_joined %>%
#   mutate(across(where(is.numeric), ~.x*split_count))
# psite_trait_agg <- psite_traits_weighted %>%
#   group_by(sample_id) %>%
#   dplyr::summarise(across(where(is.numeric), sum, na.rm = TRUE))
# psite_PA_trait_agg <- psite_PA_traits_weighted %>%
#   group_by(sample_id) %>%
#   dplyr::summarise(across(where(is.numeric), sum, na.rm = TRUE))

#change the tibbles to data frames
site_traits_weighted<-as.data.frame(site_traits_weighted)
site_PA_traits_weighted<-as.data.frame(site_PA_traits_weighted)
psite_traits_weighted<-as.data.frame(psite_traits_weighted)
psite_PA_traits_weighted<-as.data.frame(psite_PA_traits_weighted)
# psite_trait_agg<-as.data.frame(psite_trait_agg)
# site_trait_agg<-as.data.frame(site_trait_agg)
# p_P_site_trait_agg<-as.data.frame(psite_PA_trait_agg)
# site_PA_trait_agg<-as.data.frame(site_PA_trait_agg)

#assign a status to each set of sites
site_traits_weighted$status='Reference'
site_PA_traits_weighted$status<-'Reference'
psite_traits_weighted$status='Probabilistic'
psite_PA_traits_weighted$status='Probabilistic'

#bind the reference and probabilistic sites by rows
All_site_traits_weighted<-rbind(site_traits_weighted,psite_traits_weighted)
All_PA_site_traits_weighted<-rbind(site_PA_traits_weighted,psite_PA_traits_weighted)
row.names(All_site_traits_weighted)<-All_site_traits_weighted$sample_id
row.names(All_PA_site_traits_weighted)<-All_PA_site_traits_weighted$sample_id
#All_site_trait_small<-All_site_traits_weighted[,-2]

#remove the columns we do not need anymore (sampleId,status)
#because these are going into the Gower matrix
All_site_trait_small<-All_site_traits_weighted[,-c(1,ncol(All_site_traits_weighted))]

All_PA_site_trait_small<-All_PA_site_traits_weighted[,-c(1,ncol(All_PA_site_traits_weighted))]
#run the gower dissimilarity matrix on the trait datasets
sit_gow<-cluster::daisy(All_site_trait_small,metric='gower')
site_PA_traits_weighted=na.omit(site_PA_traits_weighted)
sit_PA_gow<-cluster::daisy(site_PA_traits_weighted[,3:ncol(site_PA_traits_weighted)],metric='gower')
sit_PA_NMDS=vegan::vegdist(sitepa)
#row.names(sit_PA_gow)=site_PA_traits_weighted$sample_id

#compute PCOA scores
bc<-ape::pcoa(sit_gow)
BC_scores<-as.data.frame(bc$vectors[,1:2])
Pbc<-ape::pcoa(sit_PA_gow)
pBC_scores<-as.data.frame(Pbc$vectors[,1:2])
pBC_scores$status=site_PA_traits_weighted$status
BC_scores$status=All_site_traits_weighted$status
pBC_scores$status=All_PA_site_traits_weighted$status
row.names(BC_scores)=All_site_traits_weighted$sample_id
row.names(pBC_scores)=All_site_traits_weighted$sample_id
#plot the results, color them by status, and draw 80% confidence ellipses
ggplot(BC_scores,aes(x=Axis.1,y=Axis.2,fill=status))+geom_point(pch=21)+
         scale_fill_manual(name='Status',
                            values=c('Reference' = 'purple4',
                                     'Probabilistic' = 'yellow'))+
  ggtitle('Site types in trait space\n(top 5 traits weighted by abundance)')+
  stat_ellipse(level=0.8,data = subset(BC_scores, status == "Reference"), color = "purple4", size = 1,type='t') +
  stat_ellipse(level=0.8,data = subset(BC_scores, status == "Probabilistic"), color = "yellow", size = 1,type='t')


savp(10,8, 'C://Users//andrew.caudillo.BUGLAB-I9//Box//NAMC//Research Projects//AIM//IncreaserDecreaser_OE//sites_abund_top5_traitspace.png')
ggplot(pBC_scores,aes(x=Axis.1,y=Axis.2,fill=status))+geom_point(pch=21)+
  scale_fill_manual(name='Status',
                    values=c('Reference' = 'purple4',
                             'Probabilistic' = 'yellow'))+
  ggtitle('Site types in trait space\n(top 5 traits weighted by presence)')+
  stat_ellipse(level=0.8,data = subset(pBC_scores, status == "Reference"), color = "purple4", size = 1,type='t') +
  stat_ellipse(level=0.8,data = subset(pBC_scores, status == "Probabilistic"), color = "yellow", size = 1,type='t')

savp(10,8, 'C://Users//andrew.caudillo.BUGLAB-I9//Box//NAMC//Research Projects//AIM//IncreaserDecreaser_OE//Updated_w_modelObj//sites_PA_top5_traitspace.png')



ref_counts=aggregate(split_count~sample_id,data=ref_a[ref_a$OTU %in% taxa_notraits_rare$taxon==F,],FUN=sum)
p_counts=aggregate(split_count~sample_id,data=prob_a[prob_a$OTU %in% taxa_notraits_rare$taxon==F,],FUN=sum)
ref_rich=aggregate(split_count~sample_id,data=ref_PA[ref_PA$OTU %in% taxa_notraits_rare$taxon==F,],FUN=sum)
p_rich=aggregate(split_count~sample_id,data=prob_PA[prob_PA$OTU %in% taxa_notraits_rare$taxon==F,],FUN=sum)
par(mfrow=c(2,1))
hist(ref_counts$split_count,xlab='',main='Reference')
hist(p_counts$split_count,main='Probabilistic',xlab='split count')



graphics.off()
boxplot(ref_counts$split_count,at=1,xlim=c(0,3),col='yellow3',ylab='Split Count')
boxplot(p_counts$split_count,at=2,add=T,col='purple4')
axis(side=1,at=c(1,2),labels=c('Reference','Probabilistic'))
savp(10,8, 'C://Users//andrew.caudillo.BUGLAB-I9//Box//NAMC//Research Projects//AIM//IncreaserDecreaser_OE//sites_abund_boxes.png')

boxplot(ref_rich$split_count,at=1,xlim=c(0,3),col='yellow3',ylab='Richness',ylim=c(0,max(p_rich$split_count)))
boxplot(p_rich$split_count,at=2,add=T,col='purple4')
axis(side=1,at=c(1,2),labels=c('Reference','Probabilistic'))
savp(10,8, 'C://Users//andrew.caudillo.BUGLAB-I9//Box//NAMC//Research Projects//AIM//IncreaserDecreaser_OE//sites_rich_boxes.png')



####NMDS plot####
bug_dat=rbind(ref_a,prob_a)
bug_agg=aggregate(split_count ~ sample_id + Status + OTU, data = bug_dat, sum)
bug_wide=reshape(bug_agg,
        idvar=c('sample_id','Status'),
        timevar='OTU',
        direction='wide',
        )
bug_wide[is.na(bug_wide)] <- 0
names(bug_wide) <- sub("^split_count\\.", "", names(bug_wide))


statuses_wide=bug_wide$Status;row.names(bug_wide)<-bug_wide$sample_id;bug_wide=bug_wide[,names(bug_wide) %in% c('sample_id','Status')==F]
bug_PA_wide=as.data.frame(ifelse(bug_wide > 0, 1, 0))
bug_log_wide<-log(bug_wide+1)

PA3=PA[,names(PA)[1:(ncol(PA)-1)] %in% row.names(traits_only2)]

#Gower_dist_abun<-cluster::daisy(bug_log_wide,metric='gower')
Gower_dist_abun=vegan::vegdist(bug_log_wide,'bray')

#assign the site names
row.names(Gower_abun_mat)<-row.names(bug_log_wide)
Gower_abun_MDS=vegan::metaMDS(Gower_abun_mat)
PA3=PA3[rowSums(PA3[,-ncol(PA3)])>0,] #only present at failed sites, so artifact.
#Gower_dist_PA<-cluster::daisy(bug_PA_wide,metric='gower')#type=list(asymm=1:ncol(bug_PA_wide)))
Vegdist_PA<-vegan::vegdist(PA3[,1:(ncol(PA3)-1)],'bray',binary=T)
Veg_PA_mat<-as.matrix(Vegdist_PA)
#assign the site names
row.names(Veg_PA_mat)<-row.names(PA3)
PA_MDS=vegan::metaMDS(Vegdist_PA,k=3)

#abun_scores=as.data.frame(vegan::scores(Gower_abun_MDS))
PA_scores=as.data.frame(vegan::scores(PA_MDS))


PA_scores$Status=PA3$status

NMDS_cal_preds=read.csv("C://Users//andrew.caudillo.BUGLAB-I9//Box//NAMC//Research Projects//AIM//IncreaserDecreaser_OE//Updated_w_modelObj//Cal_predictors.csv")
NMDS_prob_preds=read.csv("C://Users//andrew.caudillo.BUGLAB-I9//Box//NAMC//Research Projects//AIM//IncreaserDecreaser_OE//Updated_w_modelObj//Prob_predictors_prednew.csv")
#prob_comids=read.csv('C://Users//andrew.caudillo.BUGLAB-I9//Box//NAMC//Research Projects//AIM//IncreaserDecreaser_OE//prob_comids.csv')
#NMDS_prob_preds$sampleId=row.names(NMDS_prob_preds)
NMDS_prob_preds=NMDS_prob_preds[NMDS_prob_preds$sampleId %in% failed_sites$sampleId==F,]
NMDS_prob_preds=plyr::join(NMDS_prob_preds,prob_comids,by='sampleId')

#row.names(NMDS_prob_preds)<-NMDS_prob_preds$sampleId

#row.names(NMDS_prob_preds)==row.names(PA)[PA$status=='Probabilistic']
cal_pred_dat=read.csv('C://Users//andrew.caudillo.BUGLAB-I9//Box//NAMC//Research Projects//AIM//IncreaserDecreaser_OE//Updated_w_modelObj//WW_cal_preds.csv')
NMDS_cal_preds<-plyr::join(NMDS_cal_preds,cal_pred_dat,by='COMID','left')
names(NMDS_cal_preds)[2]='sampleId'
env_dat=rbind(NMDS_cal_preds[,c('sampleId','ElevCat','Precip8110Ws','Tmean8110Ws','WsAreaSqKm','COMID')],NMDS_prob_preds[,c('sampleId','ElevCat','Precip8110Ws','Tmean8110Ws','WsAreaSqKm','COMID')])
AgUrb_dat=read.csv('C://Users//andrew.caudillo.BUGLAB-I9//Box//NAMC//Research Projects//AIM//IncreaserDecreaser_OE//all_site_UrbAg.csv')
env_dat=plyr::join(env_dat,AgUrb_dat[,c('COMID','PCtCrop','AvgUrb')],by='COMID','left')
env_dat=env_dat[env_dat$sampleId %in% site_PA_traits_weighted$sample_id,]

env_dat3=env_dat[env_dat$sampleId %in% row.names(PA3),]
row.names(PA3)==env_dat3$sampleId
nmds_fit=envfit(PA_MDS,env_dat3[,!names(env_dat3)%in% c('COMID','sampleId')])


# nmds_fitA=envfit(PA_scores[,c(1,2)],scale(env_dat3[,names(env_dat3) %in% c('ElevCat','Precip8110Ws',
#                                                                   'Tmean8110Ws','WsAreaSqKm',
#                                                                   'PCtCrop',
#                                                                   'AvgUrb')]),na.rm=T)
env_vect<-as.data.frame(scores(nmds_fit,display = 'vectors'))

env_vect$env <- rownames(env_vect)
colnames(env_vect)[1:2] <- c("NMDS1", "NMDS2")
mult=vegan::ordiArrowMul(nmds_fit) *0.2
PA_scores2=PA_scores[PA_scores$NMDS1<34 & PA_scores$NMDS2 < 2,]
env_vect[,1:2] <- env_vect[,1:2] /
  sqrt(rowSums(env_vect[,1:2]^2))
env_vect[,1:2] <- env_vect[,1:2] * .5
env_vect[,1:2] <- env_vect[,1:2] * mult
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
      TRUE             ~ 0
    )
  )


# ggplot(abun_scores,aes(x=NMDS1,y=NMDS3,fill=Status))+geom_point(pch=21)+
#   scale_fill_manual(name='Status',
#                     values=c('Reference' = 'purple4',
#                              'Probabilistic' = 'yellow'))+ggtitle('Sites in OTU space (PA)')+
#   stat_ellipse(level=0.8,data = subset(abun_scores, Status == "Reference"), color = "purple4", size = 1,type='t') +
#   stat_ellipse(level=0.8,data = subset(abun_scores, Status == "Probabilistic"), color = "yellow", size = 1) #+stat_ellipse(level=0.8,type='t')+scale_color_manual(values=c('Reference'='red','Probabilistic' = 'blue'))
# savp(10,8, 'C://Users//andrew.caudillo.BUGLAB-I9//Box//NAMC//Research Projects//AIM//IncreaserDecreaser_OE//sites_OTUspace_abun_ellipse.png')
#windows(10,8)
ggplot(PA_scores2,aes(x=NMDS1,y=NMDS2,fill=Status))+geom_point(pch=21,size=3)+
  scale_fill_manual(name='Status',
                    values=c('Reference' = 'blue',
                             'Probabilistic' = 'orange3'))+
  geom_segment(data=env_vect,
               aes(x=0,y=0,xend=NMDS1,yend=NMDS2),#                arrow=arrow(length=unit(0.25,'cm')),
               linewidth=1,
               alpha=0.5,inherit.aes = F, arrow=arrow(length=unit(.25,"centimeters")))+
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
stat_ellipse(level=0.8,data = subset(PA_scores, Status == "Reference"), color = "blue", size = 1,type='t') +
  stat_ellipse(level=0.8,data = subset(PA_scores, Status == "Probabilistic"), color = "orange3", size = 1)+
  theme_classic()+
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.title = element_text(size = 25),
    axis.title.x = element_text(margin=margin(t=20)),
    axis.title.y  = element_text(margin=margin(r=20)),
    strip.text = element_text(size = 14),
    legend.text = element_text(size = 15),
    legend.title = element_text(size = 15),
    legend.key.size = unit(1, "lines"),
    axis.ticks=element_blank(),
    legend.position = 'inside',
    legend.position.inside = c(max(PA_scores2$NMDS1)+0.28,0.2)
    )#+stat_ellipse(level=0.8,type='t')

savp(10,8, 'C://Users//andrew.caudillo.BUGLAB-I9//Box//NAMC//Research Projects//AIM//IncreaserDecreaser_OE//Updated_w_modelObj//sites_OTUspace_NMDS_axes_1_2_wAgUrb_260513.png')


ordisurf(PA_scores[,c(1,2)], env_dat$PCtCrop,k=4,
         pch=21,bg=ifelse(PA$status=='Reference','purple4','yellow'),
         col='black',
         main='PctCrop',
         lwd=3)
legend('topleft',
       pch=rep(21,2),
       pt.bg=c('purple4','yellow'),
       bty='n',
       leg=c('Reference',
             'Probabilistic'),
       cex=0.8)

savp(10,8,'C://Users//andrew.caudillo.BUGLAB-I9//Box//NAMC//Research Projects//AIM//IncreaserDecreaser_OE//Updated_w_modelObj//PA_ordi_Crop.png')
Pgower_scores=data.frame(Pbc$vectors)
#gower_scores=data.frame(Pco1=gower_cmd[,1],Pco2=gower_cmd[,2])
Pgower_scores$taxon=site_PA_traits_weighted$sample_id

Pbc_cmd=cmdscale(sit_PA_gow,k=3)
env_fut=envfit(ord=Pgower_scores[,c(1,2)],env=site_PA_traits_weighted[,3:ncol(site_PA_traits_weighted)],eig=T)
ptrait_vect<-as.data.frame(env_fut$vectors$arrows)
ptrait_vec=as.data.frame(scores(env_fut,display='vectors'))
ptrait_vect$trait <- rownames(ptrait_vect)
colnames(ptrait_vect)[1:2] <- c("PCo1", "PCo2")
ptrait_vect=ptrait_vect[ptrait_vect$trait %in% paste0('PA_',top5s),]
#ptrait_vect[,1:2]<-ptrait_vect[,1:2]*-1
ptrait_vect$trait=c('Poor Synch. Emerge',
                   'Predator',
                   'Swimmer',
                   'Univoltine',
                   'Crawler')
pPCOA_Scores<-as.data.frame(Pbc$vectors)
#assign taxa names
ptop_trait_arrows=ptrait_vect$trait
pPCOA_Scores$sample_id<-site_PA_traits_weighted$sample_id

p_sites_and_coords=plyr::join(site_PA_traits_weighted,pPCOA_Scores)



point_max <- max(abs(Pgower_scores[,1:2]))
arrow_max <- max(abs(ptrait_vect[,1:2]))
arrow_multiplier <- (point_max / arrow_max) * 0.4  # 0.8 keeps them inside the points

# Apply scaling
ptrait_vect$PCo1 <- ptrait_vect$PCo1 * arrow_multiplier
ptrait_vect$PCo2 <- ptrait_vect$PCo2 * arrow_multiplier
names(p_sites_and_coords)[2]<-'Status'
# --- 2. Plotting with corrected scales ---
ggplot(p_sites_and_coords, aes(x=Axis.1, y=Axis.2, fill=Status)) +
  geom_point(pch=21,size=3) +
  stat_ellipse(level=0.8,data = subset(p_sites_and_coords, Status == "Reference"), color = "blue", size = 1,type='t') +
  stat_ellipse(level=0.8,data = subset(p_sites_and_coords, Status == "Probabilistic"), color = "orange3", size = 1,type='t')+# Added alpha for fill
  # FIX: Change fill to color to match aes(color=status)
  scale_fill_manual(name='Status',
                     values=c('Reference' = 'blue',
                              'Probabilistic' = 'orange3')) +
  geom_segment(data=ptrait_vect,
               aes(x=0, y=0, xend=PCo1, yend=PCo2),
               arrow=arrow(length=unit(0.25, 'cm')),
               linewidth=1,
               inherit.aes = FALSE,
               color="black") +

  geom_text_repel(
    data = ptrait_vect,
    aes(
      x = PCo1,
      y = PCo2,
      label = trait
    ),
    # Use your logic to "nudge" the starting position
    nudge_y = ifelse(ptrait_vect$trait == 'Univoltine',0.02,0),
     nudge_x = ifelse(ptrait_vect$trait == 'Crawler', 0.02,0) ,
    #                  ifelse(ptrait_vect$trait == 'Swimmer', 0.5, ifelse(ptrait_vect$vect=='Poor Synch. Emerge',0.8,0))),
    segment.color = "transparent", # Hides the segments entirely
    size = 6,
    inherit.aes = FALSE,
    box.padding = 0.2,    # Extra "breathing room" around labels
    point.padding = 0.1,  # Distance from the data point
    min.segment.length = 0, # Always draw lines to labels if they move far
    max.overlaps = Inf    # Forces all labels to show
  )+# Set a fixed color so it doesn't look for 'status'
  theme_classic()+
  labs(x='PCoA Axis 1',y='PCoA Axis 2')+
  theme(#axis.text.x = element_text(size = 14),
        axis.ticks = element_blank(),
        axis.text=element_blank(),
      #  axis.text.y = element_text(size = 14),
        axis.title = element_text(size = 25),
        legend.position = 'inside',
        legend.position.inside = c(0.15,0.9),
      axis.title.x = element_text(margin=margin(t=20)),
      axis.title.y  = element_text(margin=margin(r=20)),
      legend.text = element_text(size = 15),
      legend.title = element_text(size = 15))#sures 1 unit on X = 1 unit on Y
savp(10,8, 'C://Users//andrew.caudillo.BUGLAB-I9//Box//NAMC//Research Projects//AIM//IncreaserDecreaser_OE//Updated_w_modelObj//sites_OTUspace_PA_ellipse_traits_axes1_2_260512.png')





#curious to see what taxa appear in each quadrant of an ordination?
top5traits=traits[traits$attribute_name %in% names(top5),]
okay=top5traits %>%
  pivot_wider(id_cols = OTU,
              names_from = attribute_name,
              values_from = attribute_value,
              )
