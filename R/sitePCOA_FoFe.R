#PcoA for sites using traits as variables
#this chunk is readind in the data,
#subseting, etc
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
traits<-read.csv('C://Users//andrew.caudillo.BUGLAB-I9//Box//NAMC//Research Projects//AIM//IncreaserDecreaser_OE//Updated_w_modelObj//OTU_21_traits.csv')
OTU21<-read.csv('C://Users//andrew.caudillo.BUGLAB-I9//Box//NAMC//Research Projects//AIM//IncreaserDecreaser_OE//OTU21.csv')
traits<-traits[traits$OTU %in% OTU21$taxonOTU,]
PA=PA[,names(PA) %in% c('status',traits$OTU)]

ref_a<-ref_a[ref_a$OTU %in% OTU21$taxonOTU,]
ref_PA<-ref_a
ref_PA[,3]<-ifelse(ref_a[,3]>0,1,0)
prob_a<-prob_a[prob_a$OTU %in% OTU21$taxonOTU,]
prob_PA<-prob_a
prob_PA[,3]<-ifelse(prob_a[,3]>0,1,0)
traits<-traits[,-1]
library(tidyverse)
#drop the OTU column from the trait table
#trait_table<-trait_table[,-2]
#make the trait table wide format (since the O/E table is wide)

wide_traits<-as.data.frame(traits |>
                             pivot_wider(names_from = attribute_name, values_from = attribute_value,
                                         values_fn = list) |>
                             unnest(cols = everything() ))
wide_traits=wide_traits[!duplicated(wide_traits$OTU),]
#add the OTU names to the dataset

#set the top5 traits identified from the Gower matrix
top5s<- c("Develop_slow_season",    "Adult_exit_absent",      "Survive_desiccation_no", "Occurance_drift_rare" ,
"Shape_not_streamline")
#subset out traits we won't be using,
#such as affinities, secondary feeding, etc.
wide_traits<- wide_traits %>%
  select(-matches(c("_af",'HBI','_sec_','_PH')))
#convert the T/F codes to 1 or 0 for numeric operations
taxa_traits_num<- lapply(wide_traits[,2:ncol(wide_traits)], function(x) {
  ifelse(x == "TRUE", 1, ifelse(x == "FALSE", 0, NA))
})
#now just combine the lapplied object into one data frame
taxa_traits_num<-as.data.frame(do.call(cbind,taxa_traits_num))
#set OTU so we can join
taxa_traits_num$OTU=wide_traits$OTU

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
sit_PA_gow<-cluster::daisy(site_PA_traits_weighted[,3:ncol(site_PA_traits_weighted)],metric='gower')

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
bug_dat=bug_dat[bug_dat$OTU %in% taxa_notraits_rare$taxon==F,]
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

#Gower_dist_abun<-cluster::daisy(bug_log_wide,metric='gower')
Gower_dist_abun=vegan::vegdist(bug_log_wide,'bray')
Gower_abun_mat<-as.matrix(Gower_dist_abun)
#assign the site names
row.names(Gower_abun_mat)<-row.names(bug_log_wide)
Gower_abun_MDS=vegan::metaMDS(Gower_abun_mat)

#Gower_dist_PA<-cluster::daisy(bug_PA_wide,metric='gower')#type=list(asymm=1:ncol(bug_PA_wide)))
Gower_dist_PA<-vegan::vegdist(bug_PA_wide,'bray',binary=T)
Gower_PA_mat<-as.matrix(Gower_dist_PA)
#assign the site names
row.names(Gower_PA_mat)<-row.names(bug_log_wide)
Gower_PA_MDS=vegan::metaMDS(Gower_PA_mat)

abun_scores=as.data.frame(vegan::scores(Gower_abun_MDS))
PA_scores=as.data.frame(vegan::scores(Gower_PA_MDS))

abun_scores$Status=statuses_wide
PA_scores$Status=statuses_wide

outlier=abun_scores[abun_scores$NMDS1 < -0.75,]
abun_scores<-abun_scores[abun_scores$NMDS1 > -0.75,]
PA_scores<-PA_scores[PA_scores$NMDS1 > -0.75,]
ggplot(abun_scores,aes(x=NMDS1,y=NMDS2,fill=Status))+geom_point(pch=21)+
  scale_fill_manual(name='Status',
                    values=c('Reference' = 'purple4',
                             'Probabilistic' = 'yellow'))+ggtitle('Sites in OTU space (PA)')+
  stat_ellipse(level=0.8,data = subset(abun_scores, Status == "Reference"), color = "purple4", size = 1,type='t') +
  stat_ellipse(level=0.8,data = subset(abun_scores, Status == "Probabilistic"), color = "yellow", size = 1) #+stat_ellipse(level=0.8,type='t')+scale_color_manual(values=c('Reference'='red','Probabilistic' = 'blue'))
savp(10,8, 'C://Users//andrew.caudillo.BUGLAB-I9//Box//NAMC//Research Projects//AIM//IncreaserDecreaser_OE//sites_OTUspace_abun_ellipse.png')

ggplot(PA_scores,aes(x=NMDS1,y=NMDS2,fill=Status))+geom_point(pch=21)+
  scale_fill_manual(name='Status',
                    values=c('Reference' = 'purple4',
                             'Probabilistic' = 'yellow'))+ggtitle('Sites in OTU space (PA)')+
stat_ellipse(level=0.8,data = subset(PA_scores, Status == "Reference"), color = "purple4", size = 1,type='t') +
  stat_ellipse(level=0.8,data = subset(PA_scores, Status == "Probabilistic"), color = "yellow", size = 1)#+stat_ellipse(level=0.8,type='t')
savp(10,8, 'C://Users//andrew.caudillo.BUGLAB-I9//Box//NAMC//Research Projects//AIM//IncreaserDecreaser_OE//sites_OTUspace_PA_ellipse.png')

