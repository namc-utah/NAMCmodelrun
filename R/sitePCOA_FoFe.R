#PcoA for sites
abun<-read.csv("C://Users//andrew.caudillo.BUGLAB-I9//Box//NAMC//Research Projects//AIM//IncreaserDecreaser_OE//abundances.csv")
abun$sample_id<-as.character(abun$sample_id)
abun<-abun[,-1]
PA<-ifelse(abun>0,1,0)
failed_sites<-read.csv('C://Users//andrew.caudillo.BUGLAB-I9//Box//NAMC//Research Projects//AIM//IncreaserDecreaser_OE//failed_sites.csv')
ref_a<-abun[abun$Status=='Reference',]
ref_a<-ref_a[ref_a$sample_id %in% c(116830,132243,132807, 185022)==F,]
prob_a<-abun[abun$Status=='Probabilistic',]
prob_a<-prob_a[prob_a$sample_id %in% failed_sites$sampleId==F,]
traits<-read.csv('C://Users//andrew.caudillo.BUGLAB-I9//Box//NAMC//Research Projects//AIM//IncreaserDecreaser_OE//All_bench_taxa_atrributes2.csv')
OTU21<-read.csv('C://Users//andrew.caudillo.BUGLAB-I9//Box//NAMC//Research Projects//AIM//IncreaserDecreaser_OE//OTU21.csv')
traits<-traits[traits$OTU %in% OTU21$taxonOTU,]
ref_a<-ref_a[ref_a$OTU %in% OTU21$taxonOTU,]
ref_PA<-ref_a
ref_PA[,3]<-ifelse(ref_a[,3]>0,1,0)
prob_a<-prob_a[prob_a$OTU %in% OTU21$taxonOTU,]
prob_PA<-prob_a
prob_PA[,3]<-ifelse(prob_a[,3]>0,1,0)
traits<-traits[,-1]
#drop the OTU column from the trait table
#trait_table<-trait_table[,-2]
#make the trait table wide format (since the O/E table is wide)
if(0){
wide_traits<-as.data.frame(traits |>
                             pivot_wider(names_from = attribute_name, values_from = attribute_value,
                                         values_fn = list) |>
                             unnest(cols = everything() ))
#add the OTU names to the dataset

#rename "OTU" to be "Taxon" so we can join the tables
top5<- c("Swim_strong",           "Develop_slow_season",   "Habit_prim_Clinger",    "Rheophily_abbrev_eros",
         "Lifespan_short"    )
wide_traits<- wide_traits %>%
  select(-matches(c("_af",'HBI','_sec_','_PH')))
#forcing everything to be binary
wide_traits[wide_traits=='FALSE']<-as.numeric(0)
wide_traits[wide_traits=='TRUE']<-as.numeric(1)
OTUs<-wide_traits$OTU
wide_traits<-wide_traits[,-1]
taxa_traits_num <- wide_traits %>%
  mutate(across(where(is.character), ~ as.numeric(.x))) %>%   # converts "0","1","NA" to numeric 0,1,NA
  mutate(across(where(is.logical), ~ as.numeric(.x)))
taxa_traits_num$OTU=OTUs
taxa_traits_num<-taxa_traits_num[,names(taxa_traits_num) %in% c('OTU', top5)]
taxa_traits_num<-taxa_traits_num[taxa_traits_num$OTU %in% taxa_notraits_rare$taxon==F,]
ref_a_samps<-ref_a$sample_id
#ref_a<-ref_a[,-2]

joined<-plyr::join(ref_a,taxa_traits_num,type='left')
ref_P_joined<-plyr::join(ref_PA,taxa_traits_num,type='left')
p_joined<-plyr::join(prob_a,taxa_traits_num,type='left')
p_P_joined<-plyr::join(prob_PA,taxa_traits_num,type='left')

library(dplyr)
site_traits_weighted <- joined %>%
  mutate(across(where(is.numeric), ~.x*split_count))
site_PA_traits_weighted<-ref_P_joined%>%
  mutate(across(where(is.numeric), ~.x*split_count))
site_trait_agg <- site_traits_weighted %>%
  group_by(sample_id) %>%
  summarise(across(where(is.numeric), sum, na.rm = TRUE))

site_PA_trait_agg <- site_PA_traits_weighted %>%
  group_by(sample_id) %>%
  summarise(across(where(is.numeric), sum, na.rm = TRUE))
nrow(site_trait_agg)

psite_traits_weighted <- p_joined %>%
  mutate(across(where(is.numeric), ~.x*split_count))
psite_PA_traits_weighted<-p_P_joined %>%
  mutate(across(where(is.numeric), ~.x*split_count))
psite_trait_agg <- psite_traits_weighted %>%
  group_by(sample_id) %>%
  summarise(across(where(is.numeric), sum, na.rm = TRUE))
psite_PA_trait_agg <- psite_PA_traits_weighted %>%
  group_by(sample_id) %>%
  summarise(across(where(is.numeric), sum, na.rm = TRUE))

psite_trait_agg<-as.data.frame(psite_trait_agg)
site_trait_agg<-as.data.frame(site_trait_agg)
p_P_site_trait_agg<-as.data.frame(psite_PA_trait_agg)
site_PA_trait_agg<-as.data.frame(site_PA_trait_agg)
site_trait_agg$status='Reference'
site_PA_trait_agg$status<-'Reference'
psite_trait_agg$status='Probabilistic'
p_P_site_trait_agg$status='Probabilistic'

All_site_trait_agg<-rbind(site_trait_agg,psite_trait_agg)
All_PA_site_trait_agg<-rbind(site_PA_trait_agg,p_P_site_trait_agg)
row.names(All_site_trait_agg)<-All_site_trait_agg$sample_id
row.names(All_PA_site_trait_agg)<-All_PA_site_trait_agg$sample_id
#All_site_trait_small<-All_site_trait_agg[,-2]

statuses<-All_site_trait_agg$status

All_site_trait_small<-All_site_trait_agg[,-c(1,2,ncol(All_site_trait_agg))]
  #as.data.frame(lapply(All_site_trait_small, factor))
#All_site_trait_U<-unique(All_site_trait_small)
All_PA_site_trait_small<-All_PA_site_trait_agg[,-c(1,2,ncol(All_PA_site_trait_agg))]
sit_gow<-cluster::daisy(All_site_trait_small,metric='gower')
sit_PA_gow<-cluster::daisy(All_PA_site_trait_small,metric='gower')


bc<-ape::pcoa(sit_gow)
BC_scores<-as.data.frame(bc$vectors[,1:2])
Pbc<-ape::pcoa(sit_PA_gow)
pBC_scores<-as.data.frame(Pbc$vectors[,1:2])
BC_scores$status=All_site_trait_agg$status
pBC_scores$status=All_PA_site_trait_agg$status
row.names(BC_scores)=All_site_trait_agg$sample_id
row.names(pBC_scores)=All_site_trait_agg$sample_id
ggplot(BC_scores,aes(x=Axis.1,y=Axis.2,fill=status))+geom_point(pch=21)+
         scale_fill_manual(name='Status',
                            values=c('Reference' = 'purple4',
                                     'Probabilistic' = 'yellow'))+
  ggtitle('Site types in top 5 trait space\n(traits weighted by abundance)')
  #stat_ellipse(level=0.95,size=1)


ggplot(pBC_scores,aes(x=Axis.1,y=Axis.2,fill=status))+geom_point(pch=21)+
  scale_fill_manual(name='Status',
                    values=c('Reference' = 'purple4',
                             'Probabilistic' = 'yellow'))+
  ggtitle('Site types in top 5 trait space\n(traits weighted by presence)')
  #stat_ellipse(level=0.95,size=1)
savp(10,8, 'C://Users//andrew.caudillo.BUGLAB-I9//Box//NAMC//Research Projects//AIM//IncreaserDecreaser_OE//sites_traitspacetop5.png')
vegan::adonis2(sit_gow~status,data=BC_scores)
bd<-vegan::betadisper(sit_gow,BC_scores$status)
anova(bd)
}
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

