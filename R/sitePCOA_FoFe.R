#PcoA for sites using traits as variables
#this chunk is readind in the data,
#subseting, etc
# top5s=c("Emerge_synch_abbrev_Poorly",  "Habit_prim_Swimmer",          "Crawl_rate_high",
#         "Voltinism_abbrev_Univoltine", "Feed_prim_abbrev_PR")
if(0){
abun<-read.csv("C://Users//andrew.caudillo.BUGLAB-I9//Box//NAMC//Research Projects//AIM//IncreaserDecreaser_OE//abundances.csv")
abun$sample_id<-as.character(abun$sample_id)
abun<-abun[,-1]
}
odd_sites_samps=c("173029", "184847", "185060", "187112", "187223", "187301", "187304",
                  "187318", "187780", "190728","190768")
refs=Os
Probs=ProbOs
refs$status='Reference'
Probs$status='Probabilistic'
PA=rbind(refs,Probs)#PA<-ifelse(abun>0,1,0)
failed_sites<-read.csv('C://Users//andrew.caudillo.BUGLAB-I9//Box//NAMC//Research Projects//AIM//IncreaserDecreaser_OE//failed_sites.csv')
PA=PA[row.names(PA) %in% failed_sites$sampleId==F,]
PA=PA[row.names(PA) %in% odd_sites_samps==F,]

#ref_a<-abun[abun$Status=='Reference',]
#ref_a<-ref_a[ref_a$sample_id %in% c(116830,132243,132807, 185022)==F,]
#prob_a<-abun[abun$Status=='Probabilistic',]
#prob_a<-prob_a[prob_a$sample_id %in% failed_sites$sampleId==F,]
# traits<-read.csv('C://Users//andrew.caudillo.BUGLAB-I9//Box//NAMC//Research Projects//AIM//IncreaserDecreaser_OE//Updated_w_modelObj//OTU_21_traits.csv')
# OTU21<-read.csv('C://Users//andrew.caudillo.BUGLAB-I9//Box//NAMC//Research Projects//AIM//IncreaserDecreaser_OE//OTU21.csv')
traits=maxxed_OTUs
traits=traits[traits$OTU %in% c('Farula','Stactobiella')==F,]
#traits<-traits[traits$OTU %in% OTU21$taxonOTU,]
traits=traits[!duplicated(traits$OTU),]
PA_all=PA
PA=PA[,names(PA) %in% c('status',traits$OTU)]
PA_all=PA
ref_a<-PA[PA$status=='Reference',]

library(tidyverse)

taxa_traits_num=traits

PA_long=PA %>%
  tibble::rownames_to_column(var='sample_id') %>%
  pivot_longer(cols=-c(sample_id,status),
                    names_to = 'OTU',values_to = 'PA')
joined<-plyr::join(PA_long,taxa_traits_num,type='left')


#set the trait columns
#for dplyr to calculate the weighted traits
trait_cols=5:ncol(joined)
trait_cols=names(joined)[trait_cols]

library(dplyr)

#row.names(taxa_traits_num)<-taxa_traits_num$OTU;taxa_traits_num=taxa_traits_num[,names(taxa_traits_num) %in% 'OTU'==F]


 site_PA_traits=joined %>%
   group_by(sample_id, status) %>%
   dplyr::summarise(

     across(
       all_of(trait_cols),
       ~ as.integer(any(.x == 1 & PA == 1, na.rm = TRUE)),
       .names = "PA_{.col}"
     ),

     .groups = "drop"
   )



site_PA_traits<-as.data.frame(site_PA_traits)

#assign a status to each set of sites

site_PA_traits$status<-PA$status



site_PA_NMDS=vegan::vegdist(site_PA_traits[,3:ncol(site_PA_traits)],'bray',binary=T)


Vegdist_PA<-vegan::vegdist(PA[,-ncol(PA)],'bray',binary=T)
#vegdist_PA=vegan::vegdist(PA_all[,-ncol(PA_all)],'bray',binary=T)
Veg_PA_mat<-as.matrix(Vegdist_PA)
#assign the site names
row.names(Veg_PA_mat)<-row.names(PA)
#row.names(Veg_PA_mat)<-row.names(PA_all)
set.seed(99)
PA_MDS=vegan::metaMDS(Vegdist_PA,k=3)
#abun_scores=as.data.frame(vegan::scores(Gower_abun_MDS))
PA_scores=as.data.frame(vegan::scores(PA_MDS))
row.names(PA_scores)<-row.names(PA)
#PA_scores$Status=PA_all$status
PA_scores$Status=PA$status

NMDS_cal_preds=read.csv("C://Users//andrew.caudillo.BUGLAB-I9//Box//NAMC//Research Projects//AIM//IncreaserDecreaser_OE//Updated_w_modelObj//Cal_predictors.csv")
NMDS_prob_preds=read.csv("C://Users//andrew.caudillo.BUGLAB-I9//Box//NAMC//Research Projects//AIM//IncreaserDecreaser_OE//Updated_w_modelObj//Prob_predictors_prednew.csv")
prob_comids=read.csv('C://Users//andrew.caudillo.BUGLAB-I9//Box//NAMC//Research Projects//AIM//IncreaserDecreaser_OE//prob_comids.csv')
#NMDS_prob_preds$sampleId=row.names(NMDS_prob_preds)

NMDS_prob_preds=NMDS_prob_preds[NMDS_prob_preds$sampleId %in% failed_sites$sampleId==F,]
NMDS_prob_preds=NMDS_prob_preds[NMDS_prob_preds$sampleId %in% row.names(PA),]
NMDS_prob_preds=plyr::join(NMDS_prob_preds,prob_comids,by='sampleId')

#row.names(NMDS_prob_preds)<-NMDS_prob_preds$sampleId

#row.names(NMDS_prob_preds)==row.names(PA)[PA$status=='Probabilistic']
cal_pred_dat=read.csv('C://Users//andrew.caudillo.BUGLAB-I9//Box//NAMC//Research Projects//AIM//IncreaserDecreaser_OE//Updated_w_modelObj//WW_cal_preds.csv')
NMDS_cal_preds<-plyr::join(NMDS_cal_preds,cal_pred_dat,by='COMID','left')
names(NMDS_cal_preds)[2]='sampleId'
NMDS_cal_preds=NMDS_cal_preds[NMDS_cal_preds$sampleId %in% row.names(PA),]
env_dat=rbind(NMDS_cal_preds[,c('sampleId','ElevCat','Precip8110Ws','Tmean8110Ws','WsAreaSqKm','COMID')],NMDS_prob_preds[,c('sampleId','ElevCat','Precip8110Ws','Tmean8110Ws','WsAreaSqKm','COMID')])
AgUrb_dat=read.csv('C://Users//andrew.caudillo.BUGLAB-I9//Box//NAMC//Research Projects//AIM//IncreaserDecreaser_OE//all_site_UrbAg.csv')
env_dat=plyr::join(env_dat,AgUrb_dat[,c('COMID','PCtCrop','AvgUrb')],by='COMID','left')

env_dat=env_dat[env_dat$sampleId %in% row.names(PA),]


env_dat3=env_dat[env_dat$sampleId %in% row.names(PA_scores),]
row.names(PA_scores)==env_dat3$sampleId


#env_dat4=env_dat3[env_dat3$sampleId %in% row.names(PA_scores),]
set.seed(99)
nmds_fit=envfit(PA_scores[,1:2],env_dat3[,!names(env_dat3)%in% c('COMID','sampleId')])
# trait fitting

# nmds_fitA=envfit(PA_scores[,c(1,2)],scale(env_dat3[,names(env_dat3) %in% c('ElevCat','Precip8110Ws',
#                                                                   'Tmean8110Ws','WsAreaSqKm',
#                                                                   'PCtCrop',
#                                                                   'AvgUrb')]),na.rm=T)
env_vect<-as.data.frame(scores(nmds_fit,display = 'vectors'))

env_vect$env <- rownames(env_vect)
colnames(env_vect)[1:2] <- c("NMDS1", "NMDS2")
#mult=vegan::ordiArrowMul(nmds_fit)
#PA_scores2=PA_scores[PA_scores$NMDS1<34 & PA_scores$NMDS2 < 2,]
#env_vect[,1:2] <- env_vect[,1:2] /
 # sqrt(rowSums(env_vect[,1:2]^2))
env_vect[,1:2] <- env_vect[,1:2]*1.5
#env_vect[,1:2] <- env_vect[,1:2] * mult
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
      env == "Mean Precip"   ~0.5,
      env== 'Mean Temp' ~ -0.4,
      TRUE                   ~  0
    ),

    nudge_y = case_when(
      env=='% Urban' ~ 0.2,
      env=='% Mean Temp' ~ 0.5,
      env == "Catchment Elev"~ 0.2,
      env == "% Crops" ~ 0.028,
      env == 'Mean Precip' ~ -0.5,
      TRUE             ~ 0
    )
  )


## Traits in OTU space

all_site_PA_traits=read.csv('C://Users//andrew.caudillo.BUGLAB-I9//Box//NAMC//Research Projects//AIM//IncreaserDecreaser_OE//Updated_w_modelObj//all_site_trait_PAs.csv')
all_site_PA_traits=all_site_PA_traits[,1:26]
row.names(all_site_PA_traits)<-all_site_PA_traits$sampleId
all_site_PA_traits=all_site_PA_traits[,-ncol(all_site_PA_traits)]
row.names(all_site_PA_traits)==row.names(PA_scores)

all_site_PA_traits=all_site_PA_traits[row.names(PA_scores),]
set.seed(99)
trait_fit=envfit(PA_scores[,c(1,2)],all_site_PA_traits)

env_vect<-as.data.frame(scores(trait_fit,display = 'vectors'))

env_vect$env <- rownames(env_vect)
colnames(env_vect)[1:2] <- c("NMDS1", "NMDS2")

Fisher_Traits=c('Rheophily_abbrev_depo',
                'Habit_prim_Burrower',
                'Habit_prim_Climber',
                'Habit_prim_Swimmer',
                'Feed_prim_abbrev_CG',
                'Lifespan_very_short',
                'Emerge_synch_abbrev_Poorly',
                'Habit_prim_Sprawler',
                'Female_disp_abbrev_High')
Fishtrait_vect=env_vect[ env_vect$env %in% Fisher_Traits,]
env_vect=env_vect[env_vect$env %in% names(top5),]
#ptrait_vect$trait=row.names(ptrait_vect)
Fishtrait_vect$trait=c('Poor Synch. Emerge',
                       'Collector\nGatherer',
                       'Burrower',
                       'Climber',
                       'Sprawler',
                       'Swimmer',
                       'Rheophilic',
                       'Short Life',
                       'Female Dispersal'

)

Fishtrait_vect
#ptrait_vec=as.data.frame(scores(env_fut,display='vectors'))
#ptrait_vect[,1:2]<-ptrait_vect[,1:2]*-1
# env_vect$env=c('Poor Synch.\nEmerge',
#                     'Clinger',
#                     'Slow\nDevelopment',
#                     'Short life',
#                     'No Drift')
env_vect$env=c('Strong Flyer',
               'Poor Synch.\nEmerge',
               'Swimmer',
               'Slow Development',
               'Crawler')
#mult=vegan::ordiArrowMul(nmds_fit)
#PA_scores2=PA_scores[PA_scores$NMDS1<34 & PA_scores$NMDS2 < 2,]
#env_vect[,1:2] <- env_vect[,1:2] /
# sqrt(rowSums(env_vect[,1:2]^2))
env_vect[,1:2] <- env_vect[,1:2]*1.5
#env_vect[,1:2] <- env_vect[,1:2] * mult
# env_vect$env=c('Catchment Elev',
#                'Mean Precip',
#                'Mean Temp','Basin Area',
#                '% Crops',
#                '% Urban')

# env_vect <- env_vect %>%
#   mutate(
#     nudge_x = case_when(
#       env == "Basin Area"    ~  0.1,
#       env == "Catchment Elev"~ 0.1,
#       env == "% Urban"       ~  0.1,
#       env == "Mean Precip"   ~0.5,
#       env== 'Mean Temp' ~ -0.4,
#       TRUE                   ~  0
#     ),
#
#     nudge_y = case_when(
#       env=='% Urban' ~ 0.2,
#       env=='% Mean Temp' ~ 0.5,
#       env == "Catchment Elev"~ 0.2,
#       env == "% Crops" ~ 0.028,
#       env == 'Mean Precip' ~ -0.5,
#       TRUE             ~ 0
#     )
#   )
# ggplot(abun_scores,aes(x=NMDS1,y=NMDS3,fill=Status))+geom_point(pch=21)+
#   scale_fill_manual(name='Status',
#                     values=c('Reference' = 'purple4',
#                              'Probabilistic' = 'yellow'))+ggtitle('Sites in OTU space (PA)')+
#   stat_ellipse(level=0.8,data = subset(abun_scores, Status == "Reference"), color = "purple4", size = 1,type='t') +
#   stat_ellipse(level=0.8,data = subset(abun_scores, Status == "Probabilistic"), color = "yellow", size = 1) #+stat_ellipse(level=0.8,type='t')+scale_color_manual(values=c('Reference'='red','Probabilistic' = 'blue'))
# savp(10,8, 'C://Users//andrew.caudillo.BUGLAB-I9//Box//NAMC//Research Projects//AIM//IncreaserDecreaser_OE//sites_OTUspace_abun_ellipse.png')
#windows(10,8)

ggplot(PA_scores[PA_scores$NMDS2 < 2 & PA_scores$NMDS1 < 30,], aes(x=NMDS1, y=NMDS2, fill=Status)) +
  geom_point(pch=21, size=3) +
  scale_fill_manual(name='Status',
                    values=c('Reference' = 'blue',
                             'Probabilistic' = 'orange3')) +
  geom_segment(data=Fishtrait_vect,
               aes(x=0, y=0, xend=NMDS1, yend=NMDS2),
               linewidth=1,
               alpha=0.9, inherit.aes = FALSE,
               arrow=arrow(length=unit(.25, "centimeters"))) +
  geom_text_repel(
    data = Fishtrait_vect,
    aes(x = NMDS1, y = NMDS2, label = trait),
    size = 6,
    inherit.aes = FALSE,
    box.padding = 0.2,
    point.padding = 0.1,
    min.segment.length = 0,
    max.overlaps = Inf
  ) +
  stat_ellipse(level=0.8, data = subset(PA_scores, Status == "Reference"), color = "blue", size = 1, type='t') +
  stat_ellipse(level=0.8, data = subset(PA_scores, Status == "Probabilistic"), color = "orange3", size = 1) +
  theme_classic() +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.title = element_text(size = 25),

    # FIX APPLIED HERE: Added ggplot2:: prefix
    axis.title.x = element_text(margin = ggplot2::margin(t = 20)),
    axis.title.y = element_text(margin = ggplot2::margin(r = 20)),

    strip.text = element_text(size = 14),
    legend.text = element_text(size = 22),
    legend.title = element_text(size = 22),
    legend.key.size = unit(1, "lines"),
    axis.ticks = element_blank(),
    legend.position = 'inside',
    legend.position.inside = c(0.85, 0.2)
  )


savp(10,8, 'C://Users//andrew.caudillo.BUGLAB-I9//Box//NAMC//Research Projects//AIM//IncreaserDecreaser_OE//Updated_w_modelObj//sites_OTUspace_NMDS_Fishtraits260518.png')

#
# ordisurf(PA_scores[,c(1,2)], env_dat$PCtCrop,k=4,
#          pch=21,bg=ifelse(PA$status=='Reference','purple4','yellow'),
#          col='black',
#          main='PctCrop',
#          lwd=3)
# legend('topleft',
#        pch=rep(21,2),
#        pt.bg=c('purple4','yellow'),
#        bty='n',
#        leg=c('Reference',
#              'Probabilistic'),
#        cex=0.8)

savp(10,8,'C://Users//andrew.caudillo.BUGLAB-I9//Box//NAMC//Research Projects//AIM//IncreaserDecreaser_OE//Updated_w_modelObj//PA_ordi_Crop.png')

set.seed(99)
site_PA_traitsMDS=vegan::metaMDS(site_PA_traits[,3:ncol(site_PA_traits)],binary=T)
site_PA_traitscores=as.data.frame(scores(site_PA_traitsMDS)$sites)

#gower_scores=data.frame(Pco1=gower_cmd[,1],Pco2=gower_cmd[,2])
site_PA_tait_Weighted_Scores$sample_id=site_PA_traits_weighted$sample_id
set.seed(99)
#
# ptrait_vect<-as.data.frame(env_fut$vectors$arrows)
# ptrait_vect$trait=row.names(ptrait_vect)
# Fishtrait_vect=ptrait_vect[substr(ptrait_vect$trait,4,nchar(ptrait_vect$trait)) %in% Fisher_Traits,]
# ptrait_vect=ptrait_vect[substr(ptrait_vect$trait,4,nchar(ptrait_vect$trait)) %in% names(top5),]
# #ptrait_vect$trait=row.names(ptrait_vect)
# Fishtrait_vect$trait=c('Collector\nGatherer',
#                        'Female Dispersal',
#                        'Burrower',
#                        'Climber',
#                        'Sprawler',
#                        'Large Body',
#                        'Slow\nDevelopment',
#                        'Rare Drifter'
#
# )
#
# Fishtrait_vect
# ptrait_vec=as.data.frame(scores(env_fut,display='vectors'))
# ptrait_vect$trait <- rownames(ptrait_vect)
# colnames(ptrait_vect)[1:2] <- c("NMDS1", "NMDS2")
# ptrait_vect=ptrait_vect[ptrait_vect$trait %in% paste0('PA_',names(top5)),]
# #ptrait_vect[,1:2]<-ptrait_vect[,1:2]*-1
# ptrait_vect$trait=c('Poor Synch.\nEmerge',
#                     'Clinger',
#                     'Slow\nDevelopment',
#                     'Short life',
#                     'Streamlined')
#
# p_sites_and_coords=plyr::join(site_PA_traits_weighted,site_PA_tait_Weighted_Scores)
#
#
#
# point_max <- max(abs(site_PA_tait_Weighted_Scores[,1:2]))
# arrow_max <- max(abs(site_PA_tait_Weighted_Scores[,1:2]))
# arrow_multiplier <- (point_max / arrow_max) * 0.8 # 0.8 keeps them inside the points
#
# # Apply scaling
# ptrait_vect$NMDS1 <- ptrait_vect$NMDS1 * arrow_multiplier
# ptrait_vect$NMDS2<- ptrait_vect$NMDS2 * arrow_multiplier
# names(p_sites_and_coords)[2]<-'Status'
#
# nudge_x_vector <- ifelse(trait_vect$trait == 'Slow\nDevelopment', 0.3,
#                          ifelse(trait_vect$trait == 'Multivoltine/\nBivoltine', -0.4,
#                                 ifelse(trait_vect$trait=='Poor Synch.\nEmerge',-0.1,
#                                        ifelse(trait_vect$trait=='Swimmer',-0.4,
#                                               0))))
# nudge_y_vector <- ifelse(trait_vect$trait == 'Slow\nDevelopment', 0.2,
#                          ifelse(trait_vect$trait == 'Multivoltine/\nBivoltine', -0.7,
#                                 ifelse(trait_vect$trait=='Poor Synch.\nEmerge',0.4,
#                                        ifelse(trait_vect$trait == 'Univoltine',-0.1,
#                                               ifelse(trait_vect$trait=='Swimmer',-0.2,0)))))
# Fishnudge_x_vector <- ifelse(Fishtrait_vect$trait == 'Slow\nDevelopment', 0.5,
#                              ifelse(Fishtrait_vect$trait == 'Female Dispersal', -0.5,
#                                     ifelse(Fishtrait_vect$trait=='Large Body',0.4,
#                                            ifelse(Fishtrait_vect$trait=='Climber',0.5,
#                                                   ifelse(Fishtrait_vect$trait=='Collector\nGatherer',-0.2,
#                                                          ifelse(Fishtrait_vect$trait=='Burrower',-0.3,
#                                                                 ifelse(Fishtrait_vect$trait=='Rare Drifter',0.5,0)))))))
# Fishnudge_y_vector <- ifelse(Fishtrait_vect$trait == 'Sprawler', -0.1,
#                              ifelse(Fishtrait_vect$trait == 'Female Dispersal', -0.2,
#                                     ifelse(Fishtrait_vect$trait=='Collector\nGatherer',0.9,
#                                            ifelse(Fishtrait_vect$trait=='Gills',-0.5,
#                                                   ifelse(Fishtrait_vect$trait=='Slow\nDevelopment',0.5,
#                                                          ifelse(Fishtrait_vect$trait=='Burrower',0.5,0))))))
#
#
# # --- 2. Plotting with corrected scales ---
# ggplot(p_sites_and_coords, aes(x=NMDS1, y=NMDS2, fill=Status)) +
#   geom_point(pch=21,size=3) +
#   stat_ellipse(level=0.8,data = subset(p_sites_and_coords, Status == "Reference"), color = "blue", size = 1,type='t') +
#   stat_ellipse(level=0.8,data = subset(p_sites_and_coords, Status == "Probabilistic"), color = "orange3", size = 1,type='t')+# Added alpha for fill
#   # FIX: Change fill to color to match aes(color=status)
#   scale_fill_manual(name='Status',
#                      values=c('Reference' = 'blue',
#                               'Probabilistic' = 'orange3')) +
#   geom_segment(data=Fishtrait_vect,
#                aes(x=0, y=0, xend=NMDS1, yend=NMDS2),
#                arrow=arrow(length=unit(0.25, 'cm')),
#                linewidth=1,
#                inherit.aes = FALSE,
#                color="black") +
#
  #  geom_text_repel(
  #    data =Fishtrait_vect,
  #    aes(
  #      x = NMDS1,
  #      y = NMDS2,
  #      label = trait
  #    ),
  #    # Use your logic to "nudge" the starting position
  # #   nudge_y = Fishnudge_y_vector,
  # #   nudge_x = Fishnudge_x_vector,
  #   #                  ifelse(Fishtrait_vect$trait == 'Swimmer', 0.5, ifelse(Fishtrait_vect$vect=='Poor Synch. Emerge',0.8,0))),
  #   segment.color = "black", # Hides the segments entirely
  #    size = 6,
  #    inherit.aes = FALSE,
  #    box.padding = 0.2,    # Extra "breathing room" around labels
  #    point.padding = 0.1,  # Distance from the data point
  #    min.segment.length = 0, # Always draw lines to labels if they move far
  #    max.overlaps = Inf    # Forces all labels to show
  # )+# Set a fixed color so it doesn't look for 'status'
  # theme_classic()+
  # labs(x='NMDS1',y='NMDS2')+
  # theme(#axis.text.x = element_text(size = 14),
  #       axis.ticks = element_blank(),
  #       axis.text=element_blank(),
  #     #  axis.text.y = element_text(size = 14),
  #       axis.title = element_text(size = 25),
  #       legend.position = 'inside',
  #       legend.position.inside = c(0.15,0.9),
  #     axis.title.x = element_text(margin=margin(t=20)),
  #     axis.title.y  = element_text(margin=margin(r=20)),
  #     legend.text = element_text(size = 22),
  #     legend.title = element_text(size = 22))#sures 1 unit on X = 1 unit on Y
# savp(10,8, 'C://Users//andrew.caudillo.BUGLAB-I9//Box//NAMC//Research Projects//AIM//IncreaserDecreaser_OE//Updated_w_modelObj//sites_traitspace_PA_ellipse_Fishtraits_axes1_2_260515_nolabs.png')
#
#
#
# set.seed(99)
# PA_trait_fit=envfit(ord=site_PA_traitscores[,c(1,2)],env=site_PA_traits[,3:ncol(site_PA_traits)],eig=T)
# site_PA_traitscores$sample_id=site_PA_traits$sample_id
#
#
#
#
# PA_trait_vect<-as.data.frame(env_fut$vectors$arrows)
# PA_trait_vect$trait=row.names(PA_trait_vect)
# Fishtrait_vect=PA_trait_vect[substr(PA_trait_vect$trait,4,nchar(PA_trait_vect$trait)) %in% Fisher_Traits,]
# PA_trait_vect=PA_trait_vect[substr(PA_trait_vect$trait,4,nchar(PA_trait_vect$trait)) %in% names(top5),]
# #PA_trait_vect$trait=row.names(PA_trait_vect)
# Fishtrait_vect$trait=c('Collector\nGatherer',
#                        'Female Dispersal',
#                        'Burrower',
#                        'Climber',
#                        'Sprawler',
#                        'Large Body',
#                        'Slow\nDevelopment',
#                        'Rare Drifter'
#
# )
#
# Fishtrait_vect
# ptrait_vec=as.data.frame(scores(env_fut,display='vectors'))
# PA_trait_vect$trait <- rownames(PA_trait_vect)
# colnames(PA_trait_vect)[1:2] <- c("NMDS1", "NMDS2")
# PA_trait_vect=PA_trait_vect[PA_trait_vect$trait %in% paste0('PA_',names(top5)),]
# #PA_trait_vect[,1:2]<-PA_trait_vect[,1:2]*-1
# PA_trait_vect$trait=c('Poor Synch.\nEmerge',
#                     'Clinger',
#                     'Slow\nDevelopment',
#                     'Short life',
#                     'Streamlined')
#
# pa_sites_and_coords=plyr::join(site_PA_traits,site_PA_traitscores)
# names(pa_sites_and_coords)[2]<-'Status'
#
#
# point_max <- max(abs(site_PA_tait_Weighted_Scores[,1:2]))
# arrow_max <- max(abs(site_PA_tait_Weighted_Scores[,1:2]))
# arrow_multiplier <- (point_max / arrow_max) * 0.8 # 0.8 keeps them inside the points
#
# # Apply scaling
# PA_trait_vect$NMDS1 <- PA_trait_vect$NMDS1 * arrow_multiplier
# PA_trait_vect$NMDS2<- PA_trait_vect$NMDS2 * arrow_multiplier
# names(p_sites_and_coords)[2]<-'Status'
#
# nudge_x_vector <- ifelse(trait_vect$trait == 'Slow\nDevelopment', 0.3,
#                          ifelse(trait_vect$trait == 'Multivoltine/\nBivoltine', -0.4,
#                                 ifelse(trait_vect$trait=='Poor Synch.\nEmerge',-0.1,
#                                        ifelse(trait_vect$trait=='Swimmer',-0.4,
#                                               0))))
# nudge_y_vector <- ifelse(trait_vect$trait == 'Slow\nDevelopment', 0.2,
#                          ifelse(trait_vect$trait == 'Multivoltine/\nBivoltine', -0.7,
#                                 ifelse(trait_vect$trait=='Poor Synch.\nEmerge',0.4,
#                                        ifelse(trait_vect$trait == 'Univoltine',-0.1,
#                                               ifelse(trait_vect$trait=='Swimmer',-0.2,0)))))
# Fishnudge_x_vector <- ifelse(Fishtrait_vect$trait == 'Slow\nDevelopment', 0.5,
#                              ifelse(Fishtrait_vect$trait == 'Female Dispersal', -0.5,
#                                     ifelse(Fishtrait_vect$trait=='Large Body',0.4,
#                                            ifelse(Fishtrait_vect$trait=='Climber',0.5,
#                                                   ifelse(Fishtrait_vect$trait=='Collector\nGatherer',-0.2,
#                                                          ifelse(Fishtrait_vect$trait=='Burrower',-0.3,
#                                                                 ifelse(Fishtrait_vect$trait=='Rare Drifter',0.5,0)))))))
# Fishnudge_y_vector <- ifelse(Fishtrait_vect$trait == 'Sprawler', -0.1,
#                              ifelse(Fishtrait_vect$trait == 'Female Dispersal', -0.2,
#                                     ifelse(Fishtrait_vect$trait=='Collector\nGatherer',0.9,
#                                            ifelse(Fishtrait_vect$trait=='Gills',-0.5,
#                                                   ifelse(Fishtrait_vect$trait=='Slow\nDevelopment',0.5,
#                                                          ifelse(Fishtrait_vect$trait=='Burrower',0.5,0))))))
#
#
# # --- 2. Plotting with corrected scales ---
# ggplot(pa_sites_and_coords, aes(x=NMDS1, y=NMDS2, fill=Status)) +
#   geom_point(pch=21,size=3) +
#   stat_ellipse(level=0.8,data = subset(p_sites_and_coords, Status == "Reference"), color = "blue", size = 1,type='t') +
#   stat_ellipse(level=0.8,data = subset(p_sites_and_coords, Status == "Probabilistic"), color = "orange3", size = 1,type='t')+# Added alpha for fill
#   # FIX: Change fill to color to match aes(color=status)
#   scale_fill_manual(name='Status',
#                     values=c('Reference' = 'blue',
#                              'Probabilistic' = 'orange3')) +
#   geom_segment(data=PA_trait_vect,
#                aes(x=0, y=0, xend=NMDS1, yend=NMDS2),
#                arrow=arrow(length=unit(0.25, 'cm')),
#                linewidth=1,
#                inherit.aes = FALSE,
#                color="black") +
#
#   #  geom_text_repel(
#   #    data =PA_trait_vect,
#   #    aes(
#   #      x = NMDS1,
#   #      y = NMDS2,
#   #      label = trait
#   #    ),
#   #    # Use your logic to "nudge" the starting position
#   # #   nudge_y = Fishnudge_y_vector,
#   # #   nudge_x = Fishnudge_x_vector,
#   #   #                  ifelse(Fishtrait_vect$trait == 'Swimmer', 0.5, ifelse(Fishtrait_vect$vect=='Poor Synch. Emerge',0.8,0))),
#   #   segment.color = "black", # Hides the segments entirely
#   #    size = 6,
#   #    inherit.aes = FALSE,
#   #    box.padding = 0.2,    # Extra "breathing room" around labels
#   #    point.padding = 0.1,  # Distance from the data point
#   #    min.segment.length = 0, # Always draw lines to labels if they move far
#   #    max.overlaps = Inf    # Forces all labels to show
#   # )+# Set a fixed color so it doesn't look for 'status'
#   theme_classic()+
#   labs(x='NMDS1',y='NMDS2')+
#   theme(#axis.text.x = element_text(size = 14),
#     axis.ticks = element_blank(),
#     axis.text=element_blank(),
#     #  axis.text.y = element_text(size = 14),
#     axis.title = element_text(size = 25),
#     legend.position = 'inside',
#     legend.position.inside = c(0.15,0.9),
#     axis.title.x = element_text(margin=margin(t=20)),
#     axis.title.y  = element_text(margin=margin(r=20)),
#     legend.text = element_text(size = 22),
#     legend.title = element_text(size = 22))#sures 1 unit on X = 1 unit on Y
# savp(10,8, 'C://Users//andrew.caudillo.BUGLAB-I9//Box//NAMC//Research Projects//AIM//IncreaserDecreaser_OE//Updated_w_modelObj//sites_traitspace_PresAbs_ellipse_traits_axes1_2_260515_nolabs.png')
