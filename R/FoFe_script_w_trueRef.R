##### Increaser and Decreaser O/E using the WestWide model
#Executed by Andrew Caudillo
#Process developed by Charles P. Hawkins and Jennifer Courtwright
#along with Hawkins and Yuan and others.

#This script will:
#1 assign all the sites used in this assessment to an ecoregion
#that is set by the Bureau of Land Management, BLM (specifically,
#Western Rivers and Streams Assessment, WRSA)

#2. Assign a dominant land use to each site via StreamCat
#3. Determine the response to change (Increaser or Decreaser)
#For each taxon used in the WestWide model
#4. Join the Increaser/Decreaser table to a taxon/trait table
#for further inspection

#load package for regular expression usage and
#setting column order in datasets
rm(list=ls())
library(data.table)
library(ggrepel)
library(dplyr)
library(tidyverse)
library(tidyr)
savp = function(W,H,fn) {
  dev.copy(dev=png,file=fn,wi=W,he=H,un="in",res=650)
  dev.off()
}
val_bugs=read.csv('C://Users//andrew.caudillo.BUGLAB-I9//Box//NAMC//Research Projects//AIM//IncreaserDecreaser_OE//Updated_w_modelObj//WW_validation_data.csv')
row.names(val_bugs)<-val_bugs$Site;val_bugs=val_bugs[,-c(1,2)]
val_preds=read.csv('C://Users//andrew.caudillo.BUGLAB-I9//Box//NAMC//Research Projects//AIM//IncreaserDecreaser_OE//Updated_w_modelObj//val_preds.csv')
row.names(val_preds)<-val_preds$SiteCode;val_preds=val_preds[,-c(1,2)]
val_Pcs=read.csv('C://Users//andrew.caudillo.BUGLAB-I9//Box//NAMC//Research Projects//AIM//IncreaserDecreaser_OE//Updated_w_modelObj//Val_Pcs_summedTany.csv')
row.names(val_Pcs)=val_Pcs$sampleId;val_Pcs=val_Pcs[,-c(1,2)]
#taxa_notraits_rare=read.csv('C://Users//andrew.caudillo.BUGLAB-I9//Box//NAMC//Research Projects//AIM//IncreaserDecreaser_OE//taxa_to_drop.csv')
#nonrare=read.csv('C://Users//andrew.caudillo.BUGLAB-I9//Box//NAMC//Research Projects//AIM//IncreaserDecreaser_OE//nonrare_taxa.csv')

Os<-read.csv("C://Users//andrew.caudillo.BUGLAB-I9//Box//NAMC//Research Projects//AIM//IncreaserDecreaser_OE//Updated_w_modelObj//WestWide_Ref_PA_260209.csv")
row.names(Os)<-Os$X;Os=Os[,-1]
Os=Os[,names(Os) %in% c('Carabidae','Curculionidae')==F]
Os=Os%>%
  mutate(Tanypodinae = ifelse(Tanypodinae+Tanypodinae.1>1,1,Tanypodinae+Tanypodinae.1)) %>%
  select(-Tanypodinae.1)
Es<-read.csv("C://Users//andrew.caudillo.BUGLAB-I9//Box//NAMC//Research Projects//AIM//IncreaserDecreaser_OE//Updated_w_modelObj//Ref_Pcs_orig_WW_Tany_summed.csv")
row.names(Es)<-Es$sampleId;Es=Es[,-c(1,2)]
Es=Es[,names(Es) %in% c('Carabidae','Curculionidae')==F]
#ensure that the column order matches across the tables
#i.e., the taxa names
Os=Os[,names(Es)]
Os=Os[row.names(Os) %in% row.names(Es),]
#probabilistic predictions
Pes<-read.csv("C://Users//andrew.caudillo.BUGLAB-I9//Box//NAMC//Research Projects//AIM//IncreaserDecreaser_OE//Updated_w_modelObj//Prob_Pcs_orig_WW_droppedTany.csv")
row.names(Pes)=Pes$sampleId;Pes=Pes[,-c(1,2)]

#FailedPes=read.csv('C://Users//andrew.caudillo.BUGLAB-I9//Box//NAMC//Research Projects//AIM//IncreaserDecreaser_OE//Updated_w_modelObj//WestWide_failed_Pcs.csv')

ProbOs<-read.csv('C://Users//andrew.caudillo.BUGLAB-I9//Box//NAMC//Research Projects//AIM//IncreaserDecreaser_OE//MRF_Prob_Os.csv')
row.names(ProbOs)<-ProbOs$X;ProbOs=ProbOs[,-c(1,2)]


ProbOs=ProbOs[,names(ProbOs) %in% c('Carabidae','Curculionidae')==F]
Pes=Pes[,names(Pes) %in% c('Carabidae','Curculionidae')==F]
Pes=Pes[match(names(ProbOs), names(Pes))]
#write.csv(ProbOs,'C://Users//andrew.caudillo.BUGLAB-I9//Box//NAMC//Research Projects//AIM//IncreaserDecreaser_OE//MRF_Prob_Os.csv')

failed_sites<-read.csv('C://Users//andrew.caudillo.BUGLAB-I9//Box//NAMC//Research Projects//AIM//IncreaserDecreaser_OE//failed_sites.csv')
FailedO=ProbOs[row.names(ProbOs) %in% failed_sites$sampleId,]
FailedPes=Pes[row.names(Pes) %in% failed_sites$sampleId,]
ProbOs=ProbOs[row.names(ProbOs) %in% failed_sites$sampleId==F,]
Pes=Pes[row.names(Pes) %in% failed_sites$sampleId==F,]
#dropping sampleId and assigning that to the row names




N_Tax=ncol(Es)
refE<-colSums(Es)
ProbE<-colSums(Pes)
refO<-colSums(Os)
ProbO<-colSums(ProbOs)
refFoFe=refO / refE
probFoFe= ProbO / ProbE
FailedE<-colSums(FailedPes)
FailedO<-colSums(FailedO)
FailedE=FailedE[match(names(FailedO), names(FailedE))]
failedFoFe=FailedO/FailedE


Ref_Results=data.frame(taxon=names(Os),
                       Fo=refO,
                       Fe=refE,
                       ratio = refFoFe)



P_Results=data.frame(taxon=names(Os),
                     Fo=ProbO,
                     Fe=ProbE,
                     ratio = probFoFe)

F_Results=data.frame(taxon=names(FailedO),
                     Fo=FailedO,
                     Fe=FailedE,
                     ratio=failedFoFe)
P_Results$taxon[P_Results$taxon=='Drunella_coloradensis_flavilinea']<-'Drunella coloradensis/flavilinea'
F_Results$taxon[F_Results$taxon=='Drunella_coloradensis_flavilinea']<-'Drunella coloradensis/flavilinea'
#a few reference taxa are increasers, but this could be due to the cutoff.
#The 4 taxa (Cinygma, Hesperoconopa, Limnophila, and Paraperla) have a ratio
# > 1.2 but < 1.3.
Ref_Results$Regional_Response=ifelse(Ref_Results$ratio >1.2, "Increaser",
                                     ifelse( Ref_Results$ratio > 0.8 & Ref_Results$ratio < 1.2, 'Neutral',"Decreaser"))
Ref_Results$status='Reference'
P_Results$Regional_Response=ifelse(P_Results$ratio >1.2, "Increaser",
                                   ifelse( P_Results$ratio > 0.8 & P_Results$ratio < 1.2, 'Neutral',"Decreaser"))
P_Results$status='Probabilistic'
F_Results$Regional_Response=ifelse(F_Results$ratio >1.2, "Increaser",
                                   ifelse( F_Results$ratio > 0.8 & F_Results$ratio < 1.2, 'Neutral',"Decreaser"))
F_Results$status='Failed'
All_Results<-rbind(Ref_Results,P_Results,F_Results)
P_Results$diff <- P_Results$Fe - P_Results$Fo #make a difference column
# bottom_right <- head(P_Results[order(-P_Results$diff), ], 3) #get the two most decreasesrs
# top_left <- head(P_Results[order(P_Results$diff), ], 3) #get the two most inc
extreme_diffs <- c('Callibaetis','Laccobius','Paracloeodes','Cambaridae','Sciomyzidae')
P_Results$col=ifelse(P_Results$taxon %in% extreme_diffs,
                     'red',NA)
highlight_dat=P_Results[which(P_Results$col=='red'),]
my_colors <- setNames(
  c('red', 'purple', 'orange', 'dodgerblue','blue','yellow2'),
  extreme_diffs[1:6]
)
# F_Results$diff <- F_Results$Fe - F_Results$Fo #make a difference column
# Fbottom_right <- head(F_Results[order(-F_Results$diff), ], 3) #get the two most decreasesrs
# Ftop_left <- head(F_Results[order(F_Results$diff), ], 3) #get the two most inc
# Fextreme_diffs <- rbind(Ftop_left, Fbottom_right)
# F_Results$col=ifelse(F_Results$taxon %in% c(extreme_diffs$taxon),
#                      'red',NA)
# Fhighlight_dat=F_Results[which(F_Results$col=='red'),]
# my_F_colors <- setNames(
#   c('red', 'purple', 'orange', 'dodgerblue','blue','yellow2'),
#   extreme_diffs$taxon[1:6]
# )
P=ggplot(data=Ref_Results,aes(y=Fo,x=Fe))+geom_point()+
  geom_abline(intercept = 0,slope = 1,col='red')+
  annotate("text",
           x = -Inf, y = Inf, # Position at top-left corner
           label = 'Reference',
           hjust = 0, vjust = 1, # Justify text relative to corner
           size = 5, color = "black")+
  lims(x=c(0,350),y=c(0,max(P_Results$Fo)))
Fa=ggplot(data=F_Results,aes(y=Fo,x=Fe,label=taxon))+geom_point()+
  #geom_point(data=P_Results[P_Results$Regional_Response=='Neutral',],color='grey50')+
  #geom_text_repel(box.padding = 0.5, max.overlaps = Inf) +
  #geom_point(data = P_Results[P_Results$Regional_Response=='Decreaser',],color='dodgerblue')+
  #geom_point(data = P_Results[P_Results$Regional_Response=='Increaser',],color='orange2')+
  geom_abline(intercept = 0,slope = 1,col='red')+
  geom_point(data=F_Results[F_Results$taxon %in% extreme_diffs,],aes(x=Fe,y=Fo,colour = taxon))+
  scale_color_manual(values=my_colors)+
  coord_cartesian(clip = "off") +
  # geom_text_repel(
  #                  point.size=NA)+

  annotate("text",
           x = -Inf, y = Inf, # Position at top-left corner
           label = 'Failed',
           hjust = 0, vjust = 1, # Justify text relative to corner
           size = 5, color = "black")+
  theme(legend.position = "bottom")

gridExtra::grid.arrange(P,Fa,ncol=1)

savp(10,8, 'C://Users//andrew.caudillo.BUGLAB-I9//Box//NAMC//Research Projects//AIM//IncreaserDecreaser_OE//Updated_w_modelObj//Fo_vs_Fe_failed_colored_260218.png')

write.csv(All_Results,'C://Users//andrew.caudillo.BUGLAB-I9//Box//NAMC//Research Projects//AIM//IncreaserDecreaser_OE//Updated_w_modelObj//regional_responses_WW_260223.csv')

#assigning sites to EcoRegions
#using the shapefile on NAMC's GIS folder,
#perform a spatial join
if(0){
  EcoRegions=sf::st_read("C://Users//andrew.caudillo.BUGLAB-I9//Box//NAMC WATS Department Files//GIS//GIS_Stats//CONUS//ecoregion//hybrid_level_III_ecoregions.shp")
  site_coords<-read.csv('C://Users//andrew.caudillo.BUGLAB-I9//Box//NAMC//Research Projects//AIM//IncreaserDecreaser_OE//site_coordinates.csv')
  #setting up points with coordinate system
  site_coords<-sf::st_as_sf(site_coords,coords=c('lon','lat'),crs=4269)
  site_coords<-sf::st_transform(site_coords,sf::st_crs(EcoRegions))
  #the spatial join
  site_coords_EcoR<-sf::st_intersection(site_coords,EcoRegions)
  #check to see if there NAs, etc.
  table(site_coords_EcoR$EPA_hybird)

  groups<-split(site_coords_EcoR,site_coords_EcoR$EPA_hybird)
  groups<-do.call(rbind, groups)
  groups$EPA_hybird
  groups$EPA_hybird=ifelse(groups$EPA_hybird=='Eastern Xeric Plateaus',
                           'Eastern Xeric Plateaus','Other')
  #StreamCat Site Land Use binning
  #using waterbody code from NAMC's database,
  #assign a dominant land use type (i.e., anthropogenic use)
  #to each site via the highest percentage in its watershed
  #since NAMC's database only has StreamCat data up to 2016,
  #the StreamCat API was used for this task, which has data up to 2019.
  #2019 anthropogenic values were used for this assessment
  #read in the output of the API
  sites<-read.csv('C://Users//andrew.caudillo.BUGLAB-I9//Box//NAMC//Research Projects//AIM//IncreaserDecreaser_OE//COMIDS_and_samples.csv')
  sites<-sites[sites$sampleId %in% c(116830,132243,132807,185022)==F,]
  SCDat=read.csv('C://Users//andrew.caudillo.BUGLAB-I9//Box//NAMC//Research Projects//AIM//IncreaserDecreaser_OE//SC_Dat_anthro_natural.csv')
  SCDat<-plyr::join(SCDat,sites,by='COMID','right')

  SCDat<-SCDat[SCDat$sampleId %in% row.names(Pes),]
  SCDat$Rangeland=SCDat$pctshrb2019ws
  SCDat$Urban=rowSums(SCDat[,9:12])
  SCDat$Agriculture=rowSums(SCDat[,7:8])
  #determined cutoffs based on distribution
  #of StreamCat metrics
  #these BLM sites are not really urbanized or agriculture influenced
  #so seeing a lot of Mixed makes sense.
  classify_watershed = function(urban, ag, grass){
    if (any(is.na(c(urban, ag, grass)))) return(NA_character_)

    # Priority order: Urban > Agriculture > Rangeland > Mixed
    if (urban > 10 & ag < 20) {
      return("Urban")
    } else if (urban < 10 & ag > 20) {
      return("Agriculture")
    } else if (grass > 66) {
      return("Rangeland")
    } else {
      return("Mixed")
    }
  }





  SCDat$DomLandUse= mapply(classify_watershed,
                           grass=SCDat$Rangeland,
                           ag=SCDat$Agriculture,
                           urban=SCDat$Urban)
  #read in the COMID/sampleID matrix
  AIM_human=read.csv('C://Users//andrew.caudillo.BUGLAB-I9//Box//NAMC//Research Projects//AIM//IncreaserDecreaser_OE//AIM_human_influence.csv')
  NRSA_human=read.csv('C://Users//andrew.caudillo.BUGLAB-I9//Box//NAMC//Research Projects//AIM//IncreaserDecreaser_OE//NRSA_human_influence.csv')

  all_human_influence=rbind(AIM_human,NRSA_human)
  all_human_influence<-all_human_influence[all_human_influence$sampleId %in% c(116830,132243,132807, 185022)==F,]
  #join the two datasets
  sites[duplicated(sites$sampleId),]
  sites<-plyr::join(sites,all_human_influence,by='sampleId','left')

  #function to rescale the weighted values from 0 to 1 for easy comparison
  #regardless of weight
  normalize_component <- function(x) x / max(x, na.rm = TRUE)
  #normalize the values in a new column for each
  sites$LoggingOperations_rescale=normalize_component(sites$LoggingOperations)
  sites$OilGas_rescale=normalize_component(sites$OilGas)
  sites$PastureHayFence_rescale=normalize_component(sites$PastureHayFence)
  #Classify the local land use based on proportion
  classify_local <- function(p) {
    if (is.na(p)) return(NA)
    if (p > 2/3) return("high")
    if (p > 0.50) return("moderate")
    if (p > 1/3) return("low")
    return("negligible")
  }
  #now classify the sites using the function
  sites$LoggingOperations_class<-sapply(sites$LoggingOperations_rescale,classify_local)
  sites$OilGas_class<-sapply(sites$OilGas_rescale,classify_local)
  sites$PastureHayFence_class<-sapply(sites$PastureHayFence_rescale,classify_local)
  get_stressor_label <- function(classes) {
    # classes is a named vector, e.g. c(grazing="high", oilgas=NA, urban="low")
    classes <- classes[!is.na(classes)]  # drop any NAs

    # If no non-NA classes
    if (length(classes) == 0) return("none")

    # Priority order
    if (any(classes == "high")) {
      highs <- names(classes[classes == "high"])
      return(if(length(highs) == 1) paste('high_',highs )else "mixed")
    }

    if (any(classes == "moderate")) {
      mods <- names(classes[classes == "moderate"])
      return(if(length(mods) == 1) paste('moderate_',mods) else "mixed")
    }

    lows <- names(classes[classes == "low"])
    if (length(lows) >= 2) return("mixed")
    if (length(lows) == 1) return(paste('low_',lows))

    return("none")  # if only negligible or nothing left
  }
  #create a column that is the dominant local stressor
  sites <-sites %>%
    rowwise() %>%
    mutate(dominant_field_land_use = get_stressor_label(c(
      Grazing=PastureHayFence_class,
      Oilgas=OilGas_class,
      Logging=LoggingOperations_class
    )))
  #joining all site infor together now
  sites<-plyr::join(sites[sites$sampleId %in% SCDat$sampleId,],SCDat,by='sampleId','left')
  sites<-plyr::join(sites,sf::st_drop_geometry(site_coords_EcoR[,c('sampleId','EPA_hybird')]),by='sampleId','left')
  #assign the sampleIds as the row names, then drop
  #sampleID and COMID, as they are no longer needed.
  row.names(sites)<-sites$sampleId;sites<-sites[,-c(1,2)]

  #create a land use data frame for lookup if needed
  Land_use<-data.frame(sampleId=row.names(sites),
                       EcoRegion=sites$EPA_hybird,
                       Land_use=sites$DomLandUse,
                       Local_land_use=sites$dominant_field_land_use)
  write.csv(Land_use,'C://Users//andrew.caudillo.BUGLAB-I9//Box//NAMC//Research Projects//AIM//IncreaserDecreaser_OE//site_LandUses_251107.csv',row.names = F)
  Land_use<-read.csv('C://Users//andrew.caudillo.BUGLAB-I9//Box//NAMC//Research Projects//AIM//IncreaserDecreaser_OE//site_LandUses_250903.csv')

  #getting land use-specific ratios

  #There are no reference sites within region 1, so ratio will always be NaN.


  #this section calcualtes the Pc for each probabilistic and
  #reference site for each ecoregion.
  #It also assigns land use to sites within each ecoregion
}
Xeric_P=c(171583, 171664, 171668, 171669, 171777, 173398, 173404,
          173407, 173408, 173409, 178667, 178671, 178676, 184764,
          184768, 185060, 185063, 185699, 185702, 187162, 187182,
          187607,187662, 187779, 187780, 187849, 190728, 190762,
          190916, 191507,
          #these are the reference sites used in the westwide model
          "87087", "HAWK-216", "HAWK-225", "HAWK-226", "HAWK-227", "HAWK-228", "HAWK-241",
          "HAWK-86",  "HAWK-89", "HAWK-92",  "WE-1036",  "WE-1047",  "WE-1055",  "WE-1079",
          "WE-1085",  "WE-21",    "WE-47",    "WE-866",  "WE-889",   "WE-905")
if(1){
  #add another section for Prob Os
  Pesgrp2<-Pes[row.names(Pes) %in%Xeric_P,]
  Pesgrp2$Ecoregion='Eastern Xeric Plateaus'
  Pesgrp3<-Pes[row.names(Pes) %in% Xeric_P==F,]
  Pesgrp3$Ecoregion='Other'


  ProbOsgrp2<-ProbOs[row.names(ProbOs) %in% Xeric_P,]
  ProbOsgrp2$Ecoregion='Eastern Xeric Plateaus'
  ProbOsgrp3<-ProbOs[row.names(ProbOs) %in% Xeric_P==F,]
  ProbOsgrp3$Ecoregion='Other'


  Osgrp2<-Os[row.names(Os) %in% Xeric_P,]
  Osgrp2$Ecoregion='Eastern Xeric Plateaus'
  Osgrp3<-Os[row.names(Os) %in% Xeric_P==F,]
  Osgrp3$Ecoregion='Other'

  Esgrp2<-Es[row.names(Es) %in% Xeric_P,]

  Esgrp3<-Es[row.names(Es) %in% Xeric_P==F,]




  L_grp2<-list(ProbOsgrp2,Pesgrp2)
  L_grp3<-list(ProbOsgrp3,Pesgrp3)



  Grp2_ratio=colSums(L_grp2[[1]][1:ncol(L_grp2[[1]])-1])/colSums((L_grp2[[2]][1:ncol(L_grp2[[2]])-1]))
  Grp3_ratio=colSums(L_grp3[[1]][1:ncol(L_grp3[[1]])-1])/colSums((L_grp2[[2]][1:ncol(L_grp2[[2]])-1]))


  Land_use_tidy<-Land_use
  names(Land_use_tidy)[1]<-'ID'

  P2_tidy <- Pesgrp2 %>%
    rownames_to_column("ID")
  P3_tidy <- Pesgrp3 %>%
    rownames_to_column("ID")

  Po2_tidy <- ProbOsgrp2 %>%
    rownames_to_column("ID")
  Po3_tidy <- ProbOsgrp3 %>%
    rownames_to_column("ID")




  P2=plyr::join(P2_tidy,Land_use_tidy[,c('Land_use','Local_land_use','ID')],by='ID')
  P3=plyr::join(P3_tidy,Land_use_tidy[,c('Land_use','Local_land_use','ID')],by='ID')


  Po2=plyr::join(Po2_tidy,Land_use_tidy[,c('Land_use','Local_land_use','ID')],by='ID')
  Po3=plyr::join(Po3_tidy,Land_use_tidy[,c('Land_use','Local_land_use','ID')],by='ID')



  #create final table of probabilistic sites
  All_Prob_E<-do.call(rbind, list(P2,P3))#P4,P5,P6,P7,P8,P9,P10))
  All_Prob_O<-do.call(rbind,list(Po2,Po3))#Po4,Po5,Po6,Po7,Po8,Po9,Po10))
  All_Prob_E$Status='PExpected'
  All_Prob_O$Status='PObserved'
  All_Prob_O=All_Prob_O[All_Prob_O$ID!=185022,]
  All_Prob_E=All_Prob_E[All_Prob_E$ID!=185022,]

  All_Prob_sites<-rbind(All_Prob_O,All_Prob_E)
  All_Prob_sites<-All_Prob_sites[All_Prob_sites$ID %in% failed_sites$sampleId==F,]

  # Get column names of taxa
  taxa_cols <- colnames(All_Prob_sites)[2:(ncol(All_Prob_sites)-4)]

  # Aggregate by ecoregion and land use
  n_sites_tbl= All_Prob_sites %>%
    group_by(Ecoregion, Land_use) %>%
    dplyr::summarise(n_sites =n(), .groups='drop')

  n_sites_local_tbl<-All_Prob_sites %>%
    group_by(Ecoregion, Local_land_use) %>%
    dplyr::summarise(n_sites =n(), .groups='drop')
  #create final table of all reference sites
  Esgrp2$Ecoregion=Osgrp2$Ecoregion
  Esgrp3$Ecoregion=Osgrp3$Ecoregion



  E2_tidy <- Esgrp2 %>%
    rownames_to_column("ID")
  E3_tidy <- Esgrp3 %>%
    rownames_to_column("ID")

  O2_tidy <- Osgrp2 %>%
    rownames_to_column("ID")
  O3_tidy <- Osgrp3 %>%
    rownames_to_column("ID")

  All_Os<-do.call(rbind,list(O2_tidy,O3_tidy))#O4_tidy,O5_tidy,O6_tidy,O7_tidy,O8_tidy,O9_tidy,O10_tidy))
  All_Os$Status='rObserved'
  All_ref=do.call(rbind, list(E2_tidy,E3_tidy))#,E4_tidy,E5_tidy,E6_tidy,E7_tidy,E8_tidy,E9_tidy,E10_tidy))
  All_ref$Status='rExpected'

  All_ref=rbind(All_Os,All_ref)


  All_data=dplyr::bind_rows(All_ref,All_Prob_sites)

  Os_lookup<-data.frame(ID=All_Os$ID,
                        O=NA)

  otaxa_start <- 2
  otaxa_end <- ncol(All_Os) - 3  # because last two columns are ecoregion & land_use



  row.names(All_Os)<-All_Os$ID;All_Os<-All_Os[,-1]

  #get taxonomy columns
  #we are now going to combine all the data into the final table!
  #the final table will show land use, ecoregion, etc.
  #along with the response of the taxa.

  taxa_start=2
  #this can change based on how many
  #extra columns are tacked onto the end
  taxa_end=ncol(All_data)-4
  taxa_end=ncol(All_Prob_sites)-4
  #change our wide data into long format


  df_long <- All_data %>%
    pivot_longer(
      cols=taxa_start:taxa_end,
      names_to='taxon',
      values_to='value'
    )



  pFo= aggregate(
    value ~ Ecoregion + taxon,
    data = df_long[df_long$Status == "PObserved", ],
    FUN = sum,
    na.rm = TRUE
  )
  pFe= aggregate(
    value ~ Ecoregion + taxon,
    data = df_long[df_long$Status == "PExpected", ],
    FUN = sum,
    na.rm = TRUE
  )
  rFo=aggregate(
    value ~ Ecoregion + taxon,
    data = df_long[df_long$Status == "rObserved", ],
    FUN = sum,
    na.rm = TRUE
  )

  rFe=aggregate(
    value ~ Ecoregion + taxon,
    data = df_long[df_long$Status == "rExpected", ],
    FUN = sum,
    na.rm = TRUE
  )

  Ref_Ecoregions=plyr::join(rFo,rFe,by=c('Ecoregion','taxon'))
  names(Ref_Ecoregions)[3:4]<-c('Fo','Fe')
  Ref_Ecoregions$Rratio = Ref_Ecoregions$Fo / Ref_Ecoregions$Fe
  Ref_Ecoregions$Rresponse=ifelse(Ref_Ecoregions$Rratio > 1.2,'Increaser',
                                  ifelse(Ref_Ecoregions$Rratio > 0.8 & Ref_Ecoregions$Rratio < 1.2, 'Neutral',
                                         'Decreaser'))

  P_Ecoregions=plyr::join(pFo,pFe,by=c('Ecoregion','taxon'))
  names(P_Ecoregions)[3:4]<-c('pFo','pFe')
  P_Ecoregions$Pratio = P_Ecoregions$pFo / P_Ecoregions$pFe
  P_Ecoregions$Presponse=ifelse(P_Ecoregions$Pratio > 1.2,'Increaser',
                                ifelse(P_Ecoregions$Pratio > 0.8 & P_Ecoregions$Pratio < 1.2, 'Neutral',
                                       'Decreaser'))

  All_ratios=Ref_Ecoregions %>% left_join(P_Ecoregions, by=c('Ecoregion','taxon'))
  write.csv(All_ratios,'C://Users//andrew.caudillo.BUGLAB-I9//Box//NAMC//Research Projects//AIM//IncreaserDecreaser_OE//Updated_w_modelObj//all_ratios_260210.csv')
}
#Trait info
#read in the traits from NAMC's database
#translationID=21
#contains bench ID and OTU
#this needs to only be run once, as the table does not change.
if(0){
  trait_table=read.csv('C://Users//andrew.caudillo.BUGLAB-I9//Box//NAMC//Research Projects//AIM//IncreaserDecreaser_OE//Updated_w_modelObj//OTU_21_traits.csv')
  #trait_table<-read.csv('C://Users//andrew.caudillo.BUGLAB-I9//Box//NAMC//Research Projects//AIM//IncreaserDecreaser_OE//All_bench_taxa_atrributes2.csv')
  OTU21<-read.csv('C://Users//andrew.caudillo.BUGLAB-I9//Box//NAMC//Research Projects//AIM//IncreaserDecreaser_OE//OTU21.csv')
  trait_table<-trait_table[trait_table$OTU %in% OTU21$taxonOTU,]
  #trait_table<-trait_table[trait_table$OTU %in% All_ratios$taxon,]
  trait_table<-trait_table[trait_table$OTU %in% unique(All_Results$taxon),]
  #subset the OTUs and their mapping taxa
  OTUs<-trait_table[,1:2]
  #remove duplicated records for joining later
  OTUs<-OTUs[!duplicated(OTUs),]
  #drop the OTU column from the trait table
  #trait_table<-trait_table[,-2]
  #make the trait table wide format (since the O/E table is wide)
  wide_traits<-as.data.frame(trait_table |>
                               pivot_wider(names_from = attribute_name, values_from = attribute_value,
                                           values_fn = list) |>
                               unnest(cols = everything() ))
  #add the OTU names to the dataset
  wide_trait_table_OTUs<-plyr::join(wide_traits,OTUs)
  #rename "OTU" to be "Taxon" so we can join the tables
  names(wide_trait_table_OTUs)[names(wide_trait_table_OTUs)=='OTU']<-'taxon'
  wide_trait_table_OTUs<- wide_trait_table_OTUs %>%
    select(-matches(c("_af",'HBI','_sec_','_PH')))
  #forcing everything to be binary
  wide_trait_table_OTUs[wide_trait_table_OTUs=='FALSE']<-0
  wide_trait_table_OTUs[wide_trait_table_OTUs=='TRUE']<-1

  #this yields a many to one join, so there will be duplicates.
  ratios_w_traits=plyr::join(All_Results,wide_trait_table_OTUs,by='taxon','left')
  #ratios_w_traits=plyr::join(All_ratios,wide_trait_table_OTUs,by='taxon','left')
  #creating filter column, removing duplicates, and then dropping that column.
  ratios_w_traits<-ratios_w_traits[,names(ratios_w_traits) %in% 'scientific_name'==F]
  #ratios_w_traits<-ratios_w_traits[!duplicated(ratios_w_traits$taxon),]
  #dropping taxa that only occur at 10 or more sites

write.csv(ratios_w_traits,'C://Users//andrew.caudillo.BUGLAB-I9//Box//NAMC//Research Projects//AIM//IncreaserDecreaser_OE//Updated_w_modelObj//Regional_responses_and_traits_260210.csv',row.names = F)
  write.csv(ratios_w_traits,'C://Users//andrew.caudillo.BUGLAB-I9//Box//NAMC//Research Projects//AIM//IncreaserDecreaser_OE//Spatial_scales_and_traits_MRF_251114.csv',row.names = F)
}

#read in table that came from the trait section
finaltable3<-read.csv('C://Users//andrew.caudillo.BUGLAB-I9//Box//NAMC//Research Projects//AIM//IncreaserDecreaser_OE//Updated_w_modelObj//Regional_responses_and_traits_260210.csv')
#final_table3<-read.csv('C://Users//andrew.caudillo.BUGLAB-I9//Box//NAMC//Research Projects//AIM//IncreaserDecreaser_OE//Spatial_scales_and_traits_MRF_251114.csv')

a<-ggplot(All_Results[All_Results$status=='Reference',],aes(x=Fe,y=Fo))+geom_point()+geom_abline(color='red',linetype=2)+ggtitle('Reference')+labs(y='Fo',x='Fe')
b<-ggplot(P_Results[P_Results$status=='Probabilistic',],aes(x=Fe,y=Fo))+geom_point()+geom_abline(color='red',linetype=2)+ggtitle('Probabilistic')+labs(y='Fo',x='Fe')
gridExtra::grid.arrange(a,b,ncol=1)
savp(10,8, 'C://Users//andrew.caudillo.BUGLAB-I9//Box//NAMC//Research Projects//AIM//IncreaserDecreaser_OE//FovsFe_MRF.png')
boxplot(All_Results$ratio[All_Results$status=='Reference'],at=1,xlim=c(0,4),ylim=c(0,40),col='purple',ylab='Fo/Fe ratio')
boxplot(All_Results$ratio[All_Results$status=='Probabilistic'],at=2,add=T,col='yellow3')
boxplot(All_Results$ratio[All_Results$status=='Failed'],at=3,add=T,col='red')
Failed_to_hi=All_Results[which(All_Results$status=='Failed' & All_Results$ratio>=10),]
Failed_to_hi<-Failed_to_hi[is.finite(Failed_to_hi$ratio),]
Failed_to_hi<-Failed_to_hi[order(Failed_to_hi$ratio, decreasing = F),]
points(x=rep(3,nrow(Failed_to_hi)-1),y=Failed_to_hi$ratio[1:(nrow(Failed_to_hi)-1)],pch=21,bg=c('red','blue',
                                                                                                'black','orange',
                                                                                                'purple','dodgerblue'))
legend('topleft',
       leg=c('Menetus',
             'Cambaridae',
             'Stenelmis',
             'Erioptera',
             'Crangonyx',
             'Caecidotea',
             'Paracloeodes'),
       pt.bg=c('red','blue','black','orange',
               'purple','dodgerblue'),
       bty='n',
       pch=rep(21,6),
       ncol=3,
       cex=0.8)
axis(side=1, labels = c('Ref.','Prob','Failed'),at=c(1,2,3))
savp(10,8, 'C://Users//andrew.caudillo.BUGLAB-I9//Box//NAMC//Research Projects//AIM//IncreaserDecreaser_OE//Ref_Prob_Failed_box251125_Colored.png')

par(mfrow=c(2,1))
hist(log(All_Results$Fo[All_Results$status=='Reference'] - All_Results$Fe[All_Results$status=='Reference']), xlab='',main='Ref.')

hist(log(All_Results$Fo[All_Results$status=='Probabilistic'] - All_Results$Fe[All_Results$status=='Probabilistic']),xlab='log(Fo-Fe)',main='Prob.')
savp(10,8, 'C://Users//andrew.caudillo.BUGLAB-I9//Box//NAMC//Research Projects//AIM//IncreaserDecreaser_OE//Ref_Prob_hist_MRF.png')
#these are taxa for which we have essentially 0 metrics for, so we cannot use them in the
#Dissimilarity matrix, for they will result in NAs in the final product.



traits_only=finaltable3[,c(1,7:ncol(finaltable3))]
#traits_only <- final_table3[,c(2,13:ncol(final_table3))]



traits_only<-traits_only[!duplicated(traits_only$taxon),]
row.names(traits_only)<-traits_only$taxon
traits_only<-traits_only[,-1]


traits_only2<-as.data.frame(lapply(traits_only,factor))
#traits_only2<-as.data.frame(lapply(traits_only[,2:ncol(traits_only)],factor))
traits_num<-as.data.frame(lapply(traits_only,as.numeric))
row.names(traits_only2)<-row.names(traits_only)
# Step 2: compute correlation matrix
cor_mat <- cor(traits_num, use = "pairwise.complete.obs")

find_redundant_pairs <- function(cor_mat, threshold = 0.9) {
  redundant <- which(abs(cor_mat) > threshold & abs(cor_mat) < 1, arr.ind = TRUE)
  pairs <- data.frame(
    trait1 = rownames(cor_mat)[redundant[,1]],
    trait2 = colnames(cor_mat)[redundant[,2]],
    correlation = cor_mat[redundant]
  ) %>%
    filter(trait1 < trait2)  # avoid duplicates (A,B) and (B,A)
  return(pairs)
}

redundant_pairs <- find_redundant_pairs(cor_mat, threshold = 0.8)

redundant_pairs
#removing redundant traits that confuse Gower
#binary variables = mutually exclusive traits
#must remove those traits before running a dissimilarity distance
traits_only2<-traits_only2[,names(traits_only2) %in% c('Develop_fast_season',
                                                       'Resp_abbrev_Gills',
                                                       'Occurence_drift_common',
                                                       'Swim_strong',
                                                       'Attach_free_range',
                                                       'Attach_both',
                                                       'Female_disp_abbrev_High','Swim_none',
                                                       'Attach_free_range_both',
                                                       'Emerge_season_1_Spring',
                                                       'Emerge_season_2_Summer'
                                                       )==F]

#calculate gower distance on just traits
traits_only3<- lapply(traits_only2[,2:ncol(traits_only2)], \(x) as.integer(as.factor(x)) - 1L)
traits_only3<-as.data.frame(do.call(cbind,traits_only3))
traits_only3[is.na(traits_only3)]<-0
# traits_imputed <- traits_only3 %>%
#   mutate(across(where(is.integer), ~replace_na(.,0,)))
# traits_imputed_num <- traits_only2 %>%
#   mutate(across(everything(), ~ {
#     # 1. Convert factor to character, then to numeric
#     # 2. Replace NA with 0
#     x <- as.numeric(as.character(.))
#     replace_na(x, 0)
#   }))
# traits_imputed<-as.data.frame(lapply(traits_imputed,factor))
# traits_imputed_num=lapply(traits_imputed, \(x) as.integer(as.factor(x)) - 1L)
# traits_num=as.data.frame(do.call(cbind,traits_imputed_num))
Gower_dist<-FD::gowdis(traits_only3)

library(vegan)
library(gawdis)
#this is the top variable identification, modified from Bello et al. 2020
#this gets the PCoA
pcoaaxes=dudi.pco(cailliez(Gower_dist), scannf = F,nf=267) #all axes
#and its distance
gowdis.PCoA<-dist(pcoaaxes$li)

#now we are comparing the gower distance to the of the trait to the PCoA
cors.pcoa<-vector()
for(i in 1:dim(traits_only3)[2]){
  cors.pcoa[i]<-mantel(gowdis.PCoA, gowdis(traits_only3[,i,drop=F]))$statistic
}

names(cors.pcoa)<-names(traits_only3)
ranked_mantel=sort(cors.pcoa,decreasing = T)
top5=ranked_mantel[1:5]
top5

#convert the Gower distance into a matrix for readability
Gower_mat<-as.matrix(Gower_dist)
#assign the taxa names
row.names(Gower_mat)<-row.names(traits_only3)
#perform PCoA
PCOA<-ape::pcoa(Gower_mat)
#scree plot. Axes 1 and 2 explain ~50% of variation
eigenvalues=PCOA$values$Eigenvalues
prop_var=eigenvalues/sum(eigenvalues)*100
cum_var=cumsum(prop_var)
barplot(prop_var[1:40],
        names.arg = 1:40,
        main = "PCoA Scree Plot: Variance Explained per Axis",
        xlab = "PCoA Axis",
        ylab = "Variance Explained (%)",
        col = "steelblue",
        las = 2)
savp(10,8,'C://Users//andrew.caudillo.BUGLAB-I9//Box//NAMC//Research Projects//AIM//IncreaserDecreaser_OE//PCoA_Scree.png')
#cmd scale traits for layering?
gower_cmd=cmdscale(Gower_dist,45,eig=T)

gower_scores=data.frame(PCOA$vectors)
#gower_scores=data.frame(Pco1=gower_cmd[,1],Pco2=gower_cmd[,2])
gower_scores$taxon=row.names(Gower_dist)
row.names(gower_scores)=row.names(traits_only3)

fit=envfit(gower_scores[,c(1,2)],traits_only3,na.rm=T)
trait_vect<-as.data.frame(fit$vectors$arrows)
trait_vec=as.data.frame(scores(fit, display = 'vectors'))
trait_vect$trait <- rownames(trait_vect)
colnames(trait_vect)[1:2] <- c("PCo1", "PCo2")
mult=vegan::ordiArrowMul(trait_vec,display='species')
arrow_multiplier <- 0.75 * max(abs(gower_scores[,c(1,2)])) / max(abs(trait_vect[,1:2]))
trait_vect[,1:2] <- trait_vect[,1:2] * arrow_multiplier
trait_vect=trait_vect[trait_vect$trait %in% names(top5),]
trait_vect[,1:2]<-trait_vect[,1:2]*-1

PCOA_Scores<-as.data.frame(PCOA$vectors)
#assign taxa names
top_trait_arrows=trait_vect$trait
PCOA_Scores$taxon<-row.names(traits_only)
#renaming vectors for easier interpretation

trait_vect$trait=c('Fast develop',
                   'Adults enter/exit',
                   'Drift',
                   'Desicc. Resist.',
                   'Streamlined')




#join the responses with the PCoA scores
Responses_and_scores=plyr::join(All_Results,PCOA_Scores,by='taxon','left')
All_responses_and_scores=plyr::join(All_ratios,PCOA_Scores,by='taxon','left')
#plotting the responses by site type
#be sure to adjust axes as needed
Responses_and_scores$status2=factor(Responses_and_scores$status,levels=c('Reference','Probabilistic'))
ggplot(data=Responses_and_scores[Responses_and_scores$status!='Failed' & is.finite(Responses_and_scores$ratio) & Responses_and_scores$Fo>0,],
       aes(x=Axis.1, y=Axis.2,color=Regional_Response))+
  geom_point(size=3,alpha=0.7)+
  stat_ellipse(level=0.8)+
  geom_segment(data=trait_vect,
               aes(x=0,y=0,xend=PCo1,yend=PCo5),
               arrow=arrow(length=unit(0.25,'cm')),
               linewidth=1,
               alpha=0.5,inherit.aes = F)+
  # geom_text(data = trait_vect,
  #          aes(x = PCo1, y = PCo2, label = trait),
  #         hjust = 0.4, vjust = .9, color="black", size = 3.5) +
  geom_label(
    data = trait_vect,
    aes(x = PCo1, y = PCo2, label = trait),
    hjust = ifelse(names(trait_vect)[1]=='PCo2',0.2, .75),
    vjust = 1,
    size = 3,
    fill = "white",      # background color
    color = "red",       # text color
    label.size = 0.2,    # border thickness; set to 0 to remove outline
    inherit.aes = FALSE
  )+
    # geom_text(data=Responses_and_scores[Responses_and_scores$status!='Failed' & is.finite(Responses_and_scores$ratio),],
    #          aes(x=Axis.2,y=Axis.3,label=taxon),
    #           color='black',
    #           hjust=0,vjust=0, size=1)+
  scale_color_manual(values=c("Increaser"="orange",
                              "Decreaser"="dodgerblue3",
                              "Neutral"='red',
                              "Not expected"='black'))+
  facet_wrap('status2')+guides(color=guide_legend(title='Response'))+
  theme(legend.position = "top")
#guides(color = guide_legend(nrow = 1))

savp(10,8,'C://Users//andrew.caudillo.BUGLAB-I9//Box//NAMC//Research Projects//AIM//IncreaserDecreaser_OE//Updated_w_modelObj//regional_responses_traitspace_260219_Axes2_3_flipped.png')


#same but for prob sites.
if(0){
ggplot(All_responses_and_scores[is.finite(All_responses_and_scores$Pratio) & All_responses_and_scores$pFo > 0,],
       aes(x = Axis.1, y = Axis.2, color = Presponse)) +
  geom_point(size = 3, alpha = 0.7) +
  stat_ellipse(level = 0.8) +
  scale_color_manual(values = c(
    "Increaser" = "orange",
    "Decreaser" = "dodgerblue3",
    "Neutral" = "red",
    "Not expected" = "black"
  )) +
  facet_wrap(~Ecoregion) +
  guides(color = guide_legend(title = 'Response'))+
  theme(legend.position = "top")+
  geom_segment(data=trait_vect,
               aes(x=0,y=0,xend=PCo1,yend=PCo2),#                arrow=arrow(length=unit(0.25,'cm')),
               linewidth=1,
               alpha=0.5,inherit.aes = F, arrow=arrow(length=unit(.25,"centimeters")))+
  # geom_text(data = trait_vect,
  #           aes(x = PCo1, y = PCo2, label = trait),
  #           hjust = 0.5, vjust = .5, color="black", size = 3.5)
  geom_label(
    data = trait_vect,
    aes(x = PCo1, y = PCo2, label = trait),
    hjust = ifelse(names(trait_vect)[1]=='PCo2',0.2, .56),
    vjust = 1,
    size = 3,
    fill = "white",      # background color
    color = "red",       # text color
    label.size = 0.2,    # border thickness; set to 0 to remove outline
    inherit.aes = FALSE
  )
savp(10,8,'C://Users//andrew.caudillo.BUGLAB-I9//Box//NAMC//Research Projects//AIM//IncreaserDecreaser_OE//prob_responses_traitspace_Ecoregions_260210_Axes1_2_allaxesmantel.png')
}


A<-ggplot(All_responses_and_scores[is.finite(All_responses_and_scores$Pratio) & All_responses_and_scores$pFo > 0,],
          aes(x = Axis.1, y = Axis.2, color = Presponse)) +
  geom_point(size = 3, alpha = 0.7) +
  stat_ellipse(level = 0.8) +
  scale_color_manual(values = c(
    "Increaser" = "orange",
    "Decreaser" = "dodgerblue3",
    "Neutral" = "red",
    "Not expected" = "black"
  )) +
  facet_wrap(~Ecoregion) +
  guides(color = guide_legend(title = 'Response',ncol=3))+
  geom_segment(data=trait_vect,
               aes(x=0,y=0,xend=PCo1,yend=PCo2),#                arrow=arrow(length=unit(0.25,'cm')),
               linewidth=1,
               alpha=0.5,inherit.aes = F, arrow=arrow(length=unit(.25,"centimeters")))+
  # geom_text(data = trait_vect,
  #           aes(x = PCo1, y = PCo1, label = trait),
  #           hjust = 0.5, vjust = .5, color="black", size = 3.5)
  geom_label(
    data = trait_vect,
    aes(x = PCo1, y = PCo2, label = trait),
    hjust = ifelse(names(trait_vect)[1]=='PCo1',0.2, .56),
    vjust = 1,
    size = 3,
    fill = "white",      # background color
    color = "red",       # text color
    label.size = 0.2,    # border thickness; set to 0 to remove outline
    inherit.aes = FALSE
  )


B<-ggplot(All_responses_and_scores[is.finite(All_responses_and_scores$Rratio) & All_responses_and_scores$Fo > 0,],
          aes(x = Axis.1, y = Axis.2, color = Rresponse)) +
  geom_point(size = 3, alpha = 0.7) +
  stat_ellipse(level = 0.8) +
  scale_color_manual(values = c(
    "Increaser" = "orange",
    "Decreaser" = "dodgerblue3",
    "Neutral" = "red",
    "Not expected" = "black"
  )) +
  facet_wrap(~Ecoregion) +
  guides(color = guide_legend(title = 'Response'))+
  geom_segment(data=trait_vect,
               aes(x=0,y=0,xend=PCo1,yend=PCo2),#                arrow=arrow(length=unit(0.25,'cm')),
               linewidth=1,
               alpha=0.5,inherit.aes = F, arrow=arrow(length=unit(.25,"centimeters")))+
  # geom_text(data = trait_vect,
  #           aes(x = PCo1, y = PCo1, label = trait),
  #           hjust = 0.5, vjust = .5, color="black", size = 3.5)
  geom_label(
    data = trait_vect,
    aes(x = PCo1, y = PCo2, label = trait),
    hjust = ifelse(names(trait_vect)[1]=='PCo1',0.2, .56),
    vjust = 1,
    size = 3,
    fill = "white",      # background color
    color = "red",       # text color
    label.size = 0.2,    # border thickness; set to 0 to remove outline
    inherit.aes = FALSE
  )+theme(legend.position = 'none')+
  annotate("text",
           x = -Inf, y =Inf, # Position at top-left corner
           label = 'Reference',
           hjust = 0, vjust = 1, # Justify text relative to corner
           size = 3, color = "black")
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
tax_legend=get_only_legend(A)
A<-ggplot(All_responses_and_scores[is.finite(All_responses_and_scores$Pratio) & All_responses_and_scores$pFo > 0,],
          aes(x = Axis.1, y = Axis.2, color = Presponse)) +
  geom_point(size = 3, alpha = 0.7) +
  stat_ellipse(level = 0.8) +
  scale_color_manual(values = c(
    "Increaser" = "orange",
    "Decreaser" = "dodgerblue3",
    "Neutral" = "red",
    "Not expected" = "black"
  )) +
  facet_wrap(~Ecoregion) +
  guides(color = guide_legend(title = 'Response'))+
  geom_segment(data=trait_vect,
               aes(x=0,y=0,xend=PCo1,yend=PCo2),#                arrow=arrow(length=unit(0.25,'cm')),
               linewidth=1,
               alpha=0.5,inherit.aes = F, arrow=arrow(length=unit(.25,"centimeters")))+
  # geom_text(data = trait_vect,
  #           aes(x = PCo1, y = PCo1, label = trait),
  #           hjust = 0.5, vjust = .5, color="black", size = 3.5)
  geom_label(
    data = trait_vect,
    aes(x = PCo1, y = PCo2, label = trait),
    hjust = ifelse(names(trait_vect)[1]=='PCo1',0.2, .56),
    vjust = 1,
    size = 3,
    fill = "white",      # background color
    color = "red",       # text color
    label.size = 0.2,    # border thickness; set to 0 to remove outline
    inherit.aes = FALSE
  )+theme(legend.position = 'none')+
  annotate("text",
           x = -Inf, y = Inf, # Position at top-left corner
           label = 'Probabilistic',
           hjust = 0, vjust = 1, # Justify text relative to corner
           size = 3, color = "black")
#redefine C with no legend

#define the plot, but don't plot it
cmbin_plot=gridExtra::grid.arrange(B,A,ncol=1)
#plot the final graph with shared legend
gridExtra::grid.arrange(cmbin_plot,tax_legend,heights=c(10,1))
savp(10,8, 'C://Users//andrew.caudillo.BUGLAB-I9//Box//NAMC//Research Projects//AIM//IncreaserDecreaser_OE//Updated_w_modelObj//prob_responses_traitspace_Ecoregions_260219_axes_1_2.png')

#general trend of traits across land uses at watershed scale
#plotting distinct because of the duplicated records that come with the ecoregion,
#land use, etc.



#average of traits weighted by presebnce
#regional_and_traits=plyr::join(All_Results,traits_only[,names(traits_only) %in% c(top5,'taxon')])
#write.csv(regional_and_traits,'C://Users//andrew.caudillo.BUGLAB-I9//Box//NAMC//Research Projects//AIM//IncreaserDecreaser_OE//regional_responses_and_traits.csv')
#Os=Os[,colnames(Os) %in% taxa_notraits_rare$taxon==F]
O_long=as.data.frame(t(Os))
O_long$taxon=names(Os)
#ProbOs=ProbOs[,colnames(ProbOs) %in% taxa_notraits_rare$taxon==F]
ProbOs=ProbOs[row.names(ProbOs) %in% failed_sites$sampleId==F,]
ProbOs_long=as.data.frame(t(ProbOs))
ProbOs_long$taxon=names(ProbOs)
ProbOs_long$status='Probabilistic'
O_long$status='Reference'
traits_only3$taxon=row.names(traits_only2)
O_long=plyr::join(O_long,traits_only3)#[,names(traits_only) %in% c(names(top5),'taxon')])
ProbOs_long=plyr::join(ProbOs_long,traits_only3)#[,names(traits_only) %in% c('taxon',names(top5))])
site_cols=colnames(O_long)[1:656]#(ncol(O_long)-8)]
traitcols=colnames(O_long)[659:ncol(O_long)]#(ncol(O_long)-4):ncol(O_long)]
Psite_cols=colnames(ProbOs_long)[1:349]#1:(ncol(ProbOs_long)-8)]
Ptraitcols=colnames(ProbOs_long)[352:ncol(ProbOs_long)]#(ncol(ProbOs_long)-4):ncol(ProbOs_long)]
O_long=as.data.frame(O_long %>%
                       tidyr::pivot_longer(
                         cols = tidyselect::all_of(site_cols),
                         names_to = "Site",
                         values_to = "Presence"
                       ) %>%
                       dplyr::summarise(
                         dplyr::across(
                           .cols = tidyselect::all_of(traitcols),
                           .fns  = ~ sum(.x * Presence, na.rm = TRUE) / sum(Presence, na.rm = TRUE),
                           .names = "{.col}_mean"
                         ),
                         n_present = sum(Presence, na.rm = TRUE)
                       ))
P_long=as.data.frame(ProbOs_long %>%
                       tidyr::pivot_longer(
                         cols = all_of(Psite_cols),
                         names_to = "Site",
                         values_to = "Presence"
                       ) %>%
                       dplyr::summarise(
                         dplyr::across(
                           .cols = all_of(Ptraitcols),
                           .fns  = ~ sum(.x * Presence, na.rm = TRUE) / sum(Presence, na.rm = TRUE),
                           .names = "{.col}_mean"
                         ),
                         n_present = sum(Presence, na.rm = TRUE)
                       ))


compute_site_trait_means <- function(data, site_cols, trait_cols) {
  data %>%
    tidyr::pivot_longer(
      cols = tidyselect::all_of(site_cols),
      names_to = "Site",
      values_to = "Presence"
    ) %>%
    dplyr::group_by(Site) %>%
    dplyr::summarise(
      dplyr::across(
        tidyselect::all_of(trait_cols),
        ~ sum(.x * Presence, na.rm = TRUE) /
          sum(Presence, na.rm = TRUE)
      ),
      .groups = "drop"
    )
}

compute_design_means <- function(data, site_cols, trait_cols) {
  compute_site_trait_means(data, site_cols, trait_cols) %>%
    dplyr::summarise(
      dplyr::across(
        everything(),
        ~ mean(.x, na.rm = TRUE)
      )
    )
}
omeans=compute_design_means(O_long,site_cols,trait_cols)
pmeans=compute_design_means(ProbOs_long,Psite_cols,Ptraitcols)
omeans
pmeans
pmeans=pmeans[,names(omeans)]
omeans$status='Reference'
pmeans$status='Probabilistic'


combined_means <- dplyr::bind_rows(
  omeans %>% dplyr::mutate(status = "Reference"),
  pmeans %>% dplyr::mutate(status = "Probabilistic")
) %>%
  dplyr::summarise(
    dplyr::across(
      -c(status),
      ~ mean(.x, na.rm = TRUE)
    )
  )

O_long$status='Reference'
P_long$status='Probabilistic'
wm_traits=rbind(omeans,pmeans)
wm_traits[,names(wm_traits) %in% names(top5)]
clipr::write_clip(wm_traits)



top5

criterion<-function(values) {
  values ==0
}
common_names=intersect(names(PFos_oth),names(Oth_Fo))
crit1<-PFos_oth >0
crit2<-Oth_Fo >0
matching_crit=common_names[crit1 & crit2]
matching_crit
common_names=intersect(names(X_Fo),names(Px_Fo))
crit1<-X_Fo>0
crit2<-Px_Fo >0
matching_crit=common_names[crit1 & crit2]
X_Fo[matching_crit]

#validation data
val_bugs=val_bugs[,names(val_Pcs)]
val_bugs=val_bugs[match(row.names(val_Pcs), row.names(val_bugs)),]
row.names(val_bugs)==row.names(val_Pcs)
names(val_bugs)[names(val_bugs)=='ABEDUS']='Abedus'

val_bugs=val_bugs[,names(val_bugs) %in% c('Carabidae','Curculionidae')==F]
val_Pcs=val_Pcs[,names(val_Pcs) %in% c('Carabidae','Curculionidae')==F]
names(val_Pcs)[names(val_Pcs)=='ABEDUS']='Abedus'

valFos=colSums(val_bugs)
valFes=colSums(val_Pcs)
names(valFos)==names(valFes)
val_ratio=valFos/valFes
val_diffs=c('Soyedina','Nemertea','Margaritifera','Abedus','Pedicia','Mystacides')
val_plotdat=data.frame(Fo=valFos,Fe=valFes,taxon=names(valFos))
val_plotdat$col=ifelse(row.names(val_plotdat) %in% val_diffs,
                        'red',NA)
val_highlight=val_plotdat[which(val_plotdat$col=='red'),]
xeric_V=val_bugs[row.names(val_bugs) %in% Xeric_P,]
xeric_V_pcs=val_Pcs[row.names(val_Pcs) %in% Xeric_P,]
xeric_V_Fos=colSums(xeric_V)
xeric_V_Pcs=colSums(xeric_V_pcs)
xeric_val_plotdat=data.frame(Fo=xeric_V_Fos,Fe=xeric_V_Pcs,taxon=names(valFos))
xeric_val_plotdat$col=ifelse(row.names(xeric_val_plotdat) %in% val_diffs,
                             'red',NA)
xeric_val_highlight=xeric_val_plotdat[which(xeric_val_plotdat$col=='red'),]
Oth_V=val_bugs[!(row.names(val_bugs) %in% Xeric_P),]
Oth_V_pcs=val_Pcs[!(row.names(val_Pcs) %in% Xeric_P),]
oth_v_Fos=colSums(Oth_V)
oth_v_Fes=colSums(Oth_V_pcs)
oth_val_plotdat=data.frame(Fo=oth_v_Fos,Fe=oth_v_Fes,taxon=names(valFos))
oth_val_plotdat$col=ifelse(row.names(oth_val_plotdat) %in% val_diffs,
                           'red',NA)
oth_val_highlight=oth_val_plotdat[which(oth_val_plotdat$col=='red'),]
R=ggplot(data=Ref_Results,aes(y=Fo,x=Fe))+geom_point()+
  geom_abline(intercept = 0,slope = 1,col='red')+
  annotate("text",
           x = -Inf, y = Inf, # Position at top-left corner
           label = 'Reference',
           hjust = 0, vjust = 1, # Justify text relative to corner
           size = 5, color = "black")+
  lims(x=c(0,10),y=c(0,10))

V=ggplot(data=val_plotdat,aes(y=Fo,x=Fe))+geom_point()+
  geom_point(data=val_highlight,aes(x=Fe,y=Fo,color=taxon))+
  geom_abline(intercept = 0,slope = 1,col='red')+
  annotate("text",
           x = -Inf, y = Inf, # Position at top-left corner
           label = 'Validation',
           hjust = 0, vjust = 1, # Justify text relative to corner
           size = 5, color = "black")+
  scale_color_manual(values=c('Soyedina' = 'red',
                              'Nemertea' = 'dodgerblue',
                              'Margaritifera' = 'purple',
                              'Abedus' = 'orange',
                              'Pedicia' = 'blue',
                              'Mystacides' = 'yellow4'))+
  theme(legend.position = 'bottom')+
  lims(x=c(0,10),y=c(0,10))
gridExtra::grid.arrange(R,V,ncol=1)
savp(10,8,'C://Users//andrew.caudillo.BUGLAB-I9//Box//NAMC//Research Projects//AIM//IncreaserDecreaser_OE//Updated_w_modelObj//Ref_vs_val_FoFe_n10.png' )



#validation and ref plots

A<-ggplot(data=OthEco_plotdat,aes(y=Fo,x=Fe))+geom_point()+geom_abline(intercept = 0,slope = 1,col='red')+ggtitle('Other Ecoregions')+
  xlim(0,50)+ylim(0,50)+
  annotate("text",
           x = -Inf, y = Inf, # Position at top-left corner
           label = 'Reference',
           hjust = 0, vjust = 1, # Justify text relative to corner
           size = 4, color = "black")
B<-ggplot(data=EastXer_plotdat,aes(y=Fo,x=Fe))+geom_point()+
  ylim(0,50)+xlim(0,50)+
  geom_abline(intercept = 0,slope = 1,col='red')+ggtitle('Eastern Xeric')+
  annotate("text",
           x = -Inf, y = Inf, # Position at top-left corner
           label = 'Reference',
           hjust = 0, vjust = 1, # Justify text relative to corner
           size = 4, color = "black")



C<-ggplot(data=oth_val_plotdat,aes(y=Fo,x=Fe))+geom_point()+geom_abline(intercept = 0,slope = 1,col='red')+
  geom_point(data=oth_val_highlight,aes(x=Fe,y=Fo,color=taxon))+
  scale_color_manual(values=c('Soyedina' = 'red',
                              'Nemertea' = 'dodgerblue',
                              'Margaritifera' = 'purple',
                              'Abedus' = 'orange',
                              'Pedicia' = 'blue',
                              'Mystacides' = 'yellow4'))+
  annotate("text",
           x = -Inf, y = Inf, # Position at top-left corner
           label = 'Validation',
           hjust = 0, vjust = 1, # Justify text relative to corner
           size = 4, color = "black")+
  theme(legend.position = "bottom")
D<-ggplot(data=xeric_val_plotdat,aes(y=Fo,x=Fe))+geom_point()+geom_abline(intercept = 0,slope = 1,col='red')+
  geom_point(data=xeric_val_highlight,aes(x=Fe,y=Fo,color=taxon))+
  scale_color_manual(values=c('Soyedina' = 'red',
                              'Nemertea' = 'dodgerblue',
                              'Margaritifera' = 'purple',
                              'Abedus' = 'orange',
                              'Pedicia' = 'blue',
                              'Mystacides' = 'yellow4'))+
  annotate("text",
           x = -Inf, y = Inf, # Position at top-left corner
           label = 'Validation',
           hjust = 0, vjust = 1, # Justify text relative to corner
           size = 4, color = "black")+
  theme(legend.position = "none")+
  lims(x=c(0,50),y=c(0,50))
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
C<-ggplot(data=oth_val_plotdat,aes(y=Fo,x=Fe))+geom_point()+geom_abline(intercept = 0,slope = 1,col='red')+
  geom_point(data=oth_val_highlight,aes(x=Fe,y=Fo,color=taxon))+
  scale_color_manual(values=c('Soyedina' = 'red',
                              'Nemertea' = 'dodgerblue',
                              'Margaritifera' = 'purple',
                              'Abedus' = 'orange',
                              'Pedicia' = 'blue',
                              'Mystacides' = 'yellow4'))+
  annotate("text",
           x = -Inf, y = Inf, # Position at top-left corner
           label = 'Validation',
           hjust = 0, vjust = 1, # Justify text relative to corner
           size = 4, color = "black")+
  theme(legend.position = "none")
#define the plot, but don't plot it
cmbin_plot=gridExtra::grid.arrange(A,B,C,D,ncol=2)
#plot the final graph with shared legend
gridExtra::grid.arrange(cmbin_plot,tax_legend,heights=c(10,1))
savp(10,8, 'C://Users//andrew.caudillo.BUGLAB-I9//Box//NAMC//Research Projects//AIM//IncreaserDecreaser_OE//Updated_w_modelObj//Fo_Fe_RandV_ecos_same_axes.png')
val_OE=OE_calc(val_Pcs,
        val_bugs,
        threshold=0)
mean(val_OE$OtoE)
sd(val_OE$OtoE)
#val and ref O/Es
val_xer_OE=OE_calc(xeric_V_pcs,
                   xeric_V,
                   threshold = 0)
val_oth_OE=OE_calc(Oth_V_pcs,
                   Oth_V,
                   threshold = 0)
boxplot(ref_OEs_other$OtoE,at=2,xlim=c(1,6),col='purple4')
boxplot(val_oth_OE$OtoE,at=3,add=T,col='green4')
boxplot(XerOE$OtoE,at=4,add=T,col='purple4')
boxplot(val_xer_OE$OtoE,at=5,add=T,col='green4')
mtext(text = c("Other", "Xeric"), side = 1, line = 1, at = c(2.5, 4.5), cex = 1)
mtext(text=c('Ref.','Val.','Ref.','Val.'),side=1,line=2, at=c(2,3,4,5),cex=0.7)
legend('topright',
       leg=c('Reference','Validation'),
       pch=rep(22,2),
       pt.bg=c('purple4','green4'),
       bty='n',
       cex=0.8)
savp(10,8,'C://Users//andrew.caudillo.BUGLAB-I9//Box//NAMC//Research Projects//AIM//IncreaserDecreaser_OE//Updated_w_modelObj//Ref_vs_Val_OEs_0.png')
