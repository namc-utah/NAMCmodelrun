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
library(data.table)

library(dplyr)
library(tidyverse)
library(tidyr)
savp = function(W,H,fn) {
  dev.copy(dev=png,file=fn,wi=W,he=H,un="in",res=650)
  dev.off()
}

taxa_notraits_rare=read.csv('C://Users//andrew.caudillo.BUGLAB-I9//Box//NAMC//Research Projects//AIM//IncreaserDecreaser_OE//taxa_to_drop.csv')
nonrare=read.csv('C://Users//andrew.caudillo.BUGLAB-I9//Box//NAMC//Research Projects//AIM//IncreaserDecreaser_OE//nonrare_taxa.csv')

Os<-read.csv("C://Users//andrew.caudillo.BUGLAB-I9//Box//NAMC//Research Projects//AIM//IncreaserDecreaser_OE//MRF_Ref_Os.csv")
#Os=read.csv('C://Users//andrew.caudillo.BUGLAB-I9//Box//NAMC//Research Projects//AIM//IncreaserDecreaser_OE//Observed_PA.csv')
#reference predictions
 # Os<-Os%>%
 #   mutate(Tanypodinae = pmax(Tanypodinae, Tanypodinae.1, na.rm = TRUE)) %>%
 #   select(-Tanypodinae.1)
 #Es=read.csv('C://Users//andrew.caudillo.BUGLAB-I9//Box//NAMC//Research Projects//AIM//IncreaserDecreaser_OE//Ref_Expected.csv')
 Es<-read.csv("C://Users//andrew.caudillo.BUGLAB-I9//Box//NAMC//Research Projects//AIM//IncreaserDecreaser_OE//MRF_Ref_Pcs.csv")
#  Es<-Es%>%
#    mutate(Tanypodinae = pmax(Tanypodinae, Tanypodinae.1, na.rm = TRUE)) %>%
#    select(-Tanypodinae.1)
# # #dropping 3 problem sites for which predictors could not be calculated
 Es<-Es[Es$X %in% c(116830,132243,132807, 185022)==F,]
 Os<-Os[Os$X %in% c(116830,132243,132807,185022)==F,]

# #ensure that the column order matches across the tables
# #i.e., the taxa names
 setcolorder(Os,neworder = names(Es))
#probabilistic predictions
Pes<-read.csv("C://Users//andrew.caudillo.BUGLAB-I9//Box//NAMC//Research Projects//AIM//IncreaserDecreaser_OE//MRF_Prob_Pcs.csv")
#Pes$Tanypodinae = Pes$Tanypodinae + Pes$Tanypodinae.1;Pes<-Pes[,names(Pes) %in% 'Tanypodinae.1'==F,]
failed_sites<-read.csv('C://Users//andrew.caudillo.BUGLAB-I9//Box//NAMC//Research Projects//AIM//IncreaserDecreaser_OE//failed_sites.csv')
failed<-Pes[Pes$X %in% failed_sites$sampleId,]
Pes<-Pes[Pes$X %in% failed_sites$sampleId==F,]
ProbOs<-read.csv('C://Users//andrew.caudillo.BUGLAB-I9//Box//NAMC//Research Projects//AIM//IncreaserDecreaser_OE//MRF_Prob_Os.csv')
ProbOs=ProbOs[ProbOs$X %in% Pes$X,]
#ProbOs$Tanypodinae = ProbOs$Tanypodinae + ProbOs$Tanypodinae.1;ProbOs<-ProbOs[,names(ProbOs) %in% 'Tanypodinae.1'==F,]
# setcolorder(ProbOs,neworder = names(Es))
# names(Os)==names(Pes)
# names(Os)==names(Es)
# names(Os)==names(failed)
# names(Os)==names(ProbOs)
row.names(failed)=failed$sampleId;failed<-failed[,-1]
#dropping sampleId and assigning that to the row names
 row.names(Os)<-Os$X;Os<-Os[,-1]
 row.names(Es)<-Es$X;Es<-Es[,-1]
 row.names(Pes)<-Pes$X;Pes<-Pes[,-1]
 row.names(ProbOs)<-ProbOs$X;ProbOs<-ProbOs[,-1]
N_Tax=ncol(Es)
refE<-colSums(Es)
ProbE<-colSums(Pes)
refO<-colSums(Os)
ProbO<-colSums(ProbOs)
refFoFe=refO / refE
probFoFe= ProbO / ProbE
FailedE<-colSums(failed)


Ref_Results=data.frame(taxon=names(Os),
                                Fo=refO,
                                Fe=refE,
                                ratio = refFoFe)



P_Results=data.frame(taxon=names(Os),
                     Fo=ProbO,
                     Fe=ProbE,
                     ratio = probFoFe)

#P_Results<-P_Results[P_Results$taxon %in% taxa_notraits_rare==F,]
#Ref_Results<-Ref_Results[Ref_Results$taxon %in% taxa_notraits_rare==F,]
#Ref_Results<-Ref_Results[Ref_Results$taxon %in% P_Results$taxon,]
#Ref_Results$taxon
#a few reference taxa are increasers, but this could be due to the cutoff.
#The 4 taxa (Cinygma, Hesperoconopa, Limnophila, and Paraperla) have a ratio
# > 1.2 but < 1.3.
Ref_Results$Regional_Response=ifelse(Ref_Results$ratio >1.2, "Increaser",
                                 ifelse( Ref_Results$ratio > 0.8 & Ref_Results$ratio < 1.2, 'Neutral',"Decreaser"))
Ref_Results$status='Reference'
P_Results$Regional_Response=ifelse(P_Results$ratio >1.2, "Increaser",
                                     ifelse( P_Results$ratio > 0.8 & P_Results$ratio < 1.2, 'Neutral',"Decreaser"))
P_Results$status='Probabilistic'
All_Results<-rbind(Ref_Results,P_Results)
All_Results<-All_Results[All_Results$taxon %in% nonrare$taxon==T,]
All_Results<-All_Results[All_Results$taxon %in% taxa_notraits_rare$taxon==F,]
write.csv(All_Results,'C://Users//andrew.caudillo.BUGLAB-I9//Box//NAMC//Research Projects//AIM//IncreaserDecreaser_OE//regional_responses_MRF_251107.csv')
# Dropped_Taxa_results<-Results[Results$O < 10 | Results$ProbE <10,]
#
# Dropped_Taxa_results$Regional_Ratio = Dropped_Taxa_results$ProbE / Dropped_Taxa_results$RefE
# Dropped_Taxa_results$FRatio = Dropped_Taxa_results$ProbE / Dropped_Taxa_results$failedE
# Results$Regional_Ratio=Results$ProbE/Results$RefE
# Results$failed_ratio = Results$failedE / Results$RefE
#
# Results$FResponse=ifelse(Results$failed_ratio >1.2, "Increaser",
#                         ifelse( Results$failed_ratio > 0.8 & Results$failed_ratio < 1.2, 'Neutral',"Decreaser"))
#
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
groups$EPA_hybird==ifelse(groups$EPA_hybird=='Eastern Xeric Plateaus',
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



# classify_site <- function(pct_natural, pct_agriculture, pct_urban) {
#   # Put values in a named vector
#   p <- c(Shrubland = pct_natural,
#          Agriculture = pct_agriculture,
#          Urban = pct_urban)
#
#   # Pre-calculate
#   top_class <- names(p)[which.max(p)]
#   top_val <- max(p)
#   over_33 <- names(p[p > 33.33])
#
#   # --- TIER 1: Strong dominance (>66.66%) ---
#   if (top_val > 66.66) {
#     return(paste0(top_class))
#   }
#
#   # --- TIER 2: Moderate dominance (50-66.66%) ---
#   if (top_val > 50 && top_val <= 66.66) {
#     return(paste0(top_class))
#   }
#
#   # --- TIER 4: Mixed (two or more >33%) ---
#   if (length(over_33) >= 2) {
#     return("Mixed")
#   }
#
#   # --- TIER 3: Weak dominance (33-50%), only if NOT mixed ---
#   if (top_val > 33.33 && top_val <= 50) {
#     return(paste0(top_class))
#   }
#
#   # --- Otherwise: Low dominance (all <=33%) ---
#   return(paste0(top_class))
# }

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
#sites<-sites[sites$sampleId %in% failed_sites==F,]
#row.names(O)
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
#Pesgrp1<-Pes[row.names(Pes) %in% groups[[1]]$sampleId,]
#Pesgrp1$Ecoregion=rep(unique(groups[[1]]$EPA_hybird),nrow(Pesgrp1))

#this section calcualtes the Pc for each probabilistic and
#reference site for each ecoregion.
#It also assigns land use to sites within each ecoregion
}
if(1){
#add another section for Prob Os
Pesgrp2<-Pes[row.names(Pes) %in% groups$sampleId[groups$EPA_hybird=='Eastern Xeric Plateaus'],]
Pesgrp2$Ecoregion='Eastern Xeric Plateaus'
Pesgrp3<-Pes[row.names(Pes) %in% groups$sampleId[groups$EPA_hybird!='Eastern Xeric Plateaus'],]
Pesgrp3$Ecoregion='Other'
# Pesgrp4<-Pes[row.names(Pes) %in% groups[[4]]$sampleId,]
# Pesgrp4$Ecoregion=rep(unique(groups[[4]]$EPA_hybird),nrow(Pesgrp4))
# Pesgrp5<-Pes[row.names(Pes) %in% groups[[5]]$sampleId,]
# Pesgrp5$Ecoregion=rep(unique(groups[[5]]$EPA_hybird),nrow(Pesgrp5))
# Pesgrp6<-Pes[row.names(Pes) %in% groups[[6]]$sampleId,]
# Pesgrp6$Ecoregion=rep(unique(groups[[6]]$EPA_hybird),nrow(Pesgrp6))
# Pesgrp7<-Pes[row.names(Pes) %in% groups[[7]]$sampleId,]
# Pesgrp7$Ecoregion=rep(unique(groups[[7]]$EPA_hybird),nrow(Pesgrp7))
# Pesgrp8<-Pes[row.names(Pes) %in% groups[[8]]$sampleId,]
# Pesgrp8$Ecoregion=rep(unique(groups[[8]]$EPA_hybird),nrow(Pesgrp8))
# Pesgrp9<-Pes[row.names(Pes) %in% groups[[9]]$sampleId,]
# Pesgrp9$Ecoregion=rep(unique(groups[[9]]$EPA_hybird),nrow(Pesgrp9))
# Pesgrp10<-Pes[row.names(Pes) %in% groups[[10]]$sampleId,]
# Pesgrp10$Ecoregion=rep(unique(groups[[10]]$EPA_hybird),nrow(Pesgrp10))

ProbOsgrp2<-ProbOs[row.names(ProbOs) %in% groups$sampleId[groups$EPA_hybird=='Eastern Xeric Plateaus'],]
ProbOsgrp2$Ecoregion='Eastern Xeric Plateaus'
ProbOsgrp3<-ProbOs[row.names(ProbOs) %in% groups$sampleId[groups$EPA_hybird!='Eastern Xeric Plateaus'],]
ProbOsgrp3$Ecoregion='Other'
# ProbOsgrp4<-ProbOs[row.names(ProbOs) %in% groups[[4]]$sampleId,]
# ProbOsgrp4$Ecoregion=rep(unique(groups[[4]]$EPA_hybird),nrow(ProbOsgrp4))
# ProbOsgrp5<-ProbOs[row.names(ProbOs) %in% groups[[5]]$sampleId,]
# ProbOsgrp5$Ecoregion=rep(unique(groups[[5]]$EPA_hybird),nrow(ProbOsgrp5))
# ProbOsgrp6<-ProbOs[row.names(ProbOs) %in% groups[[6]]$sampleId,]
# ProbOsgrp6$Ecoregion=rep(unique(groups[[6]]$EPA_hybird),nrow(ProbOsgrp6))
# ProbOsgrp7<-ProbOs[row.names(ProbOs) %in% groups[[7]]$sampleId,]
# ProbOsgrp7$Ecoregion=rep(unique(groups[[7]]$EPA_hybird),nrow(ProbOsgrp7))
# ProbOsgrp8<-ProbOs[row.names(ProbOs) %in% groups[[8]]$sampleId,]
# ProbOsgrp8$Ecoregion=rep(unique(groups[[8]]$EPA_hybird),nrow(ProbOsgrp8))
# ProbOsgrp9<-ProbOs[row.names(ProbOs) %in% groups[[9]]$sampleId,]
# ProbOsgrp9$Ecoregion=rep(unique(groups[[9]]$EPA_hybird),nrow(ProbOsgrp9))
# ProbOsgrp10<-ProbOs[row.names(ProbOs) %in% groups[[10]]$sampleId,]
# ProbOsgrp10$Ecoregion=rep(unique(groups[[10]]$EPA_hybird),nrow(ProbOsgrp10))


Osgrp2<-Os[row.names(Os) %in% groups$sampleId[groups$EPA_hybird=='Eastern Xeric Plateaus'],]
Osgrp2$Ecoregion='Eastern Xeric Plateaus'
Osgrp3<-Os[row.names(Os) %in% groups$sampleId[groups$EPA_hybird!='Eastern Xeric Plateaus'],]
Osgrp3$Ecoregion='Other'
# Osgrp4<-Os[row.names(Os) %in% groups[[4]]$sampleId,]
# Osgrp4$Ecoregion=rep(unique(groups[[4]]$EPA_hybird),nrow(Osgrp4))
# Osgrp5<-Os[row.names(Os) %in% groups[[5]]$sampleId,]
# Osgrp5$Ecoregion=rep(unique(groups[[5]]$EPA_hybird),nrow(Osgrp5))
# Osgrp6<-Os[row.names(Os) %in% groups[[6]]$sampleId,]
# Osgrp6$Ecoregion=rep(unique(groups[[6]]$EPA_hybird),nrow(Osgrp6))
# Osgrp7<-Os[row.names(Os) %in% groups[[7]]$sampleId,]
# Osgrp7$Ecoregion=rep(unique(groups[[7]]$EPA_hybird),nrow(Osgrp7))
# Osgrp8<-Os[row.names(Os) %in% groups[[8]]$sampleId,]
# Osgrp8$Ecoregion=rep(unique(groups[[8]]$EPA_hybird),nrow(Osgrp8))
# Osgrp9<-Os[row.names(Os) %in% groups[[9]]$sampleId,]
# Osgrp9$Ecoregion=rep(unique(groups[[9]]$EPA_hybird),nrow(Osgrp9))
# Osgrp10<-Os[row.names(Os) %in% groups[[10]]$sampleId,]
# Osgrp10$Ecoregion=rep(unique(groups[[10]]$EPA_hybird),nrow(Osgrp10))
# #

Esgrp2<-Es[row.names(Es) %in% groups$sampleId[groups$EPA_hybird=='Eastern Xeric Plateaus'],]

Esgrp3<-Es[row.names(Es) %in% groups$sampleId[groups$EPA_hybird!='Eastern Xeric Plateaus'],]

# Esgrp4<-Es[row.names(Es) %in% groups[[4]]$sampleId,]
#
# Esgrp5<-Es[row.names(Es) %in% groups[[5]]$sampleId,]
#
# Esgrp6<-Es[row.names(Es) %in% groups[[6]]$sampleId,]
#
# Esgrp7<-Es[row.names(Es) %in% groups[[7]]$sampleId,]
#
# Esgrp8<-Es[row.names(Es) %in% groups[[8]]$sampleId,]
#
# Esgrp9<-Es[row.names(Es) %in% groups[[9]]$sampleId,]
#
# Esgrp10<-Es[row.names(Es) %in% groups[[10]]$sampleId,]
#

# L_grp2<-list(Pesgrp2,Esgrp2)
# L_grp3<-list(Pesgrp3,Esgrp3)
# L_grp4<-list(Pesgrp4,Esgrp4)
# L_grp5<-list(Pesgrp5,Esgrp5)
# L_grp6<-list(Pesgrp6,Esgrp6)
# L_grp7<-list(Pesgrp7,Esgrp7)
# L_grp8<-list(Pesgrp8,Esgrp8)
# L_grp9<-list(Pesgrp9,Esgrp9)
# L_grp10<-list(Pesgrp10,Esgrp10)

L_grp2<-list(ProbOsgrp2,Pesgrp2)
L_grp3<-list(ProbOsgrp3,Pesgrp3)
# L_grp4<-list(ProbOsgrp4,Pesgrp4)
# L_grp5<-list(ProbOsgrp5,Pesgrp5)
# L_grp6<-list(ProbOsgrp6,Pesgrp6)
# L_grp7<-list(ProbOsgrp7,Pesgrp7)
# L_grp8<-list(ProbOsgrp8,Pesgrp8)
# L_grp9<-list(ProbOsgrp9,Pesgrp9)
# L_grp10<-list(ProbOsgrp10,Pesgrp10)



Grp2_ratio=colSums(L_grp2[[1]][1:ncol(L_grp2[[1]])-1])/colSums((L_grp2[[2]][1:ncol(L_grp2[[2]])-1]))
Grp3_ratio=colSums(L_grp3[[1]][1:ncol(L_grp3[[1]])-1])/colSums((L_grp2[[2]][1:ncol(L_grp2[[2]])-1]))
# Grp4_ratio=colSums(L_grp4[[1]][1:ncol(L_grp4[[1]])-1])/colSums((L_grp2[[2]][1:ncol(L_grp2[[2]])-1]))
# Grp5_ratio=colSums(L_grp5[[1]][1:ncol(L_grp5[[1]])-1])/colSums((L_grp2[[2]][1:ncol(L_grp2[[2]])-1]))
# Grp6_ratio=colSums(L_grp6[[1]][1:ncol(L_grp6[[1]])-1])/colSums((L_grp2[[2]][1:ncol(L_grp2[[2]])-1]))
# Grp7_ratio=colSums(L_grp7[[1]][1:ncol(L_grp7[[1]])-1])/colSums((L_grp2[[2]][1:ncol(L_grp2[[2]])-1]))
# Grp8_ratio=colSums(L_grp8[[1]][1:ncol(L_grp8[[1]])-1])/colSums((L_grp2[[2]][1:ncol(L_grp2[[2]])-1]))
# Grp9_ratio=colSums(L_grp9[[1]][1:ncol(L_grp9[[1]])-1])/colSums((L_grp2[[2]][1:ncol(L_grp2[[2]])-1]))
# Grp10_ratio=colSums(L_grp10[[1]][1:ncol(L_grp10[[1]])-1])/colSums((L_grp2[[2]][1:ncol(L_grp2[[2]])-1]))


Land_use_tidy<-Land_use
names(Land_use_tidy)[1]<-'ID'

P2_tidy <- Pesgrp2 %>%
  rownames_to_column("ID")
P3_tidy <- Pesgrp3 %>%
  rownames_to_column("ID")
# P4_tidy <- Pesgrp4 %>%
#   rownames_to_column("ID")
# P5_tidy <- Pesgrp5 %>%
#   rownames_to_column("ID")
# P6_tidy <- Pesgrp6 %>%
#   rownames_to_column("ID")
# P7_tidy <- Pesgrp7 %>%
#   rownames_to_column("ID")
# P8_tidy <- Pesgrp8 %>%
#   rownames_to_column("ID")
# P9_tidy <- Pesgrp9 %>%
#   rownames_to_column("ID")
# P10_tidy <- Pesgrp10 %>%
#   rownames_to_column("ID")

Po2_tidy <- ProbOsgrp2 %>%
  rownames_to_column("ID")
Po3_tidy <- ProbOsgrp3 %>%
  rownames_to_column("ID")
# Po4_tidy <- ProbOsgrp4 %>%
#   rownames_to_column("ID")
# Po5_tidy <- ProbOsgrp5 %>%
#   rownames_to_column("ID")
# Po6_tidy <- ProbOsgrp6 %>%
#   rownames_to_column("ID")
# Po7_tidy <- ProbOsgrp7 %>%
#   rownames_to_column("ID")
# Po8_tidy <- ProbOsgrp8 %>%
#   rownames_to_column("ID")
# Po9_tidy <- ProbOsgrp9 %>%
#   rownames_to_column("ID")
# Po10_tidy <- ProbOsgrp10 %>%
#   rownames_to_column("ID")




P2=plyr::join(P2_tidy,Land_use_tidy[,c('Land_use','Local_land_use','ID')],by='ID')
P3=plyr::join(P3_tidy,Land_use_tidy[,c('Land_use','Local_land_use','ID')],by='ID')
# P4=plyr::join(P4_tidy,Land_use_tidy[,c('Land_use','Local_land_use','ID')],by='ID')
# P5=plyr::join(P5_tidy,Land_use_tidy[,c('Land_use','Local_land_use','ID')],by='ID')
# P6=plyr::join(P6_tidy,Land_use_tidy[,c('Land_use','Local_land_use','ID')],by='ID')
# P7=plyr::join(P7_tidy,Land_use_tidy[,c('Land_use','Local_land_use','ID')],by='ID')
# P8=plyr::join(P8_tidy,Land_use_tidy[,c('Land_use','Local_land_use','ID')],by='ID')
# P9=plyr::join(P9_tidy,Land_use_tidy[,c('Land_use','Local_land_use','ID')],by='ID')
# P10=plyr::join(P10_tidy,Land_use_tidy[,c('Land_use','Local_land_use','ID')],by='ID')

Po2=plyr::join(Po2_tidy,Land_use_tidy[,c('Land_use','Local_land_use','ID')],by='ID')
Po3=plyr::join(Po3_tidy,Land_use_tidy[,c('Land_use','Local_land_use','ID')],by='ID')
# Po4=plyr::join(Po4_tidy,Land_use_tidy[,c('Land_use','Local_land_use','ID')],by='ID')
# Po5=plyr::join(Po5_tidy,Land_use_tidy[,c('Land_use','Local_land_use','ID')],by='ID')
# Po6=plyr::join(Po6_tidy,Land_use_tidy[,c('Land_use','Local_land_use','ID')],by='ID')
# Po7=plyr::join(Po7_tidy,Land_use_tidy[,c('Land_use','Local_land_use','ID')],by='ID')
# Po8=plyr::join(Po8_tidy,Land_use_tidy[,c('Land_use','Local_land_use','ID')],by='ID')
# Po9=plyr::join(Po9_tidy,Land_use_tidy[,c('Land_use','Local_land_use','ID')],by='ID')
# Po10=plyr::join(Po10_tidy,Land_use_tidy[,c('Land_use','Local_land_use','ID')],by='ID')



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
# Esgrp4$Ecoregion=Osgrp4$Ecoregion
# Esgrp5$Ecoregion=Osgrp5$Ecoregion
# Esgrp6$Ecoregion=Osgrp6$Ecoregion
# Esgrp7$Ecoregion=Osgrp7$Ecoregion
# Esgrp8$Ecoregion=Osgrp8$Ecoregion
# Esgrp9$Ecoregion=Osgrp9$Ecoregion
# Esgrp10$Ecoregion=Osgrp10$Ecoregion


E2_tidy <- Esgrp2 %>%
  rownames_to_column("ID")
E3_tidy <- Esgrp3 %>%
  rownames_to_column("ID")
# E4_tidy <- Esgrp4 %>%
#   rownames_to_column("ID")
# E5_tidy <- Esgrp5 %>%
#   rownames_to_column("ID")
# E6_tidy <- Esgrp6 %>%
#   rownames_to_column("ID")
# E7_tidy <- Esgrp7 %>%
#   rownames_to_column("ID")
# E8_tidy <- Esgrp8 %>%
#   rownames_to_column("ID")
# E9_tidy <- Esgrp9 %>%
#   rownames_to_column("ID")
# E10_tidy <- Esgrp10 %>%
#   rownames_to_column("ID")
# E2=plyr::join(E2_tidy,Land_use_tidy[,c('Land_use','Local_land_use','ID')],by='ID')
# E3=plyr::join(E3_tidy,Land_use_tidy[,c('Land_use','Local_land_use','ID')],by='ID')
# E4=plyr::join(E4_tidy,Land_use_tidy[,c('Land_use','Local_land_use','ID')],by='ID')
# E5=plyr::join(E5_tidy,Land_use_tidy[,c('Land_use','Local_land_use','ID')],by='ID')
# E6=plyr::join(E6_tidy,Land_use_tidy[,c('Land_use','Local_land_use','ID')],by='ID')
# E7=plyr::join(E7_tidy,Land_use_tidy[,c('Land_use','Local_land_use','ID')],by='ID')
# E8=plyr::join(E8_tidy,Land_use_tidy[,c('Land_use','Local_land_use','ID')],by='ID')
# E9=plyr::join(E9_tidy,Land_use_tidy[,c('Land_use','Local_land_use','ID')],by='ID')
# E10=plyr::join(E10_tidy,Land_use_tidy[,c('Land_use','Local_land_use','ID')],by='ID')

O2_tidy <- Osgrp2 %>%
  rownames_to_column("ID")
O3_tidy <- Osgrp3 %>%
  rownames_to_column("ID")
# O4_tidy <- Osgrp4 %>%
#   rownames_to_column("ID")
# O5_tidy <- Osgrp5 %>%
#   rownames_to_column("ID")
# O6_tidy <- Osgrp6 %>%
#   rownames_to_column("ID")
# O7_tidy <- Osgrp7 %>%
#   rownames_to_column("ID")
# O8_tidy <- Osgrp8 %>%
#   rownames_to_column("ID")
# O9_tidy <- Osgrp9 %>%
#   rownames_to_column("ID")
# O10_tidy <- Osgrp10 %>%
#   rownames_to_column("ID")

# O2=plyr::join(O2_tidy,Land_use_tidy[,c('Land_use','Local_land_use','ID')],by='ID')
# O3=plyr::join(O3_tidy,Land_use_tidy[,c('Land_use','Local_land_use','ID')],by='ID')
# O4=plyr::join(O4_tidy,Land_use_tidy[,c('Land_use','Local_land_use','ID')],by='ID')
# O5=plyr::join(O5_tidy,Land_use_tidy[,c('Land_use','Local_land_use','ID')],by='ID')
# O6=plyr::join(O6_tidy,Land_use_tidy[,c('Land_use','Local_land_use','ID')],by='ID')
# O7=plyr::join(O7_tidy,Land_use_tidy[,c('Land_use','Local_land_use','ID')],by='ID')
# O8=plyr::join(O8_tidy,Land_use_tidy[,c('Land_use','Local_land_use','ID')],by='ID')
# O9=plyr::join(O9_tidy,Land_use_tidy[,c('Land_use','Local_land_use','ID')],by='ID')
# O10=plyr::join(O10_tidy,Land_use_tidy[,c('Land_use','Local_land_use','ID')],by='ID')
All_Os<-do.call(rbind,list(O2_tidy,O3_tidy))#O4_tidy,O5_tidy,O6_tidy,O7_tidy,O8_tidy,O9_tidy,O10_tidy))
All_Os$Status='rObserved'
All_ref=do.call(rbind, list(E2_tidy,E3_tidy))#,E4_tidy,E5_tidy,E6_tidy,E7_tidy,E8_tidy,E9_tidy,E10_tidy))
All_ref$Status='rExpected'
#All_ref<-All_ref[,-ncol(All_ref)]
#All_ref$Status='Reference'
All_ref=rbind(All_Os,All_ref)

#All_data=
All_data=dplyr::bind_rows(All_ref,All_Prob_sites)

Os_lookup<-data.frame(ID=All_Os$ID,
                      O=NA)

otaxa_start <- 1
otaxa_end <- ncol(All_Os) - 2  # because last two columns are ecoregion & land_use

# O_tidy <- All_Os %>%
#   pivot_longer(
#     cols = otaxa_start:otaxa_end,
#     names_to = "taxon",
#     values_to = "occurrence"
#   ) %>%
#   group_by(Ecoregion, taxon) %>%
#   summarize(O_count = sum(occurrence, na.rm = TRUE), .groups="drop")
# }

row.names(All_Os)<-All_Os$ID;All_Os<-All_Os[,-1]
#rFos=colSums(All_Os)
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

Long_Prob<-All_Prob_O%>%
  pivot_longer(
    cols=taxa_start:taxa_end,
    names_to = 'taxon',
    values_to = 'value'
  )
df_long <- All_data %>%
  pivot_longer(
    cols=taxa_start:taxa_end,
    names_to='taxon',
    values_to='value'
  )


# O_sub=df_long[df_long$Status=='rObserved',]
# rFo=O_sub  %>%
#   group_by(taxon, Ecoregion) %>%
#   summarise(Fo = sum(value, na.rm = TRUE), .groups = "drop")
# rFe=df_long%>%
#   filter(Status=='rExpected')%>%
#   group_by(taxon, Ecoregion) %>%
#   summarise(Fe = sum(value, na.rm = TRUE), .groups = "drop")

pFo= aggregate(
  value ~ Ecoregion + Land_use + Local_land_use + taxon,
  data = df_long[df_long$Status == "PObserved", ],
  FUN = sum,
  na.rm = TRUE
)
pFe= aggregate(
  value ~ Ecoregion + Land_use + Local_land_use + taxon,
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
names(Ref_Ecoregions)[3:4]<-c('rFo','rFe')
Ref_Ecoregions$Rratio = Ref_Ecoregions$rFo / Ref_Ecoregions$rFe
Ref_Ecoregions$Rresponse=ifelse(Ref_Ecoregions$Rratio > 1.2,'Increaser',
                               ifelse(Ref_Ecoregions$Rratio > 0.8 & Ref_Ecoregions$Rratio < 1.2, 'Neutral',
                                      'Decreaser'))

P_Ecoregions=plyr::join(pFo,pFe,by=c('Ecoregion','Land_use','Local_land_use','taxon'))
names(P_Ecoregions)[5:6]<-c('pFo','pFe')
P_Ecoregions$Pratio = P_Ecoregions$pFo / P_Ecoregions$pFe
P_Ecoregions$Presponse=ifelse(P_Ecoregions$Pratio > 1.2,'Increaser',
                               ifelse(P_Ecoregions$Pratio > 0.8 & P_Ecoregions$Pratio < 1.2, 'Neutral',
                                      'Decreaser'))

All_ratios=Ref_Ecoregions %>% left_join(P_Ecoregions, by=c('Ecoregion','taxon'))
All_ratios<-All_ratios[All_ratios$taxon %in% All_Results$taxon,]
# df_long <- bind_rows(
#   All_data %>%
#     filter(Status == "PExpected") %>%   # Only Probabilistic rows
#     pivot_longer(
#       cols = taxa_start:taxa_end,           # or a pattern that matches only taxa columns
#       names_to = "taxon",
#       values_to = "Pc"
#     ),
#   All_data %>%
#     filter(Status == "Reference") %>%       # Only Reference rows
#     pivot_longer(
#       cols = taxa_start:taxa_end,
#       names_to = "taxon",
#       values_to = "Pc"
#     )
# )
# #
# # df_long <- All_data %>%
# #   #using the taxa columns,
# #   #assign the names to "taxon",
# #   #and the values to "Pc"
# #   pivot_longer(
# #     cols = taxa_start:taxa_end,   # or use all taxa column indices explicitly
# #     names_to = "taxon",
# #     values_to = "Pc"
# #   )
# #now that we have the long format,
# #now it is just a matter of grouping by
# #Ecoregion, land use, and taxon
# #then summing the Pcs.
# #start with reference sites
# ref_sums <- df_long %>%
#   filter(Status == "PObserved") %>%
#   group_by(Ecoregion, Land_use, Local_land_use, taxon) %>%
#   #there are no NAs in the dataset, but good to keep
#   #this argument, in case you have an NA.
#   #sums will be NA if present.
#   summarize(ref_p_sum = sum(value, na.rm = TRUE), .groups="drop")
# ref_sums<-ref_sums[,-c(2:3)]
# #same, but now for the probabilistic sites
# prob_sums <- df_long %>%
#   filter(Status == "Probabilistic") %>%
#   group_by(Ecoregion, Land_use, Local_land_use, taxon) %>%
#   summarize(prob_p_sum = sum(Pc, na.rm = TRUE), .groups="drop")
#
# #now, make the final table (minus traits)
# #join the two tables by hte group, type, and taxon columns
# final_table <- prob_sums %>%
#   left_join(ref_sums, by=c("Ecoregion","taxon")) %>%
#   #make a new column called ratio
#   #and then create the response
#   #based on the condition
#   #if the E' (or probabilistic E)
#   #is greater than the E,
#   #it will be an increaser.
#   #the magnitude depends on the size of E.
#   #i.e., 0.3 / 0.1 = 3.0
#   mutate(
#     ratio = ifelse(ref_p_sum > 0, prob_p_sum / ref_p_sum, NA_real_),
#     response = case_when(
#       ratio > 1.2 ~"Increaser",
#       ratio < 0.8  ~"Decreaser",
#       ratio >0.8 & ratio < 1.2 ~"Neutral",
#       is.na(ratio)~"Not expected"
#     )
#   )
#
#
#

#Trait info
#read in the traits from NAMC's database
#translationID=21
#contains bench ID and OTU
#this needs to only be run once, as the table does not change.
if(1){
trait_table<-read.csv('C://Users//andrew.caudillo.BUGLAB-I9//Box//NAMC//Research Projects//AIM//IncreaserDecreaser_OE//All_bench_taxa_atrributes2.csv')
OTU21<-read.csv('C://Users//andrew.caudillo.BUGLAB-I9//Box//NAMC//Research Projects//AIM//IncreaserDecreaser_OE//OTU21.csv')
trait_table<-trait_table[trait_table$OTU %in% OTU21$taxonOTU,]
trait_table<-trait_table[trait_table$OTU %in% All_ratios$taxon,]

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

ratios_w_traits=plyr::join(All_ratios,wide_trait_table_OTUs,by='taxon','left')
#creating filter column, removing duplicates, and then dropping that column.
ratios_w_traits<-ratios_w_traits[,names(ratios_w_traits) %in% 'scientific_name'==F]
#ratios_w_traits<-ratios_w_traits[!duplicated(ratios_w_traits$taxon),]
#dropping taxa that only occur at 10 or more sites
#should already be done, but just checking
#final_table3<-regional_tab_w_traits[regional_tab_w_traits$O>=10 | regional_tab_w_traits$ProbE>=10,]
table(ratios_w_traits$taxon)
write.csv(ratios_w_traits,'C://Users//andrew.caudillo.BUGLAB-I9//Box//NAMC//Research Projects//AIM//IncreaserDecreaser_OE//Spatial_scales_and_traits_MRF.csv',row.names = F)
}

#read in table that came from the trait section
final_table3<-read.csv('C://Users//andrew.caudillo.BUGLAB-I9//Box//NAMC//Research Projects//AIM//IncreaserDecreaser_OE//Spatial_scales_and_traits_MRF.csv')
#final_table3<-final_table3[final_table3$O > 10,]

#final_table3<-final_table3[final_table3$O > 10 |final_table3$ProbE > 10,]
#final_table3$taxon[final_table3$O / final_table3$RefE > 1.5]

a<-ggplot(All_Results[All_Results$status=='Reference',],aes(x=Fe,y=Fo))+geom_point()+geom_abline(color='red',linetype=2)+ggtitle('Reference')+labs(y='Fo',x='Fe')
b<-ggplot(P_Results[All_Results$status=='Probabilistic',],aes(x=Fe,y=Fo))+geom_point()+geom_abline(color='red',linetype=2)+ggtitle('Probabilistic')+labs(y='Fo',x='Fe')
gridExtra::grid.arrange(a,b,ncol=1)
savp(10,8, 'C://Users//andrew.caudillo.BUGLAB-I9//Box//NAMC//Research Projects//AIM//IncreaserDecreaser_OE//FovsFe_MRF.png')
boxplot(All_Results$ratio[All_Results$status=='Reference'],at=1,xlim=c(0,3),ylim=c(0,5))
boxplot(All_Results$ratio[All_Results$status=='Probabilistic'],at=2,add=T)
axis(side=1, labels = c('Ref.','Prob'),at=c(1,2))
savp(10,8, 'C://Users//andrew.caudillo.BUGLAB-I9//Box//NAMC//Research Projects//AIM//IncreaserDecreaser_OE//Ref_Prob_box_MRF.png')

par(mfrow=c(2,1))
hist(log(All_Results$Fo[All_Results$status=='Reference'] - All_Results$Fe[All_Results$status=='Reference']), xlab='',main='Ref.')

hist(log(All_Results$Fo[All_Results$status=='Probabilistic'] - All_Results$Fe[All_Results$status=='Probabilistic']),xlab='log(Fo-Fe)',main='Prob.')
savp(10,8, 'C://Users//andrew.caudillo.BUGLAB-I9//Box//NAMC//Research Projects//AIM//IncreaserDecreaser_OE//Ref_Prob_hist_MRF.png')
#these are taxa for which we have essentially 0 metrics for, so we cannot use them in the
#Dissimilarity matrix, for they will result in NAs in the final product.
#regex_patterns <- paste(c('idae','Nemata','inae','OLIGOCHAETA','PLATYHELMINTHES','OTHER_CHLOROPERLIDAE','Juga','ARACHNIDA'), collapse = "|")

# Filter out the taxa meeting the pattern above ^
final_table3<-final_table3[final_table3$taxon %in% taxa_notraits_rare$taxon==F,]

#boxplot(final_table3$failed_ratio, at =2,xlim=c(1,10),ylim=c(0,7))
#boxplot(final_table3$Regional_Ratio,at=4,add=T)
#boxplot(Dropped_Taxa_results$Regional_Ratio,at=6,add=T)
#boxplot(Dropped_Taxa_results$FRatio,at=8,add=T)
#boxplot(Ref_FoFe,at=8,add=T)
#axis(side=1, labels=c('prob. fail','prob. pass','drop.','reference'),at=c(2,4,6,8))
#savp(10,8,'C://Users//andrew.caudillo.BUGLAB-I9//Box//NAMC//Research Projects//AIM//IncreaserDecreaser_OE//FoFe_compare_all_tall.png')
 #site_pca_dat<-All_data[,(1:(ncol(All_data)-5))]




#write.csv(filtered_final,'C://Users//andrew.caudillo.BUGLAB-I9//Box//NAMC//Research Projects//AIM//IncreaserDecreaser_OE//OTU_and_response_table_w_traits.csv')
#create the response table, but remove the "overall" ratios etc.
#because now we have created the land use-specific ratios
#Response_table=plyr::join(final_table,final_table3,'right',by='taxon')
#removing the irrelevant parts now
#Response_table=Response_table[,-c(10:13)]
#renaming columns
#names(Response_table)[3]<-'Land_use'
#names(Response_table)[4]<-'Local_land_use'

# #creating the new response table, grouped by land use and local land use
# Response_table=Response_table %>%
#   left_join(n_sites_tbl, by = c("Ecoregion", "Land_use"))
# Response_table<-Response_table %>%
#   left_join(n_sites_local_tbl, by=c('Ecoregion','Local_land_use'))
#
# #creating local response table
# LResponses=Response_table%>%
#   select(Ecoregion, Local_land_use, taxon, prob_p_sum, ref_p_sum)
# #collapsing table to get Pcs for the local land uses only
# LResponses=LResponses %>%
#   group_by(Ecoregion, Local_land_use, taxon) %>%
#   summarise(
#     Of = sum(prob_p_sum, na.rm = TRUE),
#     Ef = first(ref_p_sum),
#     .groups = "drop"
#   )
# #removing irrelevant columns
# Response_table<-Response_table[,-c(4,(ncol(Response_table)),(ncol(Response_table)-1))]
# #same as the LResponses, but for watershed level
# Responses=Response_table %>%
#   group_by(Ecoregion, Land_use, taxon) %>%
#   summarise(
#     Of = sum(prob_p_sum, na.rm = TRUE),
#     Ef = first(ref_p_sum),   # if duplicated, otherwise take first()
#     .groups = "drop"
#   )
# #this section is compute variance, which can affect the ratio, but
# #not needed for regional scale
# if(1){
# # Response_table <- Response_table %>%
# #   mutate(
# #     p_ref = ref_p_sum / n_sites.x,
# #     var = n_sites.x * p_ref * (1 - p_ref),
# #     z = (prob_p_sum - n_sites.x * p_ref) / sqrt(var),
# #     Response = case_when(
# #       z > 1.96  ~ "increaser",
# #       z < -1.96 ~ "decreaser",
# #       TRUE      ~ "neutral"
# #     )
# #   )
# # Response_table <- Response_table %>%
# #   mutate(
# #     Lp_ref = ref_p_sum / n_sites.y,
# #     Lvar = n_sites.y * Lp_ref * (1 - Lp_ref),
# #     Lz = (prob_p_sum - n_sites.y * Lp_ref) / sqrt(Lvar),
# #     LResponse = case_when(
# #       Lz > 1.96  ~ "increaser",
# #       Lz < -1.96 ~ "decreaser",
# #       TRUE      ~ "neutral"
# #     )
# #   )
# # Response_table_clean=Response_table[!is.nan(Response_table$z),]
# # Response_table_clean=Response_table_clean[!is.nan(Response_table_clean$Lz),]
# #Response_table_clean=Response_table[,-c(6:7)]
# # Response_table_clean=Response_table_clean %>%
# #   relocate(c('var','z','Lvar','Lz','Response','LResponse'), .after = ref_p_sum)
# }
#isolating the traits only for Gower distance
traits_only <- final_table3[,c(2,13:ncol(final_table3))]

#taxa_to_drop<-wide_trait_table_OTUs$taxon[wide_trait_table_OTUs$taxon %in% traits_only$taxon==F]

traits_only<-traits_only[!duplicated(traits_only$taxon),]
row.names(traits_only)<-traits_only$taxon
#traits_only<-traits_only[traits_only$taxon %in% Ref_Results$taxon,]



traits_only2<-as.data.frame(lapply(traits_only[,2:ncol(traits_only)],factor))
traits_num<-as.data.frame(lapply(traits_only,as.numeric))
row.names(traits_only2)<-traits_only$taxon
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
                                                    'AdultFlyingStrength_abbrev_Weak',
                                                    'Emerge_synch_abbrev_Poorly',
                                                    'Shape_not_streamline','Survive_desiccation_no','Adult_exit_absent',
                                                    'Female_disp_abbrev_High','Swim_none',
                                                    'Lifespan_long')==F]
#'Swim_strong','Occurance_drift_common')==F]
#calculate gower distance on just traits
Gower_dist<-cluster::daisy(traits_only2,metric='gower')
#function that will find variable contribution
#since Gower distance does not readily give that info
gower_var_contrib <- function(var) {
  # if trait has < 2 non-missing values, then skip
  if (sum(!is.na(traits_only2[[var]])) < 2) return(NA_real_)

  # single-variable Gower distance
  d_var <- cluster::daisy(traits_only2[var], metric = "gower")

  # flatten and compute correlation, ignoring NA pairs
  cor(as.vector(Gower_dist), as.vector(d_var), use = "pairwise.complete.obs")
}

#get names of traits
vars <- names(traits_only2)
# correlate each variable’s dissimilarity with the full one
var_contrib <- sapply(vars, gower_var_contrib)
var_contrib=sort(var_contrib, decreasing = TRUE)
#these top 5 contribute the most to the dissimilarity of traits
top5<-names(var_contrib)[1:5]
top5


imp_traits <- traits_only2[,names(traits_only2) %in% top5]

#convert the Gower distance into a matrix for readability
 Gower_mat<-as.matrix(Gower_dist)
 #assign the taxa names
 row.names(Gower_mat)<-row.names(traits_only)
 #perform PCoA
 PCOA<-ape::pcoa(Gower_mat)

#smaller Gower
top5_Gower_dist<-cluster::daisy(imp_traits,metric='gower')
top5_Gower_mat<-as.matrix(top5_Gower_dist)
row.names(top5_Gower_mat)<-row.names(traits_only)
top5_PCOA<-ape::pcoa(top5_Gower_dist)

top5_scores<-as.data.frame(PCOA$vectors[,c(1,2)])
#assign scores to its own df
PCOA_Scores<-as.data.frame(PCOA$vectors[,c(1,2)])
#assign taxa names
PCOA_Scores$taxon<-row.names(traits_only)
top5_scores$taxon<-row.names(traits_only)

Responses_and_scores=plyr::join(All_Results,PCOA_Scores,by='taxon','left')
All_responses_and_scores=plyr::join(All_ratios,PCOA_Scores,by='taxon','left')
Responses_and_scores=Responses_and_scores[Responses_and_scores$taxon %in% row.names(traits_only2),]
ggplot(Responses_and_scores,
       aes(Axis.1, Axis.2, color=Regional_Response)) +
  geom_point(size=3, alpha=0.7) +
  scale_color_manual(values=c("Increaser"="orange",
                              "Decreaser"="dodgerblue3",
                              "Neutral"='red',
                              "Not expected"='black'))+
  facet_wrap('status')+guides(color=guide_legend(title='Response'))

savp(10,8,'C://Users//andrew.caudillo.BUGLAB-I9//Box//NAMC//Research Projects//AIM//IncreaserDecreaser_OE//regional_responses_traitspace_MRF.png')
#creating separate table for just local land use
#with a response of Increaser or Decreaser
# LResponses=LResponses %>%
# mutate(
#   ratio = ifelse(Ef > 0, Of / Ef, NA_real_),
#   response = case_when(
#     ratio > 1.2 ~"Increaser",
#     ratio < 0.8  ~"Decreaser",
#     ratio >0.8 & ratio < 1.2 ~"Neutral",
#     is.na(ratio)~"Not expected"
#   )
# )
# #doing the same for WS level land use
# Responses=Responses %>%
#   mutate(
#     ratio = ifelse(Ef > 0, Of / Ef, NA_real_),
#     response = case_when(
#       ratio > 1.2 ~"Increaser",
#       ratio < 0.8  ~"Decreaser",
#       ratio >0.8 & ratio < 1.2 ~"Neutral",
#       is.na(ratio)~"Not expected"
#     )
#   )
# #create metadata to join to PCoA for both land uses
# WSmetadata_all = Responses[,c('taxon','response','Ecoregion','Land_use')]
# Lmetadata_all=LResponses[,c('taxon','response','Ecoregion','Local_land_use')]
# #now joining the tables together for each spatial scale
# pcoa_plot_dat=left_join(PCOA_Scores,WSmetadata_all,by='taxon')
# Lpcoa_plot_dat=left_join(PCOA_Scores,Lmetadata_all,by='taxon')


#general trend of traits across land uses at watershed scale
#plotting distinct because of the duplicated records that come with the ecoregion,
#land use, etc.
All_responses_and_scores<-All_responses_and_scores[All_responses_and_scores$Local_land_use!='none',]
All_responses_and_scores$Local_land_use_Grp=ifelse(All_responses_and_scores$Local_land_use %like% 'Graz','Grazing',
                                                   'Logging')

Local_responses_and_scores = All_responses_and_scores %>%
  group_by(taxon,Local_land_use_Grp) %>%

  summarise(PFO=sum(pFo),
         PFE=sum(pFe),
         across(-c(pFo,pFe),first))

 # select(PFO,PFE,taxon,Axis.1,Axis.2,Local_land_use_Grp)

Local_responses_and_scores$PRATIO=Local_responses_and_scores$PFO / Local_responses_and_scores$PFE
Local_responses_and_scores$PRESPONSE=ifelse(Local_responses_and_scores$PRATIO >1.2, "Increaser",
                                                ifelse( Local_responses_and_scores$PRATIO > 0.8 & Local_responses_and_scores$PRATIO < 1.2, 'Neutral',"Decreaser"))


Watershed_responses_and_scores = All_responses_and_scores %>%
  group_by(taxon,Land_use) %>%
  summarise(PFO=sum(pFo),
            PFE=sum(pFe),
            across(-c(pFo,pFe),first))

Watershed_responses_and_scores$PRATIO=Watershed_responses_and_scores$PFO / Watershed_responses_and_scores$PFE
Watershed_responses_and_scores$PRESPONSE=ifelse(Watershed_responses_and_scores$PRATIO >1.2, "Increaser",
                                            ifelse( Watershed_responses_and_scores$PRATIO > 0.8 & Watershed_responses_and_scores$PRATIO < 1.2, 'Neutral',"Decreaser"))
Ecos_responses_and_scores=All_responses_and_scores %>%
  group_by(taxon,Ecoregion) %>%
  summarise(PFO=sum(pFo),
            PFE=sum(pFe),
            FO=sum(rFo),
            FE=sum(rFe),
            across(-c(pFo,pFe,-rFo,-rFe),first))
Ecos_responses_and_scores$PRATIO=Ecos_responses_and_scores$PFO / Ecos_responses_and_scores$PFE
Ecos_responses_and_scores$PRESPONSE=ifelse(Ecos_responses_and_scores$PRATIO >1.2, "Increaser",
                                                ifelse( Ecos_responses_and_scores$PRATIO > 0.8 & Ecos_responses_and_scores$PRATIO < 1.2, 'Neutral',"Decreaser"))
Ecos_responses_and_scores$RATIO=Ecos_responses_and_scores$FO / Ecos_responses_and_scores$FE
Ecos_responses_and_scores$RESPONSE=ifelse(Ecos_responses_and_scores$RATIO >1.2, "Increaser",
                                           ifelse( Ecos_responses_and_scores$RATIO > 0.8 & Ecos_responses_and_scores$RATIO < 1.2, 'Neutral',"Decreaser"))
 ggplot(Ecos_responses_and_scores %>% distinct(taxon, Axis.1, Axis.2, RESPONSE,Ecoregion),
       aes(Axis.1, Axis.2, color=RESPONSE)) +
  geom_point(size=3, alpha=0.7) +
  scale_color_manual(values=c("Increaser"="orange",
                              "Decreaser"="dodgerblue3",
                              "Neutral"='red',
                              "Not expected"='red'))+
  facet_wrap(~Ecoregion)+
  labs(color='RESPONSE')+guides(color=guide_legend(title='Response'))
savp(20,10,'C://Users//andrew.caudillo.BUGLAB-I9//Box//NAMC//Research Projects//AIM//IncreaserDecreaser_OE//Ref_Ecoregion_response_MRF.png')


RP=as.data.frame(Ecos_responses_and_scores[Ecos_responses_and_scores$Ecoregion=='Rangeland Plains',])

#
# # Convert site_taxa to presence/absence
# site_taxa_pa <- (site_taxa > 0) * 1
#
# # Site × trait (unweighted)
# site_trait_pa <- site_taxa_pa %*% taxa_traits
#
# # Distance and PCoA
# dist_pa <- vegdist(site_trait_pa, method = "bray")
# pcoa_pa <- cmdscale(dist_pa, eig = TRUE, k = 2)
#
# # Combine for plotting
# plot(pcoa_w$points[,1], pcoa_w$points[,2], col="blue", pch=16,
#      xlab="PCoA1", ylab="PCoA2", main="Weighted vs. Presence/Absence")
# points(pcoa_pa$points[,1], pcoa_pa$points[,2], col="red", pch=17)
# legend("topright", legend=c("Weighted","PA"), col=c("blue","red"), pch=c(16,17))


regional_and_traits=plyr::join(All_Results,traits_only[,names(traits_only) %in% c(top5,'taxon')])
Os=Os[,colnames(Os) %in% taxa_notraits_rare$taxon==F]
O_long=as.data.frame(t(Os))
O_long$taxon=names(Os)
ProbOs=ProbOs[,colnames(ProbOs) %in% taxa_notraits_rare$taxon==F]
ProbOs=ProbOs[row.names(ProbOs) %in% failed_sites$sampleId==F,]
ProbOs_long=as.data.frame(t(ProbOs))
ProbOs_long$taxon=names(ProbOs)
ProbOs_long$status='Probabilistic'
O_long$status='Reference'
O_long=plyr::join(O_long,traits_only[,names(traits_only) %in% c(top5,'taxon')])
ProbOs_long=plyr::join(ProbOs_long,traits_only[,names(traits_only) %in% c(top5,'taxon')])
site_cols=colnames(O_long)[1:(ncol(O_long)-7)]
traitcols=colnames(O_long)[(ncol(O_long)-4):ncol(O_long)]
Psite_cols=colnames(ProbOs_long)[1:(ncol(ProbOs_long)-7)]
Ptraitcols=colnames(ProbOs_long)[(ncol(ProbOs_long)-4):ncol(ProbOs_long)]
O_long=as.data.frame(O_long %>%
                       pivot_longer(
                         cols = all_of(site_cols),
                         names_to = "Site",
                         values_to = "Presence"
                       ) %>%
                       summarise(
                         across(
                           .cols = all_of(traitcols),
                           .fns  = ~ sum(.x * Presence, na.rm = TRUE) / sum(Presence, na.rm = TRUE),
                           .names = "{.col}_mean"
                         ),
                         n_present = sum(Presence, na.rm = TRUE)
                       ))
P_long=as.data.frame(ProbOs_long %>%
                       pivot_longer(
                         cols = all_of(Psite_cols),
                         names_to = "Site",
                         values_to = "Presence"
                       ) %>%
                       summarise(
                         across(
                           .cols = all_of(Ptraitcols),
                           .fns  = ~ sum(.x * Presence, na.rm = TRUE) / sum(Presence, na.rm = TRUE),
                           .names = "{.col}_mean"
                         ),
                         n_present = sum(Presence, na.rm = TRUE)
                       ))
#write.csv(regional_and_traits,'C://Users//andrew.caudillo.BUGLAB-I9//Box//NAMC//Research Projects//AIM//IncreaserDecreaser_OE//regional_responses_and_traits.csv')

O_long$status='Reference'
P_long$status='Probabilistic'
rbind(O_long,P_long)

clipr::write_clip(rbind(O_long,P_long))
