#univariate random forest/WestWide model ! Using model object, the actual ref sites.
#compare the Fo and Fe to the probabilistic sites' Fo and Fe. Maybe a better discrimination
#between site type and more trend in increaser/decreasers?
library(tidyr)
library(dplyr)
savp = function(W,H,fn) {
  dev.copy(dev=png,file=fn,wi=W,he=H,un="in",res=650)
  dev.off()
}
#read in ref data
ben_dat=read.csv("C://Users//andrew.caudillo.BUGLAB-I9//Box//NAMC//Research Projects//AIM//IncreaserDecreaser_OE//Updated_w_modelObj//WestWide_Ref_PA_260209.csv")

names(ben_dat)[1]<-'SiteCode'
row.names(ben_dat)=ben_dat$SiteCode;ben_dat=ben_dat[,-1]
ben_dat_raw=ben_dat
#omit Carabidae and Curuclionidae - not aquatic.
ben_dat<-ben_dat[,names(ben_dat) %in% c('Carabidae','Curculionidae')==F]
ben_dat_trim=ben_dat[,names(ben_dat) %in% 'Tanypodinae.1'==F]
ben_dat=ben_dat%>%
        mutate(Tanypodinae = ifelse(Tanypodinae+Tanypodinae.1>1,1,Tanypodinae+Tanypodinae.1)) %>%
        select(-Tanypodinae.1)
Pcs_raw=read.csv('C://Users//andrew.caudillo.BUGLAB-I9//Box//NAMC//Research Projects//AIM//IncreaserDecreaser_OE//Updated_w_modelObj//Ref_Pcs_orig_WW_raw.csv')
row.names(Pcs_raw)<-Pcs_raw$sampleId;Pcs_raw=Pcs_raw[,-c(1,2)]
Pcs_trim=read.csv('C://Users//andrew.caudillo.BUGLAB-I9//Box//NAMC//Research Projects//AIM//IncreaserDecreaser_OE//Updated_w_modelObj//Ref_Pcs_orig_WW_droppedTany.csv')
row.names(Pcs_trim)=Pcs_trim$sampleId;Pcs_trim=Pcs_trim[,-c(1:2)]
Pcs_trim=Pcs_trim[,names(Pcs_trim) %in% c('Carabidae','Curculionidae')==F]
Pcs_sum=read.csv('C://Users//andrew.caudillo.BUGLAB-I9//Box//NAMC//Research Projects//AIM//IncreaserDecreaser_OE//Updated_w_modelObj//Ref_Pcs_orig_WW_Tany_summed.csv')
row.names(Pcs_sum)=Pcs_sum$sampleId;Pcs_sum=Pcs_sum[,-c(1:2)]
Pcs_sum=Pcs_sum[,names(Pcs_sum) %in% c('Carabidae','Curculionidae')==F]
reorder_index=match(rownames(Pcs_sum), rownames(ben_dat))
ben_dat=ben_dat[reorder_index,]
reorder_index=match(rownames(Pcs_raw),rownames(ben_dat_raw))
ben_dat_raw=ben_dat_raw[reorder_index,]

reorder_index=match(rownames(Pcs_trim), rownames(ben_dat_trim))
ben_dat_trim=ben_dat_trim[reorder_index,]

Fos_raw=colSums(ben_dat_raw[row.names(ben_dat_raw) %in% Xeric_P==F,])
Fos_raw=Fos_raw[sort(names(Fos_raw))]
Fes_raw=colSums(Pcs_raw[row.names(ben_dat_raw) %in% Xeric_P==F,])
Fes_raw=Fes_raw[sort(names(Fes_raw))]
raw_ratio=Fos_raw/Fes_raw
plot(Fes_raw,Fos_raw,xlab='Fe',ylab='Fo',
     main='Original Westwide Fo/Fe')
abline(0,1,col='red')
raw_OEx=OE_calc(results_data = Pcs_raw[
  (row.names(Pcs_raw) %in% Xeric_P),
  #!(names(Pcs_raw) %in% c("Carabidae", "Curculionidae"))
],
               PA=ben_dat_raw[
                 (row.names(ben_dat_raw) %in% Xeric_P),
                 #!(names(ben_dat_raw) %in% c("Carabidae", "Curculionidae"))
               ],threshold = 0.5)
sd(raw_OE$OtoE)
mean(raw_OE$OtoE)
raw_OEx=OE_calc(results_data = Pcs_raw[
  (row.names(Pcs_raw) %in% Xeric_P),
  #!(names(Pcs_raw) %in% c("Carabidae", "Curculionidae"))
],
PA=ben_dat_raw[
  (row.names(ben_dat_raw) %in% Xeric_P),
  #!(names(ben_dat_raw) %in% c("Carabidae", "Curculionidae"))
],threshold = 0.5)
#calculating Fo and Fe
Fos=colSums(ben_dat)
Fos_trim=colSums(ben_dat_trim)
Fos_trim=Fos_trim[sort(names(Fos_trim))]
Fos=Fos[sort(names(Fos))]
Fes=colSums(Pcs_sum)
Fes_trim=colSums(Pcs_trim)
Fes_trim=Fes_trim[sort(names(Fes_trim))]
Fes=Fes[sort(names(Fes))]
names(Fos)==names(Fes)
names(Fos_trim)==names(Fes_trim)
ratio=Fos/Fes
sort(ratio)
ratio_dat=data.frame(Fo=Fos,Fe=Fes)
ratio_dat[ratio_dat$Fo > 300 & ratio_dat$Fo < 400 & ratio_dat$Fe >=300,]
graphics.off()
boxplot(ratio,ylab='Fo/Fe ratio')
#savp(10,8,'C://Users//andrew.caudillo.BUGLAB-I9//Box//NAMC//Research Projects//AIM//IncreaserDecreaser_OE//MRF_OE//Fo_Fe_boxplot_251114.png')
plot(Fes,Fos,ylab='Fo',xlab='Fe',main='Reference Fo vs Fe',pch=21,bg=rgb(0,0,1,alpha=0.5))
points(Fes_trim,Fos_trim,pch=21,bg=rgb(1,0,0,alpha=0.5))
abline(0,1,lty=2,lwd=1)
#savp(10,8,'C://Users//andrew.caudillo.BUGLAB-I9//Box//NAMC//Research Projects//AIM//IncreaserDecreaser_OE//MRF_OE//Ref_FovsFe_251114.png')
Fos_trim[names(Fos_trim)=='Tanypodinae'] / Fes_trim[names(Fes_trim)=='Tanypodinae']
Fos[names(Fos)=='Tanypodinae'] / Fes[names(Fes)=='Tanypodinae']
summed_tany=OE_calc(results_data = Pcs_sum[row.names(Pcs_sum)%in% Xeric_P==F,],
        PA=ben_dat[row.names(ben_dat) %in% Xeric_P==F,],
        threshold=0.5)
dropped_tany=OE_calc(results_data = Pcs_trim[row.names(Pcs_trim)%in% Xeric_P==F,],
                     PA=ben_dat_trim[row.names(ben_dat_trim) %in% Xeric_P==F,],
                     threshold=0.5)
sd(summed_tany$OtoE)
sd(dropped_tany$OtoE)

boxplot(raw_OE$OtoE,at=1,xlim=c(0,4))
boxplot(summed_tany$OtoE,at=2,add=T)
boxplot(dropped_tany$OtoE,at=3,add=T)
axis(1,at=c(1:3),labels=c('Orig','Sum','Drop'))


test_dat=read.csv('C://Users//andrew.caudillo.BUGLAB-I9//Box//NAMC//Research Projects//AIM//IncreaserDecreaser_OE//MRF_Prob_Os.csv')
#test_dat$Tanypodinae = test_dat$Tanypodinae + test_dat$Tanypodinae.1;test_dat<-test_dat[,names(test_dat) %in% 'Tanypodinae.1'==F,]
test_dat<-test_dat[,-1]
names(test_dat)[1]<-'sampleId'
#test_trim=test_dat[,names(test_dat)!='Tanypodinae.1']
Prob_pcs=read.csv("C://Users//andrew.caudillo.BUGLAB-I9//Box//NAMC//Research Projects//AIM//IncreaserDecreaser_OE//Updated_w_modelObj//Prob_Pcs_orig_WW_droppedTany.csv")
failed_sites<-read.csv('C://Users//andrew.caudillo.BUGLAB-I9//Box//NAMC//Research Projects//AIM//IncreaserDecreaser_OE//failed_sites.csv')
failed_dat=test_dat[test_dat$sampleId %in% failed_sites$sampleId,]
test_dat<-test_dat[test_dat$sampleId %in% failed_sites$sampleId==F,]
Prob_pcs=Prob_pcs[Prob_pcs$sampleId %in% failed_sites$sampleId==F,]
#Prob_pcs_trim=Prob_pcs[,names(Prob_pcs)!='Tanypodinae.1']
# Prob_pcs<-Prob_pcs%>%
#   mutate(Tanypodinae = pmax(Tanypodinae, Tanypodinae.1, na.rm = TRUE)) %>%
#   select(-Tanypodinae.1)
test_dat<-test_dat[,names(test_dat) %in% c('Carabidae','Curculionidae')==F]
Prob_pcs=Prob_pcs[,names(Prob_pcs) %in% c('Carabidae','Curculionidae')==F]
rownames(test_dat) <- test_dat$sampleId;test_dat<-test_dat[,names(test_dat) %in% 'sampleId'==F]
rownames(Prob_pcs)<-Prob_pcs$sampleId;Prob_pcs=Prob_pcs[,-c(1,2)]
#write.csv(Pr_PA, 'C://Users//andrew.caudillo.BUGLAB-I9//Box//NAMC//Research Projects//AIM//IncreaserDecreaser_OE//MRF_Prob_Os.csv')


Xeric_P=c(171583, 171664, 171668, 171669, 171777, 173398, 173404,
          173407, 173408, 173409, 178667, 178671, 178676, 184764,
          184768, 185060, 185063, 185699, 185702, 187162, 187182,
          187607,187662, 187779, 187780, 187849, 190728, 190762,
          190916, 191507,
          #these are the reference sites used in the westwide model
          "87087", "HAWK-216", "HAWK-225", "HAWK-226", "HAWK-227", "HAWK-228", "HAWK-241",
          "HAWK-86",  "HAWK-89", "HAWK-92",  "WE-1036",  "WE-1047",  "WE-1055",  "WE-1079",
          "WE-1085",  "WE-21",    "WE-47",    "WE-866",  "WE-889",   "WE-905",
          'R6REM-5','WE-317')

ben_dat_oth=ben_dat[row.names(ben_dat) %in% Xeric_P==F,]
Pcs_sum_oth=Pcs_sum[row.names(Pcs_sum) %in% Xeric_P==F,]

ben_dat_X=ben_dat[row.names(ben_dat) %in% Xeric_P,]
Pcs_sum_X=Pcs_sum[row.names(Pcs_sum) %in% Xeric_P,]
Oth_Fo=colSums(ben_dat_oth)
Oth_Fo=Oth_Fo[sort(names(Oth_Fo))]
Oth_Fe=colSums(Pcs_sum_oth)
Oth_Fe=Oth_Fe[sort(names(Oth_Fe))]
X_Fo=colSums(ben_dat_X)
X_Fo=X_Fo[sort(names(X_Fo))]
X_Fe=colSums(Pcs_sum_X)
X_Fe=X_Fe[sort(names(X_Fe))]
Oth_ratio=Oth_Fo/Oth_Fe
X_ratio=X_Fo/X_Fe

Pr_PA_oth=test_dat[row.names(test_dat) %in% Xeric_P==F,]
pred_probs_oth=Prob_pcs[row.names(Prob_pcs) %in% Xeric_P==F,]
Pr_PA_X=test_dat[row.names(test_dat) %in% Xeric_P,]
pred_probs_X=Prob_pcs[row.names(Prob_pcs) %in% Xeric_P,]
PFos=colSums(test_dat)
PFes=colSums(Prob_pcs)
PFos=PFos[sort(names(PFos))]
PFes=PFes[sort(names(PFes))]
PFos_oth=colSums(Pr_PA_oth)
PFos_oth=PFos_oth[sort(names(PFos_oth))]
#PFos_trim=colSums(test_trim)
PFes_oth=colSums(pred_probs_oth)
PFes_oth=PFes_oth[sort(names(PFes_oth))]
Pratio=PFos/PFes

Pratio_oth=PFos_oth/PFes_oth
Pratio_oth=Pratio_oth[match(names(Fos), names(Pratio_oth))]

Prob_combined_FoFes=as.data.frame(cbind(PFos,PFes))
#Pratio_small=Pratio[names(Pratio) %in% taxa_notraits_rare$taxon==F]
#PFos_small=PFos[names(PFos) %in% taxa_notraits_rare$taxon==F]
#PFes_small=PFes[names(PFes) %in% taxa_notraits_rare$taxon==F]

#Poth_Fosmall=Poth_Fo[names(Poth_Fo) %in% taxa_notraits_rare$taxon==F]
#Poth_Fesmall=Poth_Fe[names(Poth_Fe) %in% taxa_notraits_rare$taxon==F]

Px_Fo=colSums(Pr_PA_X);Px_Fe=colSums(pred_probs_X)
names(Px_Fo)==names(Px_Fe)
Px_Fo=Px_Fo[sort(names(Px_Fo))]
Px_Fe=Px_Fe[sort(names(Px_Fe))]
#Px_Fosmall=Px_Fo[names(Px_Fo) %in% taxa_notraits_rare$taxon==F]
#Px_Fesmall=Px_Fe[names(Px_Fe) %in% taxa_notraits_rare$taxon==F]
#Poth_ratio=Poth_Fo/Poth_Fe
Px_ratio=Px_Fo/Px_Fe
#Poth_ratio_small=Poth_Fosmall/Poth_Fesmall
#Px_ratio_small=Px_Fosmall/Poth_Fesmall
names(PFos)==names(PFes)
names(PFos_oth)==names(PFes_oth)


plot(PFes,PFos,ylab='Fo',xlab='Fe',main='Prob. Fo vs Fe')

abline(0,1,col='red')


boxplot(Pratio,ylab='Fo/Fe ratio')
#savp(10,8,'C://Users//andrew.caudillo.BUGLAB-I9//Box//NAMC//Research Projects//AIM//IncreaserDecreaser_OE//MRF_OE//Updated_Ref//Prob_FoFe_ratio_260210.png')
#savp(10,8,'C://Users//andrew.caudillo.BUGLAB-I9//Box//NAMC//Research Projects//AIM//IncreaserDecreaser_OE//MRF_OE//Updated_Ref//Prob_FovsFe_box_260210.png')


boxplot(Pratio,at=2,xlim=c(1,4),ylab='Fo/Fe')
boxplot(ratio,at=3,xlim=c(1,4),add=T)

axis(side=1,at=c(2,3),labels = c('Prob','Train'))


savp(10,8,'C://Users//andrew.caudillo.BUGLAB-I9//Box//NAMC//Research Projects//AIM//IncreaserDecreaser_OE//MRF_OE//Updated_Ref//Box_compare_FoFes_260210.png')

boxplot(ratio,at=1,col='purple3',ylim=c(0,max(Pratio[is.finite(Pratio)])),xlim=c(0,3),
        ylab='Fo/Fe ratio')
boxplot(Pratio,at=2,col='yellow3',add=T)
points(x=rep(2, length(sort(Pratio[is.finite(Pratio)],decreasing = T)[1:3])),y=(sort(Pratio[is.finite(Pratio)],decreasing = T)[1:3]),bg=c('blue','red','orange'),pch=21)

legend('topright',
       leg=c('Reference',
             'Probabilistic'),
       pch=rep(22,2),
       pt.bg=c('purple4',
               'yellow3'),
       bty='n',
       cex=0.8)
legend('topleft',
       leg=c('Callibaties','Laccobius','Gyraulus'),
       pch=c(21,21),
       pt.bg=c('blue','red','orange'),
       bty='n',
       cex=0.8)
#savp(10,8,'C://Users//andrew.caudillo.BUGLAB-I9//Box//NAMC//Research Projects//AIM//IncreaserDecreaser_OE//MRF_OE//Updated_Ref//regional_boxplot_compare_260210_colored.png')
graphics.off()
boxplot(Oth_ratio,at=1,xlim=c(0,5),ylim=c(0,max(Px_ratio[is.finite(Px_ratio)])),col='purple3',ylab='Fo/Fe ratio')
boxplot(Pratio_oth,at=2,add=T,col='yellow3')
boxplot(X_ratio,at=3,add=T,col='purple3')
boxplot(Px_ratio,at=4,add=T,col='yellow3')
points(rep(2,length(sort(Pratio_oth[is.finite(Pratio_oth)],decreasing = T)[1:3])),y=(sort(Pratio_oth[is.finite(Pratio_oth)],decreasing = T)[1:3]),bg=c('blue','red','orange'),pch=21)
points(x=rep(4, length(sort(Px_ratio[is.finite(Px_ratio)],decreasing = T)[1:3])),y=(sort(Px_ratio[is.finite(Px_ratio)],decreasing = T)[1:3]),bg=c('blue','red','orange'),pch=21)
mtext(text = c("Other", "Xeric"), side = 1, line = 1, at = c(1.5, 3.5), cex = 1)
mtext(text=c('Ref.','Prob.','Ref.','Prob.'),side=1,line=2, at=c(1,2,3,4),cex=0.7)
legend('topleft',
       leg=c('Paracloeodes',
             'Callibaetis',
             'Gammarus'),
       pt.bg=c('blue','red','orange'),
       pch=rep(21,5),
       bty='n',
       cex=0.8,
       ncol=1)
# legend('topright',
#        leg=c('Reference',
#              'Probabilistic'),
#        pch=rep(22,2),
#        pt.bg=c('purple4',
#                'yellow3'),
#        bty='n',
#        cex=0.8)
savp(10,8,'C://Users//andrew.caudillo.BUGLAB-I9//Box//NAMC//Research Projects//AIM//IncreaserDecreaser_OE//MRF_OE//Updated_Ref//ecoregion_boxplot_compare_260218_colored.png')
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

OEs=OE_calc(results_data = Pcs,
        PA=ben_dat,
        threshold = 0.5)
row.names(OEs)=row.names(Pcs)
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
Oth_Fo=Oth_Fo[sort(names(Oth_Fo))]
Oth_Fe=Oth_Fe[sort(names(Oth_Fe))]

OthEco_plotdat=data.frame(Fo=Oth_Fo,Fe=Oth_Fe)
EastXer_plotdat=data.frame(Fo=X_Fo,Fe=X_Fe)
P_othplotdat=data.frame(Fo=PFos_oth,Fe=PFes_oth)
row.names(P_othplotdat)=ifelse(row.names(P_othplotdat)=='Drunella_coloradensis_flavilinea',
                               'Drunella coloradensis/flavilinea',row.names(P_othplotdat))
Px_plotdat=data.frame(Fo=Px_Fo,Fe=Px_Fe)
row.names(Px_plotdat)=ifelse(row.names(Px_plotdat)=='Drunella_coloradensis_flavilinea',
                               'Drunella coloradensis/flavilinea',row.names(Px_plotdat))

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
  xlim(0,max(PFos_oth))+ylim(0,max(PFos_oth))+
  annotate("text",
           x = -Inf, y = Inf, # Position at top-left corner
           label = 'Reference',
           hjust = 0, vjust = 1, # Justify text relative to corner
           size = 4, color = "black")
B<-ggplot(data=EastXer_plotdat,aes(y=Fo,x=Fe))+geom_point()+
  ylim(0,max(P_othplotdat$Fo))+xlim(0,max(P_othplotdat$Fe))+
  geom_abline(intercept = 0,slope = 1,col='red')+ggtitle('Eastern Xeric')+
  annotate("text",
           x = -Inf, y = Inf, # Position at top-left corner
           label = 'Reference',
           hjust = 0, vjust = 1, # Justify text relative to corner
           size = 4, color = "black")+
  lims(x=c(0,30),y=c(0,30))
#ggplot(data=combined_dat,aes(x=Fe,y=Fo))+geom_point()+geom_abline(intercept = 0,slope = 1,col='red')+
# xlim(0,max(Px_plotdat$Fe))
#savp(10,8,'C://Users//andrew.caudillo.BUGLAB-I9//Box//NAMC//Research Projects//AIM//IncreaserDecreaser_OE//MRF_OE//Updated_Ref//MRF_allsites_FoFe.png')

#getting highlight taxa for showing some inc/decs
P_othplotdat$col=ifelse(row.names(P_othplotdat) %in% c('Callibaetis','Cambaridae','Laccobius',
                                                       'Paracloeodes','Sciomyzidae'),
                        'red',NA)
Px_plotdat$col=ifelse(row.names(Px_plotdat) %in% c('Callibaetis','Cambaridae','Laccobius',
                                                   'Paracloeodes','Sciomyzidae'),
                      'red',NA)
Poth_highlight=P_othplotdat[which(P_othplotdat$col=='red'),]
Poth_highlight$taxon=row.names(Poth_highlight)
Px_highlight=Px_plotdat[which(Px_plotdat$col=='red'),]
Px_highlight$taxon=row.names(Px_highlight)
C<-ggplot(data=P_othplotdat,aes(y=Fo,x=Fe))+geom_point()+geom_abline(intercept = 0,slope = 1,col='red')+
  geom_point(data=Poth_highlight,aes(x=Fe,y=Fo,color=taxon))+
  ylim(0,max(P_othplotdat$Fo))+xlim(0,max(P_othplotdat$Fe))+
  scale_color_manual(values=c('Callibaetis' = 'red',
                              'Cambaridae' = 'dodgerblue',
                              'Laccobius' = 'purple',
                              'Paracloeodes' = 'orange',
                              'Sciomyzidae' = 'blue'))+
  annotate("text",
           x = -Inf, y = Inf, # Position at top-left corner
           label = 'Probabilistic',
           hjust = 0, vjust = 1, # Justify text relative to corner
           size = 4, color = "black")+
  theme(legend.position = "bottom")
D<-ggplot(data=Px_plotdat,aes(y=Fo,x=Fe))+geom_point()+geom_abline(intercept = 0,slope = 1,col='red')+
  ylim(0,30)+xlim(0,30)+# EastXer_plotdat$Fe))+
  geom_point(data=Px_highlight,aes(x=Fe,y=Fo,color=taxon))+
  scale_color_manual(values=c('Callibaetis' = 'red',
                              'Cambaridae' = 'dodgerblue',
                              'Laccobius' = 'purple',
                              'Paracloeodes' = 'orange',
                              'Sciomyzidae' = 'blue'))+
  annotate("text",
           x = -Inf, y = Inf, # Position at top-left corner
           label = 'Probabilistic',
           hjust = 0, vjust = 1, # Justify text relative to corner
           size = 4, color = "black")+
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
C<-ggplot(data=P_othplotdat,aes(y=Fo,x=Fe))+geom_point()+geom_abline(intercept = 0,slope = 1,col='red')+
  geom_point(data=Poth_highlight,aes(x=Fe,y=Fo,color=taxon))+
  ylim(0,max(P_othplotdat$Fo))+xlim(0,max(P_othplotdat$Fe))+
  scale_color_manual(values=c('Callibaetis' = 'red',
                              'Cambaridae' = 'dodgerblue',
                              'Laccobius' = 'purple',
                              'Paracloeodes' = 'orange',
                              'Sciomyzidae' = 'blue'))+
  annotate("text",
           x = -Inf, y = Inf, # Position at top-left corner
           label = 'Probabilistic',
           hjust = 0, vjust = 1, # Justify text relative to corner
           size = 4, color = "black")+
  theme(legend.position = "none")
#define the plot, but don't plot it
cmbin_plot=gridExtra::grid.arrange(A,B,C,D,ncol=2)
#plot the final graph with shared legend
gridExtra::grid.arrange(cmbin_plot,tax_legend,heights=c(10,1))

savp(10,8,'C://Users//andrew.caudillo.BUGLAB-I9//Box//NAMC//Research Projects//AIM//IncreaserDecreaser_OE//MRF_OE//Updated_Ref//ecoregions_FoFe_sitecompare_260224_Colored.png')
#get O/E scores for all sites / ecoregion subsets
ref_OEs=OE_calc(results_data = Pcs_sum,
                    PA=ben_dat,
                    threshold=threshold)
ref_OEs_other=OE_calc(results_data = Pcs_sum_oth,
                      PA=ben_dat_oth,
                      threshold=threshold)
pred_probs=pred_probs[match(names(PFos), names(pred_probs))]
Pr_PA=Pr_PA[match(names(PFos), names(Pr_PA))]
test_OEs=OE_calc(results_data=Pes,
                     PA=ProbOs,
                     threshold=threshold)
XerOE=OE_calc(results_data = Pcs_sum_X,
              PA=ben_dat_X,
              threshold=threshold)
OthEco_OE=OE_calc(results_data = Pcs_sum_oth,
                  PA=ben_dat_oth,
                  threshold=threshold)

POthOE=OE_calc(results_data = pred_probs_oth,
               PA=Pr_PA_oth,
               threshold = threshold)
PXOE=OE_calc(results_data = pred_probs_X,
             PA=Pr_PA_X,
             threshold = threshold)
# WW_calc=OE_calc(results_data = WW_E[,-1],
#                 PA=WW[,-1],
#                 threshold = threshold)
# WW_X_calc=OE_calc(results_data = WW_E_Xer[,-1],
#                 PA=WW_Xer[,-1],
#                 threshold = threshold)
sd(OthEco_OE$OtoE)
boxplot(OthEco_OE$OtoE,at=1,xlim=c(0,5),ylab='O/E score',col='purple4',ylim=c(0,2))
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
savp(10,8,'C://Users//andrew.caudillo.BUGLAB-I9//Box//NAMC//Research Projects//AIM//IncreaserDecreaser_OE//MRF_OE//Updated_Ref//OtoE_boxes_ecoregions_260212.png')
#savp(10,8,'C://Users//andrew.caudillo.BUGLAB-I9//Box//NAMC//Research Projects//AIM//IncreaserDecreaser_OE//MRF_OE//OtoE_boxes_Pc05_251106.png')
### This is just looking at O/E performance and metrics surrounding it
boxplot(ref_OEs$OtoE,at=1,xlim=c(0,3),col='purple4',ylim=c(0,max(test_OEs$OtoE)))
boxplot(test_OEs$OtoE,at=2,add=T,col='yellow3')
legend('topright',
       leg=c('Reference',
             'Probabilistic'),
       pch=rep(22,2),
       pt.bg=c('purple4',
               'yellow3'),
       bty='n',
       cex=0.8)
savp(10,8,'C://Users//andrew.caudillo.BUGLAB-I9//Box//NAMC//Research Projects//AIM//IncreaserDecreaser_OE//MRF_OE//Updated_Ref//OtoE_boxes_Pc0_260212.png')

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

boxplot(ratio,at=2,xlim=c(0,6),col='purple4',ylim=c(0,max(failedFoFe[is.finite(failedFoFe)])))
boxplot(Pratio,at=3,add=T,col='yellow3')
points(x=rep(3, length(sort(Pratio[is.finite(Pratio)],decreasing = T)[1:3])),y=(sort(Pratio[is.finite(Pratio)],decreasing = T)[1:3]),bg=c('blue','red','orange'),pch=21)
boxplot(failedFoFe,at=4,add=T,col='red')
points(x=rep(4, length(sort(failedFoFe[is.finite(failedFoFe)],decreasing = T)[1:3])),y=(sort(failedFoFe[is.finite(failedFoFe)],decreasing = T)[1:3]),bg=c('blue','purple','dodgerblue'),pch=21)
axis(side=1,at=c(2:4),labels=c('Ref.','Prob.','Failed'))
legend('topleft',
       leg=c('Callibaties','Laccobius','Gyraulus','Paracloeodes','Gammarus'),
       pch=rep(21,5),
       pt.bg=c('blue','red','orange','purple','dodgerblue'),
       bty='n',
       cex=0.8)

savp(10,8, 'C://Users//andrew.caudillo.BUGLAB-I9//Box//NAMC//Research Projects//AIM//IncreaserDecreaser_OE//MRF_OE//Updated_Ref//Ref_prob_failed_boxes.png')





#savp(10,8,'C://Users//andrew.caudillo.BUGLAB-I9//Box//NAMC//Research Projects//AIM//IncreaserDecreaser_OE//MRF_OE//OtoE_boxes_Pc0_251117.png')




#savp(10,8,'C://Users//andrew.caudillo.BUGLAB-I9//Box//NAMC//Research Projects//AIM//IncreaserDecreaser_OE//MRF_OE//OvsE_Pc0_251023.png')





#clipr::write_clip(MRF_WW_OEs)



#savp(10,8,'C://Users//andrew.caudillo.BUGLAB-I9//Box//NAMC//Research Projects//AIM//IncreaserDecreaser_OE//MRF_OE//MRF_OOB_per_response_box.png')


#testing model on known low OE score sites.
#optional, just testing for overfit.




#pois binomial distrib
#this gets the confidence intervals for the predictions
#quick and simple vectorized calculations
#copy/paste results into a spreadsheet with overall responses and
#should have enough info for a descriptive table.
ben_dat=ben_dat[order(names(ben_dat)),]
Pcs_sum=Pcs_sum[order(names(Pcs_sum)),]
ben_dat=ben_dat[,names(Fos)]
Pcs_sum=Pcs_sum[,names(Fos)]
Fes=Fes[names(Fos)]
names(ben_dat)==names(Fos)
test_dat=test_dat[,names(PFos)]
Prob_pcs=Prob_pcs[,names(PFos)]
PFes=PFes[names(PFos)]

ben_dat_oth=ben_dat_oth[,names(Pcs_sum_oth)]
Oth_Fo=Oth_Fo[names(Pcs_sum_oth)]
Oth_Fe=Oth_Fe[names(Pcs_sum)]

ben_dat_X=ben_dat_X[,names(Pcs_sum_X)]
X_Fo=X_Fo[names(Pcs_sum_X)]
X_Fe=X_Fe[names(Pcs_sum_X)]

CI_dat <- data.frame(
  Fe_mean = rep(NA, ncol(Pcs_sum_X)),
  LL      = rep(NA, ncol(Pcs_sum_X)),
  UL      = rep(NA, ncol(Pcs_sum_X))
)
row.names(CI_dat)=names(Pcs_sum_X)
n=nrow(Pcs_sum_X)
for(i in 1:ncol(Pcs_sum_X)){
  x=poibin::qpoibin(qq=c(0.025, 0.975),pp=Pcs_sum_X[,i])
  Fe_mean = sum(Pcs_sum_X[,i])/n
  ci_Fe=x/n
  CI_dat$LL[i]=(x[1]/n) * n
  CI_dat$UL[i]=(x[2]/n) * n
  CI_dat$Fe_mean[i]<-Fe_mean
  CI_dat$Fo[i]=X_Fo[i]
  #CI_dat$X[i]=x
}
CI_dat
clipr::write_clip(CI_dat)
