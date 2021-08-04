#read in test data

Boxdata <- ("\\\\share1.bluezone.usu.edu\\miller\\\\buglab\\OE_Modeling\\NAMC_Supported_OEmodels\\UTDEQ\\AllSeasonsModel_2015\\InputsAndResults\\CurrentRun")
setwd(Boxdata) ##I am having trouble re-setting wd so placed this section ahead of loading model, and then spelled out file paths for input and output
getwd()
prednew <- read.csv("UTDEQ_Habitat.csv",row.names="SAMPLE")
bugnew.raw <- read.csv("UTDEQ_Bugs.csv")#Calculate O/E for test data:

setwd("\\\\share1.bluezone.usu.edu\\miller\\\\buglab\\OE_Modeling\\NAMC_Supported_OEmodels\\UTDEQ\\AllSeasonsModel_2015\\")
load('Model files\\UTDEQ_15_OE_model.rdata')

library(randomForest)
library(labdsv)
library(vegan)

###########################################################################

#Matrify and resample to 400 counts if necessary- use matrify function in labdsv and rrarefy function in vegan
bugnew.matrix=matrify(bugnew.raw)
rarefy.vec=rowSums(bugnew.matrix)
rarefy.vec[rarefy.vec>400]=400#already done in NAMC export
bugnew.400cnt=rrarefy(bugnew.matrix,rarefy.vec)

#sort predictors to match bug data
predictors=prednew[row.names(bugnew.400cnt),]

#Convert bug data from matrix to dataframe for Van Sickle's O/E function to work.
bugnew.400cnt=data.frame(bugnew.400cnt)

ranfor.mod=rfmod.final
bugcal=bugcal.pa

#Make predictions

OE.test.5<-model.predict.RanFor.4.2(bugcal.pa,grps.final, preds.final,ranfor.mod, prednew,bugnew=bugnew.400cnt,Pc=0.5,Cal.OOB=FALSE)

write.csv(OE.test.5$OE.scores,'\\\\share1.bluezone.usu.edu\\miller\\\\buglab\\OE_Modeling\\NAMC_Supported_OEmodels\\UTDEQ\\AllSeasonsModel_2015\\InputsAndResults\\\\CurrentRun\\UTDEQ_OEscores.csv')


############################
#Environmental outlier test#
############################
#Uses watershed area (km2), ELVmean_WS (m), Pmax_WS (mm)


#Calculating distances between reference sites
CalPredsModelApplicability=read.csv('Z:\\buglab\\OE_Modeling\\NAMC_Supported_OEmodels\\UTDEQ\\AllSeasonsModel_2015\\Model files\\CalPredsModelApplicability.csv')


logSQ_KM=log(CalPredsModelApplicability$WsAreaSqKm)
boxplot(logSQ_KM)
ElevCat=CalPredsModelApplicability$ElevCat
boxplot(ElevCat)
logPrecip8110Ws=log(CalPredsModelApplicability$Precip8110Ws)
boxplot(logPrecip8110Ws)
Tmean8110Ws=CalPredsModelApplicability$Tmean8110Ws
boxplot(Tmean8110Ws)
LAT=CalPredsModelApplicability$LAT
LONG=CalPredsModelApplicability$LONG
outlier.preds.ref.raw=data.frame(logSQ_KM,ElevCat,logPrecip8110Ws,Tmean8110Ws, row.names=CalPredsModelApplicability$SiteID)
#outlier.preds.ref.raw=data.frame(logSQ_KM,ElevCat,logPrecip8110Ws,Tmean8110Ws,LAT,LONG, row.names=CalPredsModelApplicability$SiteID)
#standardize by min and max of ref data
outlier.preds.ref.std=matrix(nrow=dim(outlier.preds.ref.raw)[1], ncol=0)
for (n in 1:dim(outlier.preds.ref.raw)[2]){
  cur.var=outlier.preds.ref.raw[,n]
  cur.var.std=(cur.var-min(outlier.preds.ref.raw[,n]))/(max(outlier.preds.ref.raw[,n])-min(outlier.preds.ref.raw[,n]))
  outlier.preds.ref.std=cbind(outlier.preds.ref.std, cur.var.std)
}
colnames(outlier.preds.ref.std)=names(outlier.preds.ref.raw)
row.names(outlier.preds.ref.std)=row.names(outlier.preds.ref.raw)

#calculate distances between ref data and average distance of 10 closest sites:
ref.dist=as.matrix(dist(outlier.preds.ref.std, method="euclidean"))
ref.top10mean.dist=vector()
for (n in 1:dim(ref.dist)[1]){
  cur.dist=ref.dist[n,-n]
  cur.10dist=(cur.dist[order(cur.dist)])[1:10]
  mean10dist=mean(cur.10dist)
  ref.top10mean.dist=append(ref.top10mean.dist, mean10dist)
}

# #calculate distances between ref data and average distance of all sites:
# ref.dist=as.matrix(dist(outlier.preds.ref.std, method="euclidean"))
# ref.mean.dist=vector()
# for (n in 1:dim(ref.dist)[1]){
#   cur.dist=ref.dist[n,-n]
#   #cur.10dist=(cur.dist[order(cur.dist)])[1:10]
#   meandist=mean(cur.dist)
#   ref.mean.dist=append(ref.mean.dist, meandist)
# }

TestPredsModelApplicability=read.csv('Z:\\buglab\\OE_Modeling\\NAMC_Supported_OEmodels\\UTDEQ\\AllSeasonsModel_2015\\InputsAndResults\\CurrentRun\\TestPredsModelApplicability.csv')


logSQ_KM=log(TestPredsModelApplicability$WsAreaSqKm)
boxplot(logSQ_KM)
ElevCat=TestPredsModelApplicability$ElevCat
boxplot(ElevCat)
logPrecip8110Ws=log(TestPredsModelApplicability$Precip8110Ws)
boxplot(logPrecip8110Ws)
Tmean8110Ws=TestPredsModelApplicability$Tmean8110Ws
boxplot(Tmean8110Ws)
LAT=TestPredsModelApplicability$LAT
LONG=TestPredsModelApplicability$LONG
outlier.preds.test.raw=data.frame(logSQ_KM,ElevCat,logPrecip8110Ws,Tmean8110Ws, row.names=TestPredsModelApplicability$SampleID)
#outlier.preds.test.raw=data.frame(logSQ_KM,ElevCat,logPrecip8110Ws,Tmean8110Ws,LAT,LONG, row.names=TestPredsModelApplicability$SampleID)
#standardize by min and max of ref data
outlier.preds.test.std=matrix(nrow=dim(outlier.preds.test.raw)[1], ncol=0)
for (n in 1:dim(outlier.preds.test.raw)[2]){
  cur.var=outlier.preds.test.raw[,n]
  cur.var.std=(cur.var-min(outlier.preds.ref.raw[,n]))/(max(outlier.preds.ref.raw[,n])-min(outlier.preds.ref.raw[,n]))
  outlier.preds.test.std=cbind(outlier.preds.test.std, cur.var.std)
}
colnames(outlier.preds.test.std)=names(outlier.preds.ref.raw)
row.names(outlier.preds.test.std)=row.names(outlier.preds.test.raw)

#calculate distances between test and ref sites and calculate average distance of 10 closest sites:
testdist.data=rbind(outlier.preds.ref.std, outlier.preds.test.std)
test.dist=as.matrix(dist(testdist.data, method="euclidean"))
test.dist=test.dist[row.names(outlier.preds.test.std),row.names(outlier.preds.ref.std)]
test.top10mean.dist=vector()
out.flag90=vector()
sample.list = row.names(prednew)
for (n in 1:length(sample.list)){
  if(length(sample.list)>1){cur.dist=test.dist[n,]}
  if(length(sample.list)==1){cur.dist=test.dist}
  cur.10dist=(cur.dist[order(cur.dist)])[1:10]
  mean10dist=mean(cur.10dist)
  test.top10mean.dist=append(test.top10mean.dist, mean10dist)
  if (mean10dist>quantile(ref.top10mean.dist,0.90)) flag90="Yes" else flag90="No"
  out.flag90=append(out.flag90, flag90)
}

# #calculate distances between test and ref sites and calculate average distance of all sites:
# testdist.data=rbind(outlier.preds.ref.std, outlier.preds.test.std)
# test.dist=as.matrix(dist(testdist.data, method="euclidean"))
# test.dist=test.dist[row.names(outlier.preds.test.std),row.names(outlier.preds.ref.std)]
# test.mean.dist=vector()
# out.flag90=vector()
# sample.list = row.names(test.preds)
# for (n in 1:length(sample.list)){
#   if(length(sample.list)>1){cur.dist=test.dist[n,]}
#   if(length(sample.list)==1){cur.dist=test.dist}
#   #cur.10dist=(cur.dist[order(cur.dist)])[1:10]
#   meandist=mean(cur.dist)
#   test.mean.dist=append(test.mean.dist, meandist)
#   if (meandist>quantile(ref.mean.dist,0.90)) flag90="Yes" else flag90="No"
#   out.flag90=append(out.flag90, flag90)
# }

final=cbind(OE.test.5$OE.scores,out.flag90,TestPredsModelApplicability,outlier.preds.test.std)

write.csv(final,'\\\\share1.bluezone.usu.edu\\miller\\\\buglab\\OE_Modeling\\NAMC_Supported_OEmodels\\UTDEQ\\AllSeasonsModel_2015\\InputsAndResults\\CurrentRun\\UTDEQ_OEscores.csv')


####double check the results visually 
library("factoextra")
outlier.preds.ref.std2=as.data.frame(outlier.preds.ref.std)
outlier.preds.test.std2=as.data.frame(outlier.preds.test.std)
outlier.preds.ref.std2$Status="Ref"
outlier.preds.test.std2$Status=out.flag90
modelappinput=rbind(outlier.preds.test.std2,outlier.preds.ref.std2)
modelappinputactive=modelappinput[,c(1:4)]
res.pca <- prcomp(modelappinputactive, scale=FALSE)
fviz_eig(res.pca) #scree plot
groups=as.factor(modelappinput$Status)#create group variable to visualize formal test above
#create PCA plot in R with 90% confidence ellipses around each group
png(filename="UTDEQ_PCA.png")
fviz_pca_biplot(res.pca, repel = TRUE,
                col.var = "#696969", # Variables color
                col.ind = groups,
                geom=c("point"),
                addEllipses = TRUE, # Concentration ellipses
                ellipse.type = "norm",
                ellipse.level=0.90
)
dev.off()
##This is the location of the new PCA biplot
getwd()










