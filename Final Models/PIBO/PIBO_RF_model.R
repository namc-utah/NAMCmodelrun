

################################;

# PIBO model using RF.

# first, source the prediction script and also load the desired model;

library(cluster); library(Hmisc);
library(randomForest);
library(plyr);
library(dplyr);

   source("\\\\share1.bluezone.usu.edu\\miller\\\\buglab\\OE_Modeling\\NAMC_Supported_OEmodels\\PIBO\\PIBO_RF2020\\model.predict.RanFor.r")
   load("\\\\share1.bluezone.usu.edu\\miller\\\\buglab\\OE_Modeling\\NAMC_Supported_OEmodels\\PIBO\\PIBO_RF2020\\Benkendorf.RF.Model.Version1.Rdata")
   view_Rdata<- load("\\\\share1.bluezone.usu.edu\\miller\\\\buglab\\OE_Modeling\\NAMC_Supported_OEmodels\\PIBO\\PIBO_RF2020\\Benkendorf.RF.Model.Version1.Rdata")
   view_Rdata

   # User must supply a sample-by-taxa matrix of taxon presence/absence (coded as 1 or 0), for all new samples;
   # User must also supply a corresponding file of predictor data for those same samples;
   # These 2 files should have similar formats as the original taxa and predictor data sets used to build the model (see step 1 above);
   # Notes on format --
   #   A) The sample ID column in both files should be read into R as a row name (see Step 1 examples).
   #   B) Predictor data set -- Must include columns with the same names, units, etc.,
   #        as the model's predictor variables. All other columns will be ignored;
   #        Column order does not matter;
   #        Predictions and calculations of O/E will be made only for those samples that have;
   #        complete data for all model predictors.;
   #   C)  Sample-by-taxa matrix.  Sample ID's (row names) must match those of predictor data.
   #        OTU names (column names) must include all taxa for which the model can make predictions. Other taxa (columns) are ignored.
   #        The model can make predictions for all taxa in the calibration data (bugcal.pa) that were present at 1 or more calibration sites;
   #        To see a list of these taxa, do the following:
         names(bugall)[colSums(bugall)>0];
#########;
pred.test <- read.csv("\\\\share1.bluezone.usu.edu\\miller\\\\buglab\\OE_Modeling\\NAMC_Supported_OEmodels\\PIBO\\InputsAndResults_PIBO2009oe\\Current Run\\habitat.csv",header=T);
head(pred.test)
bug.test <- read.csv("\\\\share1.bluezone.usu.edu\\miller\\\\buglab\\OE_Modeling\\NAMC_Supported_OEmodels\\PIBO\\InputsAndResults_PIBO2009oe\\Current Run\\bugs.csv",header=T);
head(bug.test)
bug.otu <- read.csv("\\\\share1.bluezone.usu.edu\\miller\\\\buglab\\OE_Modeling\\NAMC_Supported_OEmodels\\PIBO\\InputsAndResults_PIBO2009oe\\Current Run\\OTUbugs1.csv",header=T);
head(bug.otu)

# merges test bug data with a null set of all OTU bug names in order for VanSickle code to work properly.
bug.test.pa <- merge(x=bug.otu, y=bug.test, all=TRUE)         
# converts na data produced from columns above into 0s;
bug.test.pa[is.na(bug.test.pa)]<- 0 

# Reorder sampleid in ascending order
bug.test.pa = bug.test.pa[order(bug.test.pa$SAMPLEID), ];
pred.test = pred.test[order(pred.test$SAMPLEID), ];
head(bug.test.pa)
              
# makes predictions for test data;
OE.final<-model.predict.RanFor(bugall, grps.final, preds.final, ranfor.mod=rf.mod, prednew=pred.test, bugnew.pa=bug.test.pa, Pc=0.5);

# Append Sampleid to results
OE.final$OE.scores = data.frame(bug.test.pa$SAMPLEID, OE.final$OE.scores)
names(OE.final$OE.scores) = c('SAMPLEID',names(OE.final$OE.scores)[2:ncol(OE.final$OE.scores)] )
OE.final$Capture.Probs = data.frame(bug.test.pa$SAMPLEID, OE.final$Capture.Probs)
names(OE.final$Capture.Probs) = c('SAMPLEID',names(OE.final$Capture.Probs)[2:ncol(OE.final$Capture.Probs)] )

head(OE.final)  
write.csv(OE.final$OE.scores,'\\\\share1.bluezone.usu.edu\\miller\\\\buglab\\OE_Modeling\\NAMC_Supported_OEmodels\\PIBO\\InputsAndResults_PIBO2009oe\\Current Run\\Stribnite_OE_Scores.csv')

#head(OE.final$Capture.Probs)  # Look at the predicted probabilties for each taxon at each site;
#write.csv(OE.final$Capture.Probs,'\\\\share1.bluezone.usu.edu\\miller\\\\buglab\\OE_Modeling\\NAMC_Supported_OEmodels\\PIBO\\InputsAndResults_PIBO2009oe\\Current Run\\CaptProbs.csv')#

##Ensuring that the bug.test.pa merge is working well  
##write.csv(bug.test.pa,'\\\\share1.bluezone.usu.edu\\miller\\\\buglab\\OE_Modeling\\NAMC_Supported_OEmodels\\PIBO\\InputsAndResults_PIBO2009oe\\Current Run\\bugtestpa.csv')


############################
#Environmental outlier test#
############################
#Uses watershed area (km2), ELVmean_WS (m), Pmax_WS (mm)


#Calculating distances between reference sites
CalPredsModelApplicability=read.csv('Z:\\buglab\\OE_Modeling\\NAMC_Supported_OEmodels\\PIBO\\ModelFiles_PIBO2009oe\\CalPredsModelApplicability.csv')


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

TestPredsModelApplicability=read.csv('Z:\\buglab\\OE_Modeling\\NAMC_Supported_OEmodels\\PIBO\\InputsAndResults_PIBO2009oe\\Current Run\\TestPredsModelApplicability.csv')


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
sample.list = row.names(pred.test)
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

final=cbind(OE.final$OE.scores,out.flag90,TestPredsModelApplicability,outlier.preds.test.std)

write.csv(final,'\\\\share1.bluezone.usu.edu\\miller\\\\buglab\\OE_Modeling\\NAMC_Supported_OEmodels\\PIBO\\InputsAndResults_PIBO2009oe\\Current Run\\PIBO_OE_Scores_modelapp.csv')


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
png(filename="PIBO_PCA.png")
fviz_pca_biplot(res.pca, repel = TRUE,
                col.var = "#696969", # Variables color
                col.ind = groups,
                geom=c("point"),
                addEllipses = TRUE, # Concentration ellipses
                ellipse.type = "norm",
                ellipse.level=0.90
)
dev.off()












