##This code is used to determine the model appropriateness for a given site.
##Uses closest 10 ref sites for analysis
##Considers elevation, watershed area, precipitation, temperature

ModelApplicability<- function(CalPredsModelApplicability, modelId, prednew){

#CalPredsModelApplicability=read.csv("C:/Users/jenni/Box/NAMC (Trip Armstrong)/OE_Modeling/NAMC_Supported_OEModels/Model Applicability/CalPredsModelApplicability.csv")

##subset calibration set to only include sites used in the appropriate model build
CalPreds <- subset(CalPredsModelApplicability,modelID==modelId)
logSQ_KM=log(CalPreds$WsAreaSqKm)
ElevCat_standardized=CalPreds$ElevCat
logPrecip8110Ws=log(CalPreds$Precip8110Ws)
Tmean8110Ws_standardized=CalPreds$Tmean8110Ws

outlier.preds.ref.raw=data.frame(logSQ_KM,ElevCat_standardized,logPrecip8110Ws,Tmean8110Ws_standardized, row.names=CalPreds$SiteID)

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


logSQ_KM=log(applicabilitypreds$WsAreaSqKm)
ElevCat_standardized=applicabilitypreds$ElevCat
logPrecip8110Ws=log(applicabilitypreds$Precip8110Ws)
Tmean8110Ws_standardized=applicabilitypreds$Tmean8110Ws

outlier.preds.test.raw=data.frame(logSQ_KM,ElevCat_standardized,logPrecip8110Ws,Tmean8110Ws_standardized, row.names=applicabilitypreds$sampleId)

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
sample.list = row.names(outlier.preds.test.raw)
for (n in 1:length(sample.list)){
  if(length(sample.list)>1){cur.dist=test.dist[n,]}
  if(length(sample.list)==1){cur.dist=test.dist}
  cur.10dist=(cur.dist[order(cur.dist)])[1:10]
  mean10dist=mean(cur.10dist)
  test.top10mean.dist=append(test.top10mean.dist, mean10dist)
  if (mean10dist>quantile(ref.top10mean.dist,0.90)) flag90="Yes" else flag90="No"
  out.flag90=append(out.flag90, flag90)
}

final=cbind(out.flag90,applicabilitypreds,outlier.preds.test.std)
final$ModelApplicability=ifelse(out.flag90=="Yes","Fail","Pass")
ModelApplicability=ifelse(out.flag90=="Yes","Fail","Pass")
return(final)
}

#
# ##the following section creates visuals and is not needed for applicability determination
# ## this does not need to be run every time
# ####double check the results visually
# library("factoextra")
# outlier.preds.ref.std2=as.data.frame(outlier.preds.ref.std)
# outlier.preds.test.std2=as.data.frame(outlier.preds.test.std)
# outlier.preds.ref.std2$Status="Ref"
# outlier.preds.test.std2$Status=out.flag90
# modelappinput=rbind(outlier.preds.test.std2,outlier.preds.ref.std2)
# modelappinputactive=modelappinput[,c(1:4)]
# res.pca <- prcomp(modelappinputactive, scale=FALSE)
# fviz_eig(res.pca) #scree plot
# groups=as.factor(modelappinput$Status)#create group variable to visualize formal test above
# #create PCA plot in R with 90% confidence ellipses around each group
# png(filename="OR_MWCF_PCA.png")
# fviz_pca_biplot(res.pca, repel = TRUE,
#                 col.var = "#696969", # Variables color
#                 col.ind = groups,
#                 geom=c("point"),
#                 addEllipses = TRUE, # Concentration ellipses
#                 ellipse.type = "norm",
#                 ellipse.level=0.90
# )
# dev.off()
