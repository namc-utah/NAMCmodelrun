#Calculate MMI for test data:


setwd("\\\\share1.bluezone.usu.edu\\miller\\buglab\\OE_Modeling\\NAMC_Supported_OEmodels\\NV\\NV_MMI\\Model Files")

## General model run

load("OE_MMI_models.rdata")
library(randomForest)
library(vegan)
Boxdata <- ("\\\\share1.bluezone.usu.edu\\miller\\\\buglab\\OE_Modeling\\NAMC_Supported_OEmodels\\NV\\NV_MMI\\InputsAndResults\\CurrentRun")

setwd(Boxdata)
getwd()
prednew <- read.csv("NV_Habitat.csv",row.names="SampleID")
new.metrics <- read.csv("NV_Bugs.csv",row.names="SampleID")


#modeling, adjusting, rescaling
predictors=prednew
new.metrics=new.metrics[row.names(predictors),]
namecheck=row.names(predictors)==row.names(new.metrics)
if (any(namecheck=="FALSE")) stop("Predictors do not match bug metrics")

####adjust metrics for natural variability and rescale
INSET.raw=new.metrics$INSET
INSET.pred=predict(INSET.rf, newdata=predictors, type="response")
INSET.adj=INSET.raw-INSET.pred
INSET.rs=100*((INSET.adj--20.23451)/(9.53285--20.23451))
INSET.rs[INSET.rs>100]=100
INSET.rs[INSET.rs<0]=0

PER_CFA.raw=new.metrics$PER_CFA
PER_CFA.pred=predict(PER_CFA.rf, newdata=predictors, type="response")
PER_CFA.adj=PER_CFA.raw-PER_CFA.pred
PER_CFA.rs=100*((PER_CFA.adj--26.84895)/(20.15404--26.84895))
PER_CFA.rs[PER_CFA.rs>100]=100
PER_CFA.rs[PER_CFA.rs<0]=0

PER_EPHEA.raw=new.metrics$PER_EPHEA
PER_EPHEA.pred=predict(PER_EPHEA.rf, newdata=predictors, type="response")
PER_EPHEA.adj=PER_EPHEA.raw-PER_EPHEA.pred
PER_EPHEA.rs=100*((PER_EPHEA.adj--37.05147)/(29.21069--37.05147))
PER_EPHEA.rs[PER_EPHEA.rs>100]=100
PER_EPHEA.rs[PER_EPHEA.rs<0]=0

NONSET.raw=new.metrics$NONSET
NONSET.pred=predict(NONSET.rf, newdata=predictors, type="response")
NONSET.adj=NONSET.raw-NONSET.pred
NONSET.rs=100*((NONSET.adj--3.733934)/(3.773724--3.733934))
NONSET.rs[NONSET.rs>100]=100
NONSET.rs[NONSET.rs<0]=0

CLINGER.raw=new.metrics$CLINGER
CLINGER.pred=predict(CLINGER.rf, newdata=predictors, type="response")
CLINGER.adj=CLINGER.raw-CLINGER.pred
CLINGER.rs=100*((CLINGER.adj--8.310543)/(5.719022--8.310543))
CLINGER.rs[CLINGER.rs>100]=100
CLINGER.rs[CLINGER.rs<0]=0

PER_PLECA.raw=new.metrics$PER_PLECA
PER_PLECA.pred=predict(PER_PLECA.rf, newdata=predictors, type="response")
PER_PLECA.adj=PER_PLECA.raw-PER_PLECA.pred
PER_PLECA.rs=100*((PER_PLECA.adj--6.212933)/(15.33518--6.212933))
PER_PLECA.rs[PER_PLECA.rs>100]=100
PER_PLECA.rs[PER_PLECA.rs<0]=0

SHDIVER.raw=new.metrics$SHDIVER #SHDIVER not adjusted by natural gradients
SHDIVER.rs=100*((SHDIVER.raw-0.8980395)/(3.18087-0.8980395))
SHDIVER.rs[SHDIVER.rs>100]=100
SHDIVER.rs[SHDIVER.rs<0]=0

#Calculate MMI scores for new samples
new.metrics.rs=data.frame(INSET.rs, PER_EPHEA.rs, SHDIVER.rs, PER_CFA.rs, PER_PLECA.rs, NONSET.rs, CLINGER.rs, row.names=row.names(new.metrics))
MMI.scores=rowSums(new.metrics.rs)/7


#Nevada ref site score mean MMI = 55.98970963
MMI_Condition=vector()
for(n in 1:length(MMI.scores)){
  if (MMI.scores[n]<=44.58171) Condition="Impaired"
  if (MMI.scores[n]<47.01497 & MMI.scores[n]>44.58171) Condition="Undetermined"
  if (MMI.scores[n]>=47.01497) Condition="Reference"
  MMI_Condition=append(MMI_Condition, Condition)}

FinalResults=cbind(MMI.scores,MMI_Condition)

write.csv(FinalResults, file = "NV_MMI_Scores.csv")


###This is the model applicability section
###Do not run below this until the code is sorted out

############################
#Environmental outlier test#
############################
#Uses watershed area (km2), ELVmean_WS (m), Pmax_WS (mm)


#Calculating distances between reference sites
logSQ_KM=log(pred.cal$SQ_KM)
ELVmean_WS=pred.cal$ELVmean_WS
logPmax_WS=log(pred.cal$Pmax_WS)
outlier.preds.ref.raw=data.frame(logSQ_KM,ELVmean_WS,logPmax_WS, row.names=row.names(pred.cal))
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


logSQ_KM=log(prednew$SQ_KM)
ELVmean_WS=prednew$ELVmean_WS
logPmax_WS=log(prednew$Pmax_WS)
outlier.preds.test.raw=data.frame(logSQ_KM,ELVmean_WS,logPmax_WS, row.names=row.names(prednew))
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
test.dist=test.dist[row.names(outlier.preds.test.std),row.names(pred.cal)]
test.top10mean.dist=vector()
out.flag90=vector()
sample.list = row.names(predictors)
for (n in 1:length(sample.list)){
  if(length(sample.list)>1){cur.dist=test.dist[n,]}
  if(length(sample.list)==1){cur.dist=test.dist}
  cur.10dist=(cur.dist[order(cur.dist)])[1:10]
  mean10dist=mean(cur.10dist)
  test.top10mean.dist=append(test.top10mean.dist, mean10dist)
  if (mean10dist>quantile(ref.top10mean.dist,0.90)) flag90="Yes" else flag90="No"
  out.flag90=append(out.flag90, flag90)
}



###THIS IS WRONG!!! THE FLAG RESULTS DO NOT MATCH THE INDEX RESULTS DUE TO THE WAY THINGS ARE SORTED.Do not run these together, run seperately or order the same as the bugs & predictors are ordered and fixed in the OE index.  
final=cbind(FinalResults,out.flag90)
##write.csv(final, file="Corrected_FinalResults.csv")

write.csv(final, file="Corrected_FinalResults.csv")

