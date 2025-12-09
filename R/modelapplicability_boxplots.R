
r=read.csv("C:/Users/jenni/Box/NAMC (Trip Armstrong)/OEModeling/NAMC_Supported_OEmodels/Model Applicability/CalPredsModelApplicability.csv")
r=subset(r,modelID==136)
r$logWsArea=log10(r$WsAreaSqKm)
windows()
par(mfrow=c(2,2))
boxplot(ElevCat~Model,xlab="",ylab="",main="Catchment Elevation(m)", dat=r)
boxplot(Precip8110Ws~Model,xlab="",ylab="",main="Mean annual precipitation (mm)",dat=r)
boxplot(Tmean8110Ws~Model,xlab="",ylab="",main="Mean annual temperature (C)",dat=r)
boxplot(logWsArea~Model,xlab="",ylab="",main="Log Watershed Area (km2)",dat=r)


