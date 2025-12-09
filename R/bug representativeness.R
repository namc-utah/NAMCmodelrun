box=read.csv("C:/Users/jenni/OneDrive - USU/Documents/modelapplic.csv")

model=unique(box$Model)

for (m in 1:12){
   png(file=paste(model[m],".png"),width=1000,height=1000,pointsize=15)
   par(mfrow=c(2,2),cex=1.2,mar=c(.5,2.5, 2, 1))
  EPAsubset=subset(box,Model==model[m])
  boxplot(EPAsubset$ElevCat,main='Catchment Elevation (m)',data=EPAsubset)
  boxplot(EPAsubset$Precip8110Ws,main='Mean annual precipitation (mm)',data=EPAsubset)
  boxplot(EPAsubset$Tmean8110Ws,main='Mean annual temperature (C)',data=EPAsubset)
  boxplot(log10(EPAsubset$WsAreaSqKm),main='Log Watershed Area (km2)',data=EPAsubset)
 dev.off()
   }
