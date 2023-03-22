samples=NAMCr::query("samples",projectId=49)
sites=NAMCr::query("sites",sampleIds=samples$sampleId)
d=dplyr::left_join(samples,sites, by="siteId")
siteIds=unlist(unique(sites$siteId))
#read in only watersheds from the mastersheds layer that are in the current box or project
existingWatersheds=sf::st_make_valid(sf::st_read(watershed_file_path, query=sprintf('SELECT * FROM %s WHERE siteId in(%s)',watershed_layer_name, inLOOP(substr(siteIds, 1, 10)))))
pred_geometry_base_path="C:/Users/jenni/Box/NAMC (Trip Armstrong)/"
SQLite_file_path="C:/NAMC_S3/StreamCat/StreamCat2022.sqlite"
watershed_file_path=paste0(pred_geometry_base_path,"GIS/Watersheds/Mastersheds/mastersheds.shp") #siteId must be in file!!!!
watershed_layer_name="mastersheds" #siteId must be in file!!!!
 inLOOP<- function(inSTR,...) {
inSTR=unlist(inSTR)
if (inSTR[1]==''){loopSTR="''"} else{
 for (i in 1:length(inSTR)){
comma=ifelse(i==length(inSTR),'',',')
 STRl=sprintf("'%s'%s",inSTR[i],comma)
if(i==1){loopSTR=STRl} else{loopSTR=paste(loopSTR,STRl)}
} }
return(loopSTR)
}


existingWatersheds=sf::st_make_valid(sf::st_read(watershed_file_path, query=sprintf('SELECT * FROM %s WHERE siteId in(%s)',watershed_layer_name, inLOOP(substr(siteIds, 1, 10)))))
points_subset=subset(sites,!(siteId %in% existingWatersheds$siteId))
write.csv(d,"d.csv")
gradients=read.csv('gradients.csv')

existingWatersheds$siteId=as.numeric(existingWatersheds$siteId)
existingWatersheds2=existingWatersheds[!duplicated(existingWatersheds$siteId),]
gradients=read.csv('gradients.csv')
total=merge(existingWatersheds2,sites,by='siteId')

total2=merge(total,gradients,by='waterbodyCode')

total2=merge(total,gradients,by='waterbodyCode',all=TRUE)
total3=dplyr::inner_join(total2,samples,by='siteId')
total2$COMID=total2$waterbodyCode
total3=total2[,c('COMID','siteId','ElevCat','WsAreaSqKm','Precip8110Ws','Tmean8110Ws')]
total3$COMID=as.numeric(total3$COMID)
sf::st_write(total3,'watersheds.shp')
samples2=samples[!duplicated(samples$siteId),]
samples2$EvaluationID=samples2$visitId
samples2$PointID=samples2$customerSiteCode
samples2=samples2[,c('siteId','EvaluationID','PointID')]
total4=merge(total3,samples2,by='siteId')
total5=subset(total4,is.na(EvaluationID)==FALSE)
sf::st_write(total5,'watersheds5.shp')
