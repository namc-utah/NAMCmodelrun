#assigning cold vs warm water InvertReg to AZ sites for AZ IBI
#Andrew Caudillo of NAMC
#different thresholds are used for cold vs warm water
#in the ADEQ IBI
#this code finds elevation for the site locations
#and assings cold or warm, based on elevation

#read in the package we need
library(elevatr)
library(NAMCr)
library(sf)
#change the boxID as needed
boxID=3793
#read in samples
samps = NAMCr::query(
  "samples",
  boxId=boxID
)
#create spatial df
spatdf<-st_as_sf(samps,coords=c("siteLongitude",'siteLatitude'))
#assing lat/long projection
spatdf<-st_set_crs(spatdf,"+proj=longlat")
#get the elevation (returns a df)
invertreg_values<-get_elev_point(spatdf)
#convert meters to feet (we could keep meters and just make a new
#conditional statement, but this is fine.)
invertreg_values$elevation<-invertreg_values$elevation*3.28084; invertreg_values$elev_units<-"feet"

#assign the cold or warm
invertreg_values$Type=ifelse(invertreg_values$elevation>5000,"cold","warm")
#append sampleId (helps to find cold vs warm sites
#on INSTAR, since we need to add the values individually on EDAS)
invertreg_values$sampleId<-samps$sampleId
#get site name, too. This is what we will use when
#assigning InvertReg on EDAS, since the ActivityID is less intuitive there.
#(i.e., filter by the site names below and then add "warm" or "cold" accordingly)
invertreg_values$site<-samps$siteName
