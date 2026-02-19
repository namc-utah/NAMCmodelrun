#desert region checker for Westwide
site_dat=NAMCr::query('sites',boxIds=10822) #get the site info
#read in the ecoregions shapefile
ecos=sf::st_read("C://Users//andrew.caudillo.BUGLAB-I9//Box//NAMC WATS Department Files//GIS//GIS_Stats//CONUS//ecoregion//hybrid_level_III_ecoregions - Copy.shp")
#convert the sites to spatial objects and convert the coordinates to the ecos object's crs
site_dat=sf::st_as_sf(site_dat,coords=c('longitude','latitude'),crs=4269)
site_dat=sf::st_transform(site_dat,sf::st_crs(ecos))

#intersect the two datasets
site_regions<-sf::st_intersection(site_dat,ecos)
site_regions<-site_regions$EPA_hybird
#look to see if the xeric ecoregion appears
#if the xeric ecoregion appears, then modelId 26 has to be used,
#perhaps along with 25, based on the total number of ecoregions.
table(site_regions)
