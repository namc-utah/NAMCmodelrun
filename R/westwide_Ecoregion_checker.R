#desert region checker for Westwide
site_dat=NAMCr::query('sites',boxIds=10822)
ecos=sf::st_read("C://Users//andrew.caudillo.BUGLAB-I9//Box//NAMC WATS Department Files//GIS//GIS_Stats//CONUS//ecoregion//hybrid_level_III_ecoregions - Copy.shp")

site_dat=sf::st_as_sf(site_dat,coords=c('longitude','latitude'),crs=4269)
site_dat=sf::st_transform(site_dat,sf::st_crs(ecos))

site_regions<-sf::st_intersection(site_dat,ecos)
site_regions<-site_regions$EPA_hybird
table(site_regions)
