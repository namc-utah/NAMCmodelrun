Crater=sf::st_read("C://Users//andrew.caudillo.BUGLAB-I9//Box//NAMC WATS Department Files//GIS//General_Layers//CraterLakeBoundary.shp")
#query AREMP to see if any sites fall into the Crater boundary:

AREMP=NAMCr::query('sites',
                   boxId=) #change this to the box you are interested in

AREMP=sf::st_as_sf(AREMP,coords=c('longitude','latitude'),crs=4269)
AREMP=sf::st_transform(AREMP,sf::st_crs(Crater))
#this will tell you which sites, if any, are in the CL boundary
Crateremp=sf::st_intersection(AREMP,Crater)




Crater_lake=NAMCr::query('sampleTaxa',
             sampleIds=c(121076,
                         145291,
                         121071,
                         145256,
                         168535,
                         121074,
                         145259,
                         168538,
                         121072,
                         121077,
                         125376,
                         145257,
                         168536,
                         145294,
                         168542,
                         121075,
                         121078,
                         125377,
                         145260,
                         168539,
                         121073,
                         145258,
                         168537,
                         168543,
                         145292,
                         168540,
                         168544,
                         104151,
                         104152,
                         104153
             ))
