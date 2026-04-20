###########################################
#' Check for packages and if necessary install into library 
#+ message = FALSE
rm(list=ls())
pkgs <- c("data.table", "dplyr","sf","spatialEco","ggplot2") 

inst <- pkgs %in% installed.packages()
if (any(inst)) install.packages(pkgs[!inst])
pkg_out <- lapply(pkgs, require, character.only = TRUE)

#defining working folder:
project_folder="C:/Users/Duchenne/Documents/safeguard/"

# LOAD SOME BASIC SHAPEFILES THAT WILL USE TO DEFINE THE AREA
# Define the UTM projection for a suitable UTM zone
utm_crs <- st_crs("+proj=utm +zone=32 +ellps=WGS84")

# Define bounding box coordinates for Spain to northern Finland in the chosen UTM zone
bbox <- st_bbox(c(xmin = -1200000, xmax = 3500000, ymin = 3500000, ymax = 8000000), crs = utm_crs)
bboxsf=st_as_sf(st_as_sfc(st_bbox(bbox)))

bioregions <- st_read ("D:/land use change/biogeographic_regions/BiogeoRegions2016.shp")
bioregions <- st_transform(bioregions, crs = utm_crs)

bioregions2 <- st_crop(bioregions, bbox)

continents <- st_read ("D:/land use change/continents/continents-of-the-world-merged.shp")

bioregions2 <- st_transform(bioregions2, crs = st_crs(continents))
bioregions2=subset(bioregions2,short_name!="outside")

ggplot()+
theme_bw()+
geom_sf(data=continents,fill="grey",color="grey")+
geom_sf(data=bioregions2,fill="black",color="black")+
coord_sf(crs = "+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +units=m +no_defs ")


