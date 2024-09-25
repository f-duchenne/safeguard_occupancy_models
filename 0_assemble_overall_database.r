###########################################
#' Check for packages and if necessary install into library 
#+ message = FALSE
rm(list=ls())
pkgs <- c("data.table", "dplyr","sf","spatialEco") 

inst <- pkgs %in% installed.packages()
if (any(inst)) install.packages(pkgs[!inst])
pkg_out <- lapply(pkgs, require, character.only = TRUE)

#defining working folder:
setwd(dir="C:/Users/Duchenne/Documents/safeguard/data")

# Loading SAFEGUARD data
dat=readRDS("Safeguard_bee_df_final_v9_20240902.rds")

#removing data from spain
dat <- dat [ !dat$DATABASE_REFERENCE_CODE_2 == "BartomeusI_Iberia_2023",]

#Loading more recent data from spain:
data_Spain <- fread("iberian_bees.csv")

# I include approximate coordinates for two villages that do not have but hold 2750 records together (Pina de Ebro and Castello de Vide).
data_Spain$Longitude[data_Spain$Locality == "Castello de Vide"] <- -7.45680
data_Spain$Latitude[data_Spain$Locality == "Castello de Vide"] <- 39.41624

# Assign coordinates for "Pina de Ebro"
data_Spain$Longitude[data_Spain$Locality == "Pina de Ebro"] <- -0.5261
data_Spain$Latitude[data_Spain$Locality == "Pina de Ebro"] <- 41.4814 

######################### FILTERS THE RECORDS TO KEEP ONLY THE ONE WITH YEAR AND COORDINATES
data_Spain_subset=subset(data_Spain, !is.na(Year) & !is.na(Longitude) & !is.na(Latitude))

dat2=subset(dat, !is.na(YEAR_2) & !is.na(LONGITUDE) & !is.na(LATITUDE))

filter1=(nrow(data_Spain)-nrow(data_Spain_subset))+(nrow(dat)-nrow(dat2)) #number of records removed

#Select the columns that are common between both datasets and give them the right names
Spain <- data_Spain_subset %>% select (Accepted_name, Year, Month, Day, Longitude,Latitude,Country,Genus,Unique.identifier)
names(Spain)<- c("TAXON","YEAR_2","MONTH_2","DAY_2","LONGITUDE","LATITUDE","COUNTRY","GENUS","UUID")

dat2 <- dat2 %>% select (TAXON,YEAR_2,MONTH_2,DAY_2,LONGITUDE,LATITUDE,COUNTRY,GENUS,UUID)

## Integrate Iberian database with the European dataset.
data2 <- rbind(dat2, Spain)

############### LOAD SOME BASIC SHAPEFILES THAT WILL USE TO DEFINE THE AREA
# Define the UTM projection for a suitable UTM zone
utm_crs <- st_crs("+proj=utm +zone=32 +ellps=WGS84")

# Define bounding box coordinates for Spain to northern Finland in the chosen UTM zone
bbox <- st_bbox(c(xmin = -1200000, xmax = 3500000, ymin = 3500000, ymax = 8000000), crs = utm_crs)
bboxsf=st_as_sf(st_as_sfc(st_bbox(bbox)))

######################### FILTERS THE RECORDS TO KEEP ONLY THE ONE IN THE DEFINED AREA
# Set the CRS of data3 to be the same as that of the map
data3 <- st_as_sf(data2,coords = c("LONGITUDE", "LATITUDE"), crs = 4326) %>% st_transform(crs = utm_crs)
  
#crop to the extent of the map
data3_filtered <- st_crop(data3, bbox)

filter2=nrow(data3)-nrow(data3_filtered)

######################### MAKE DIFFERENTE GRIDS WITH DIFFERENT RESOLUTION
bioregions <- bioregions <- st_read ("D:/land use change/biogeographic_regions/BiogeoRegions2016.shp")
bioregions <- st_transform(bioregions, crs = utm_crs)
resolutions=c(20000,50000,100000,200000)

for(i in resolutions){
	# Create hexagonal grid cells
	hex_grid <- hexagons(bboxsf, res = i)
	hex_grid[,paste0("gridID")]=1:nrow(hex_grid)
	hex_grid <- st_transform(hex_grid, crs = utm_crs)
	hex_grid=st_join(hex_grid,bioregions,largest=TRUE)
	names(hex_grid)=c("gridID","PK_UID","region","p2012","code","name","geometry")
	ge=which(names(hex_grid)=="geometry")
	names(hex_grid)[-ge]=paste0(names(hex_grid)[-ge],"_",i/1000)
	st_write(hex_grid,paste0("grid_",i/1000,"KM.shp")) #export the grid to be able to access them latter and furnish them with the paper
}

######################### MERGE DATA AND GRIDS
nbr=nrow(data3_filtered)
for(i in resolutions){
	hex_grid=st_read(paste0("grid_",i/1000,"KM.shp"),crs=utm_crs)
	data3_filtered=st_join(data3_filtered,hex_grid[,paste0(c("gridID","region"),"_",i/1000)])
	data3_filtered=data3_filtered[,-which(names(data3_filtered)=="geometry")]
	print(nrow(data3_filtered))
}


#SOME RECORDS FALL ON THE EDGE OF TWO CELLS AND ARE THUS DUPLICATED; ATTRIBUTE THEM TO ONE ARBITRARY
bidon=data3_filtered[data3_filtered$UUID %in% data3_filtered$UUID[duplicated(data3_filtered$UUID)],]
unique(bidon$region_20) #check that all concern cells have a region
unique(bidon$region_50) #check that all concern cells have a region
unique(bidon$region_100) #check that all concern cells have a region
unique(bidon$region_200) #check that all concern cells have a region

data4=data3_filtered[!duplicated(data3_filtered$UUID),] #remove the duplicated values

############ SOME GRIDCELLS ALONG THE BORDER DO NOT HAVE A REGION; REATRRIBUTE THEM
data4$region_20[is.na(data4$region_20) & !is.na(data4$region_50)]=data4$region_50[is.na(data4$region_20) & !is.na(data4$region_50)]

data4_filtered=subset(data4,!is.na(region_20) & region_20!="outside" & !is.na(region_50) & !is.na(region_100) & !is.na(region_200))

filter3=nrow(data4)-nrow(data4_filtered)

fwrite(data4_filtered,"database_clean_filtered.csv")
filters=data.frame(filters=c("year and coordinates","geographical extent","not in a region"),nb_removed=c(filter1,filter2,filter3))
fwrite(filters,"nb_records_removed_during _filtering.csv")



