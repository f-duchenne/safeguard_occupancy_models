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

############################ LOADING AND ASSEMBLING BEE DATA
dat=readRDS("Bee_DB.rds")
dat2 <- dat %>% select (scientificName,endYear,endMonth,endDay,decimalLongitude,decimalLatitude,country,genus,family,occurrenceID,datasetProvider)
names(dat2)<- c("TAXON","YEAR_2","MONTH_2","DAY_2","LONGITUDE","LATITUDE","COUNTRY","GENUS","FAMILY","UUID","DATABASE_REFERENCE_CODE_2")

## removing data from spain
dat2 <- dat2[!(dat2$DATABASE_REFERENCE_CODE_2 %in% c("Ignasi Bartomeus","Nacho Bartomeus" ,"I. Bartomeus")),]

# Loading more recent data from spain:
data_Spain <- fread("iberian_bees.csv")
taxi=fread("species_family_table,.csv")
nrow(data_Spain)
data_Spain=merge(data_Spain,taxi,by.x="Accepted_name",by.y="TAXON")
nrow(data_Spain)
## I include approximate coordinates for two villages that do not have but hold 2750 records together (Pina de Ebro and Castello de Vide).
data_Spain$Longitude[data_Spain$Locality == "Castello de Vide"] <- -7.45680
data_Spain$Latitude[data_Spain$Locality == "Castello de Vide"] <- 39.41624

## Assign coordinates for "Pina de Ebro"
data_Spain$Longitude[data_Spain$Locality == "Pina de Ebro"] <- -0.5261
data_Spain$Latitude[data_Spain$Locality == "Pina de Ebro"] <- 41.4814 

## Select the columns that are common between both datasets and give them the right names
Spain <- data_Spain %>% select (Accepted_name, Year, Month, Day, Longitude,Latitude,Country,Genus,FAMILY,Unique.identifier)
Spain$DATABASE_REFERENCE_CODE_2="BartomeusI_Iberia_2023"
names(Spain)<- c("TAXON","YEAR_2","MONTH_2","DAY_2","LONGITUDE","LATITUDE","COUNTRY","GENUS","FAMILY","UUID","DATABASE_REFERENCE_CODE_2")

# Adding data from Portugal (https://doi.org/10.15468/dl.7vqj99):
data_port=fread("0000019-250214093936778.csv")
data_port$Country="Portugal"

## Select the columns that are common between both datasets and give them the right names
Portugal <- data_port %>% select (species, year, month, day, decimalLongitude,decimalLatitude,Country,genus,family,occurrenceID,datasetKey)
names(Portugal)<- c("TAXON","YEAR_2","MONTH_2","DAY_2","LONGITUDE","LATITUDE","COUNTRY","GENUS","FAMILY","UUID","DATABASE_REFERENCE_CODE_2")

## Integrate Iberian database with the European dataset.
data1 <- rbind(dat2, Spain,Portugal)
data1=subset(data1,TAXON!="Apis mellifera")

############################ LOADING AND ASSEMBLING HOVERFLIES DATA
dath=readRDS("Hoverfly_DB.rds")
dath <- dath %>% select (scientificName,startYear,startMonth,startDay,decimalLongitude,decimalLatitude,country,genus,family,occurrenceID,datasetProvider)
names(dath)<- c("TAXON","YEAR_2","MONTH_2","DAY_2","LONGITUDE","LATITUDE","COUNTRY","GENUS","FAMILY","UUID","DATABASE_REFERENCE_CODE_2")
dath$taxo_group="hoverflies"
data1$taxo_group="bees"
## Integrate hoverflies with bees
data2 <- rbind(data1, dath)


############################ FILTERS THE RECORDS TO KEEP ONLY THE ONE WITH YEAR AND COORDINATES
data2_subset=subset(data2, !is.na(YEAR_2) & !is.na(LONGITUDE) & !is.na(LATITUDE))

filter1=(nrow(data2)-nrow(data2_subset)) #number of records removed

# LOAD SOME BASIC SHAPEFILES THAT WILL USE TO DEFINE THE AREA
# Define the UTM projection for a suitable UTM zone
utm_crs <- st_crs("+proj=utm +zone=32 +ellps=WGS84")

# Define bounding box coordinates for Spain to northern Finland in the chosen UTM zone
bbox <- st_bbox(c(xmin = -1200000, xmax = 3500000, ymin = 3500000, ymax = 8000000), crs = utm_crs)
bboxsf=st_as_sf(st_as_sfc(st_bbox(bbox)))

# FILTERS THE RECORDS TO KEEP ONLY THE ONE IN THE DEFINED AREA
# Set the CRS of data3 to be the same as that of the map
data3 <- st_as_sf(data2_subset,coords = c("LONGITUDE", "LATITUDE"), crs = 4326) %>% st_transform(crs = utm_crs)
  
#crop to the extent of the map
data3_filtered <- st_crop(data3, bbox)

filter2=nrow(data3)-nrow(data3_filtered)

######################### MAKE DIFFERENTE GRIDS WITH DIFFERENT RESOLUTION
bioregions <- bioregions <- st_read ("D:/land use change/biogeographic_regions/BiogeoRegions2016.shp")
bioregions <- st_transform(bioregions, crs = utm_crs)
resolutions=c(20000,50000,100000,200000)
plot(st_geometry(bioregions))
plot(st_geometry(st_as_sfc(bbox)),add=TRUE)

# for(i in resolutions){
	# # Create hexagonal grid cells
	# hex_grid <- hexagons(bboxsf, res = i)
	# hex_grid[,paste0("gridID")]=1:nrow(hex_grid)
	# hex_grid <- st_transform(hex_grid, crs = utm_crs)
	# hex_grid=st_join(hex_grid,bioregions,largest=TRUE)
	# names(hex_grid)=c("gridID","PK_UID","region","p2012","code","name","geometry")
	# ge=which(names(hex_grid)=="geometry")
	# names(hex_grid)[-ge]=paste0(names(hex_grid)[-ge],"_",i/1000)
	# st_write(hex_grid,paste0("grid_",i/1000,"KM.shp")) #export the grid to be able to access them latter and furnish them with the paper
# }

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


coords=as.data.table(st_coordinates(data4))
names(coords)=c("LONGITUDE","LATITUDE")
data4_filtered=cbind(as.data.table(data4),coords)

data4_filtered=subset(data4_filtered,!is.na(region_20) & region_20!="outside" & !is.na(region_50) & !is.na(region_100) & !is.na(region_200))

filter3=nrow(data4)-nrow(data4_filtered)

############ ATTRIBUTE TIME PERIOD
## Time periods:
#P1 <- 1921 - 1940 # Every 20 years.
#P2 <- 1941 - 1960
#P3 <- 1961 - 1980
#P4 <- 1981 - 2000
#P5 <- 2001 - 2020
#define time periods
data4_filtered$time_period=NA
data4_filtered$time_period[data4_filtered$YEAR_2>=1921 & data4_filtered$YEAR_2<=1940]= "1921-1940"
data4_filtered$time_period[data4_filtered$YEAR_2>=1941 & data4_filtered$YEAR_2<=1960]="1941-1960"
data4_filtered$time_period[data4_filtered$YEAR_2>=1961 & data4_filtered$YEAR_2<=1980]="1961-1980"
data4_filtered$time_period[data4_filtered$YEAR_2>=1981 & data4_filtered$YEAR_2<=2000]="1981-2000"
data4_filtered$time_period[data4_filtered$YEAR_2>=2001 & data4_filtered$YEAR_2<=2020]="2001-2020"

dataf=subset(data4_filtered,!is.na(time_period) & !is.na(region_50))

filter4=nrow(data4_filtered)-nrow(dataf)

fwrite(dataf,"database_clean_filtered.csv")
filters=data.frame(filters=c("year and coordinates","geographical extent","not in a region","not in a period"),nb_removed=c(filter1,filter2,filter3,filter4))
fwrite(filters,"nb_records_removed_during _filtering.csv")



library(rgbif)
tab=data.frame(species=unique(dat$TAXON),family=NA)
for(i in 1:nrow(tab)){
  obj=name_backbone(name=gsub("[^0-9A-Za-z///' ]","",tab$species[i],ignore.case=T),order ="Hymenoptera") #look for the given species
  if(length(obj$rank)>0 & !is.na(tab$species[i])){  #if obj is not empty
    if(obj$status=="SYNONYM"){obj=name_backbone(name=name_usage(key=obj$acceptedUsageKey)$data$scientificName)}
	tab$canon[i]=ifelse(length(obj$canonicalName)>0,obj$canonicalName,NA) 
    tab$family[i]=ifelse(length(obj$family)>0,obj$family,NA) #family name according to the GBIF
  }
}
vec=tab$species[duplicated(tab$canon)]
dat$COUNTRY[dat$TAXON %in% vec] 

taxi=rbind(taxi,tab)