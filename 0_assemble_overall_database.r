###########################################
#' Check for packages and if necessary install into library 
#+ message = FALSE
rm(list=ls())
pkgs <- c("data.table", "dplyr","sf","spatialEco") 

inst <- pkgs %in% installed.packages()
if (any(inst)) install.packages(pkgs[!inst])
pkg_out <- lapply(pkgs, require, character.only = TRUE)

#defining working folder:
project_folder="C:/Users/Duchenne/Documents/safeguard/"

############################ LOADING AND ASSEMBLING BEE DATA
datb=readRDS(paste0(project_folder,"data/Bee_DB_2025-09-03.rds"))

data1b <- datb %>% select (scientificName,endYear,endMonth,endDay,decimalLongitude,decimalLatitude,country,genus,family,occurrenceID,datasetProvider,institutionName,isPseudodata)
data1b$taxo_group="bees"

ndata=nrow(data1b)

#remove Apis mellifera
data1b=subset(data1b,scientificName!="Apis mellifera")

#### remove duplicated and blurr records in 1989

data1b=subset(data1b,isPseudodata==FALSE)
filter1=ndata-nrow(data1b)

############################ LOADING AND ASSEMBLING HOVERFLIES DATA
dath=readRDS(paste0(project_folder,"data/Hoverfly_DB_2025-09-03.rds"))
data1h <- dath %>% select (scientificName,endYear,endMonth,endDay,decimalLongitude,decimalLatitude,country,genus,family,occurrenceID,datasetProvider,institutionName)
data1h$taxo_group="hoverflies"

## Integrate hoverflies with bees
data2 <- rbind(data1b[,-which(names(data1b)=="isPseudodata")], data1h)

############################ FILTERS THE RECORDS TO KEEP ONLY THE ONE WITH YEAR AND COORDINATES
data2_subset=subset(data2, !is.na(endYear) & !is.na(decimalLongitude) & !is.na(decimalLatitude))

filter2=(nrow(data2)-nrow(data2_subset)) #number of records removed

# LOAD SOME BASIC SHAPEFILES THAT WILL USE TO DEFINE THE AREA
# Define the UTM projection for a suitable UTM zone
utm_crs <- st_crs("+proj=utm +zone=32 +ellps=WGS84")

# Define bounding box coordinates for Spain to northern Finland in the chosen UTM zone
bbox <- st_bbox(c(xmin = -1200000, xmax = 3500000, ymin = 3500000, ymax = 8000000), crs = utm_crs)
bboxsf=st_as_sf(st_as_sfc(st_bbox(bbox)))

# FILTERS THE RECORDS TO KEEP ONLY THE ONE IN THE DEFINED AREA
# Set the CRS of data3 to be the same as that of the map
data3 <- st_as_sf(data2_subset,coords = c("decimalLongitude", "decimalLatitude"), crs = 4326) %>% st_transform(crs = utm_crs)
  
#crop to the extent of the map
data3_filtered <- st_crop(data3, bbox)

filter3=nrow(data3)-nrow(data3_filtered)


######################### MAKE DIFFERENTE GRIDS WITH DIFFERENT RESOLUTION
bioregions <- st_read ("D:/land use change/biogeographic_regions/BiogeoRegions2016.shp")
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
	# st_write(hex_grid,paste0(project_folder,"data/grid_",i/1000,"KM.shp")) #export the grid to be able to access them latter and furnish them with the paper
# }

######################### MERGE DATA AND GRIDS
nbr=nrow(data3_filtered)
for(i in resolutions){
	hex_grid=st_read(paste0(project_folder,"data/grid_",i/1000,"KM.shp"),crs=utm_crs)
	data3_filtered=st_join(data3_filtered,hex_grid[,paste0(c("gridID","region"),"_",i/1000)])
	data3_filtered=data3_filtered[,-which(names(data3_filtered)=="geometry")]
	print(nrow(data3_filtered))
}


#SOME RECORDS FALL ON THE EDGE OF TWO CELLS AND ARE THUS DUPLICATED; ATTRIBUTE THEM TO ONE ARBITRARY
bidon=data3_filtered[data3_filtered$occurrenceID %in% data3_filtered$occurrenceID[duplicated(data3_filtered$occurrenceID)],]
unique(bidon$region_20) #check that all concern cells have a region
unique(bidon$region_50) #check that all concern cells have a region
unique(bidon$region_100) #check that all concern cells have a region
unique(bidon$region_200) #check that all concern cells have a region

data4=data3_filtered[!duplicated(data3_filtered$occurrenceID),] #remove the duplicated values

############ SOME GRIDCELLS ALONG THE BORDER DO NOT HAVE A REGION; REATRRIBUTE THEM
data4$region_20[is.na(data4$region_20) & !is.na(data4$region_50)]=data4$region_50[is.na(data4$region_20) & !is.na(data4$region_50)]


coords=as.data.table(st_coordinates(data4))
names(coords)=c("decimalLongitude", "decimalLatitude")
data4_filtered=cbind(as.data.table(data4),coords)

data4_filtered=subset(data4_filtered,!is.na(region_20) & region_20!="outside" & !is.na(region_50) & !is.na(region_100) & !is.na(region_200))

filter4=nrow(data4)-nrow(data4_filtered)

############ ATTRIBUTE TIME PERIOD
## Time periods:
#P1 <- 1921 - 1940 # Every 20 years.
#P2 <- 1941 - 1960
#P3 <- 1961 - 1980
#P4 <- 1981 - 2000
#P5 <- 2001 - 2020
#define time periods
data4_filtered$time_period=NA
data4_filtered$time_period[data4_filtered$endYear>=1921 & data4_filtered$endYear<=1940]= "1921-1940"
data4_filtered$time_period[data4_filtered$endYear>=1941 & data4_filtered$endYear<=1960]="1941-1960"
data4_filtered$time_period[data4_filtered$endYear>=1961 & data4_filtered$endYear<=1980]="1961-1980"
data4_filtered$time_period[data4_filtered$endYear>=1981 & data4_filtered$endYear<=2000]="1981-2000"
data4_filtered$time_period[data4_filtered$endYear>=2001 & data4_filtered$endYear<=2020]="2001-2020"

dataf=subset(data4_filtered,!is.na(time_period) & !is.na(region_50))

filter5=nrow(data4_filtered)-nrow(dataf)

fwrite(dataf,paste0(project_folder,"data/database_clean_filtered.csv"))
filters=data.frame(filters=c("year and coordinates","geographical extent","duplicated","not in a region","not in a period"),nb_removed=c(filter1,filter2,filter3,filter4,filter5))
fwrite(filters,paste0(project_folder,"data/nb_records_removed_during _filtering.csv"))



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