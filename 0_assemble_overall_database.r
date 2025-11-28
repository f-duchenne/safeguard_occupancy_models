###########################################
#' Check for packages and if necessary install into library 
#+ message = FALSE
rm(list=ls())
pkgs <- c("data.table", "dplyr","sf","spatialEco") 

inst <- pkgs %in% installed.packages()
if (any(inst)) install.packages(pkgs[!inst])
pkg_out <- lapply(pkgs, require, character.only = TRUE)

#defining working folder if needed:
#project_folder="C:/Users/Duchenne/Documents/safeguard/"
project_folder=""

############################ LOADING AND ASSEMBLING BEE DATA
datb=readRDS(paste0(project_folder,"data/raw_data/Bee_DB_2025-09-03.rds"))

data1b <- datb %>% select (scientificName,endYear,endMonth,endDay,decimalLongitude,decimalLatitude,country,genus,family,occurrenceID,datasetProvider,institutionName,isPseudodata,scientificNameAuthorship)
data1b$taxo_group="bees"

ndata=nrow(data1b)

#remove Apis mellifera
data1b=subset(data1b,scientificName!="Apis mellifera")

#### remove duplicated and blurr records in 1989

data1b=subset(data1b,isPseudodata==FALSE)
filter1=ndata-nrow(data1b)

############################ LOADING AND ASSEMBLING HOVERFLIES DATA
dath=readRDS(paste0(project_folder,"data/raw_data/Hoverfly_DB_2025-09-03.rds"))
data1h <- dath %>% select (scientificName,endYear,endMonth,endDay,decimalLongitude,decimalLatitude,country,genus,family,occurrenceID,datasetProvider,institutionName,scientificNameAuthorship)
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
bioregions <- st_read ("data/raw_data/grids_shapefiles/BiogeoRegions2016.shp")
bioregions <- st_transform(bioregions, crs = utm_crs)
resolutions=c(20000,50000,100000,200000)
plot(st_geometry(bioregions))
plot(st_geometry(st_as_sfc(bbox)),add=TRUE)

#(GRIDS ARE ALREADY DONE AND AVAILABLE IN THE DATA FOLDER, SO THE LOOP IS COMMENTED)
# for(i in resolutions){
	# # Create hexagonal grid cells
	# hex_grid <- hexagons(bboxsf, res = i)
	# hex_grid[,paste0("gridID")]=1:nrow(hex_grid)
	# hex_grid <- st_transform(hex_grid, crs = utm_crs)
	# hex_grid=st_join(hex_grid,bioregions,largest=TRUE)
	# names(hex_grid)=c("gridID","PK_UID","region","p2012","code","name","geometry")
	# ge=which(names(hex_grid)=="geometry")
	# names(hex_grid)[-ge]=paste0(names(hex_grid)[-ge],"_",i/1000)
	# st_write(hex_grid,paste0(project_folder,"data/raw_data/grids_shapefiles/grid_",i/1000,"KM.shp")) #export the grid to be able to access them latter and furnish them with the paper
# }

######################### MERGE DATA AND GRIDS
nbr=nrow(data3_filtered)
for(i in resolutions){
	hex_grid=st_read(paste0(project_folder,"data/raw_data/grids_shapefiles/grid_",i/1000,"KM.shp"),crs=utm_crs)
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

fwrite(dataf,paste0(project_folder,"data/final_and_intermediate_outputs/database_clean_filtered.csv"))
filters=data.frame(filters=c("year and coordinates","geographical extent","duplicated","not in a region","not in a period"),nb_removed=c(filter1,filter2,filter3,filter4,filter5))
fwrite(filters,paste0(project_folder,"data/final_and_intermediate_outputs/nb_records_removed_during _filtering.csv"))

########################### ASSEMBLING TRAITS

## Assemble Traits data
#this script explores the database for descriptive analysis.
pkgs <- c("data.table", "dplyr","tidyverse","rgbif") 
inst <- pkgs %in% installed.packages()
if (any(inst)) install.packages(pkgs[!inst])
pkg_out <- lapply(pkgs, require, character.only = TRUE)


#defining working folder:
#project_folder="C:/Users/Duchenne/Documents/safeguard/"
project_folder <- ""

select <- dplyr::select

Traits <- fread(paste0(project_folder,"data/raw_data/hoverfly_bee_traits_2025_04_01.csv"))
#correct a small mistake of names:
Traits$Species[Traits$Species=="Syrphus niditifrons"]="Syrphus nitidifrons"

#select useful columns
Traits2 <- Traits %>% select(Species,Order,Family,Genus, Sociality, STI_Species_temperature_index, Larval_diet_breadth, ITD_F_mm, 
								Adult_Body_size, Larval_nutrition,Flight_ability,Flight_height)
								
# ADDING DESCRIPTION DATE FOR EACH SPECIES AND KEEP ONLY SPECIES THAT WE HAVE IN OUR DATASET
datf=fread(paste0(project_folder,"data/final_and_intermediate_outputs/database_clean_filtered.csv"))
spec_li=datf[,c("scientificNameAuthorship","scientificName")]
spec_li=unique(spec_li)
spec_li$year_description=as.numeric(unlist(regmatches(spec_li$scientificNameAuthorship, gregexpr("[[:digit:]]+", spec_li$scientificNameAuthorship))))
#check that all species match
spec_li$scientificName[!(spec_li$scientificName %in% Traits2$Species)]
#merge tables
Traits2=merge(Traits2,spec_li,by.x="Species",by.y="scientificName",all.x=FALSE,all.y=TRUE)								


## Include IUCN Status.
IUCN_status <- fread(paste0(project_folder,"data/raw_data/Pollinators_IUCN.csv"))

#check that all species match
Traits2$Species[!(Traits2$Species %in% IUCN_status$TAXON)]

############ clean taxonomy of IUCN FILE
## Change some hoverfly names mannually.----
name_changes <- c(
  "Syrphus niditifrons" = "Syrphus nitidifrons",
  "Cheilosia chrysocoma" = "Cheilosia chrysocomus",
  "Cheilosia himantopa" = "Cheilosia himantopus",
  "Cheilosia urbana" = "Cheilosia urbanus",
  "Chrysosyrphus nigra" = "Chrysosyrphus niger",
  "Metasyrphus frequens" = "Eupeodes frequens",
  "Eriozona erratica" = "Megasyrphus erraticus",
  "Melangyna triangulifera" = "Meligramma triangulifera",
  "Meligramma cingulata" = "Meligramma cingulatum",
  "Meligramma Meligramma guttata" = "guttatum",
  "Xanthogramma aeginae" = "Philhelius aeginae",
  "Xanthogramma marginale" = "Philhelius marginalis",
  "Xanthogramma pedissequum" = "Philhelius pedissequus",
  "Xanthogramma pilosum" = "Philhelius pilosus",
  "Xanthogramma dives" = "Philhelius dives",
  "Xanthogramma laetum" = "Philhelius laetus",
  "Xanthogramma stackelbergi" = "Philhelius stackelbergi"
)

IUCN_status$TAXON <- ifelse(
  IUCN_status$TAXON %in% names(name_changes),
  name_changes[IUCN_status$TAXON],
  IUCN_status$TAXON
)

#automatic check against GBIF
tab=data.frame(species=Traits2$Species[!(Traits2$Species %in% IUCN_status$TAXON)],family=NA)
for(i in 1:nrow(tab)){
  obj=name_backbone(name=gsub("[^0-9A-Za-z///' ]","",tab$species[i],ignore.case=T),order ="Hymenoptera") #look for the given species
  if(length(obj$rank)>0 & !is.na(tab$species[i])){  #if obj is not empty
    if(obj$status=="SYNONYM"){obj=name_backbone(name=name_usage(key=obj$acceptedUsageKey)$data$scientificName)}
	tab$canon[i]=ifelse(length(obj$canonicalName)>0,obj$canonicalName,NA) 
    tab$family[i]=ifelse(length(obj$family)>0,obj$family,NA) #family name according to the GBIF
	tab$rank[i]=obj$rank
  }
}

IUCN_status=merge(IUCN_status,tab,by.x="TAXON",by.y="canon",all.x=TRUE,all.y=FALSE)

IUCN_status$TAXON2=ifelse(!is.na(IUCN_status$species),IUCN_status$species,IUCN_status$TAXON)

#reamining species that do not match
sort(Traits2$Species[!(Traits2$Species %in% IUCN_status$TAXON2)])

nrow(Traits2)
Traits3=merge(Traits2,IUCN_status[,c("TAXON2","European_Category")],by.x="Species",by.y="TAXON2",all.x=TRUE,all.y=FALSE)
nrow(Traits3)


### Create a new column with the updated category (IUCN_2025).
New_IUCN_2025 <- fread(paste0(project_folder,"data/raw_data/BEES_UICN_EUROPE_2025.csv"))

sort(Traits3$Species[Traits3$Order=="Hymenoptera"][!(Traits3$Species[Traits3$Order=="Hymenoptera"] %in% New_IUCN_2025$TAXON)])
#automatic check against GBIF
tab=data.frame(species=Traits3$Species[Traits3$Order=="Hymenoptera"][!(Traits3$Species[Traits3$Order=="Hymenoptera"] %in% New_IUCN_2025$TAXON)],family=NA)
for(i in 1:nrow(tab)){
  obj=name_backbone(name=gsub("[^0-9A-Za-z///' ]","",tab$species[i],ignore.case=T),order ="Hymenoptera") #look for the given species
  if(length(obj$rank)>0 & !is.na(tab$species[i])){  #if obj is not empty
    if(obj$status=="SYNONYM"){obj=name_backbone(name=name_usage(key=obj$acceptedUsageKey)$data$scientificName)}
	tab$canon[i]=ifelse(length(obj$canonicalName)>0,obj$canonicalName,NA) 
    tab$family[i]=ifelse(length(obj$family)>0,obj$family,NA) #family name according to the GBIF
	tab$rank[i]=obj$rank
  }
}

New_IUCN_2025=merge(New_IUCN_2025,tab,by.x="TAXON",by.y="canon",all.x=TRUE,all.y=FALSE)

New_IUCN_2025$TAXON2=ifelse(!is.na(New_IUCN_2025$species),New_IUCN_2025$species,New_IUCN_2025$TAXON)

names(New_IUCN_2025)[names(New_IUCN_2025)=="European category"]="European_Category_2025"

#reamining species that do not match
sort(Traits3$Species[Traits3$Order=="Hymenoptera"][!(Traits3$Species[Traits3$Order=="Hymenoptera"] %in% New_IUCN_2025$TAXON2)])


#MERGE NEW IUCN FOR BEES WITH TRAITS DATA
nrow(Traits3)
Traits4=merge(Traits3,New_IUCN_2025[,c("TAXON2","European_Category_2025")],by.x="Species",by.y="TAXON2",all.x=TRUE,all.y=FALSE)
nrow(Traits4)

Traits4$European_Category_2025[Traits4$European_Category_2025=="NA"]=NA
Traits4$European_Category_2025[Traits4$European_Category=="NA"]=NA

Traits4 <- Traits4 %>%
  mutate(Sociality_simplified = case_when(
    Sociality %in% c("eusocial", "parasocial", "parasocial_(semisocial)") ~ "social",
    Sociality == "solitary" ~ "solitary",
    Sociality %in% c("cleptoparasite", "facultative_inquiline", "obligate_inquiline") ~ "parasite",
    TRUE ~ NA_character_
  ))



Traits4 <- Traits4 %>%
  mutate(
    European_Category = factor(European_Category,
                               levels = c("DD","LC","NT","VU","EN","CR","RE")),
    Adult_Body_size = factor(Adult_Body_size, 
                             levels = c("small", "small_or_medium", "medium", 
                                        "medium_or_large", "large")),
    Adult_Body_size_num = as.numeric(Adult_Body_size),
    
    Flight_height = fct_recode(Flight_height,
                               "Arboreal" = "arboreal",
                               "Mixed" = "arboreal_and_near_ground",
                               "Near ground" = "near_ground"),
    
    Larval_nutrition = fct_recode(Larval_nutrition,
                                  "Phytophagous (bulbs)" = "phytophagous_bulbs",
                                  "Phytophagous (roots)" = "phytophagous_roots",
                                  "Saprophagous" = "saprophagous",
                                  "Saproxylic" = "saproxylic",
                                  "Zoophagous" ="zoophagous")
  )


Traits4$Larval_diet_breadth[Traits4$Larval_diet_breadth=="na"]=NA


Traits4$taxo_group="bees"
Traits4$taxo_group[Traits4$Order=="Diptera"]="hoverflies"

fwrite(Traits4,paste0(project_folder,"data/final_and_intermediate_outputs/traits_table.csv"))



