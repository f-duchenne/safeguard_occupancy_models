# Install packages
library(rgbif)
library(tidyverse)
library(data.table)
library(CoordinateCleaner)


trendsf=fread(paste0("data/final_and_intermediate_outputs/all_trends.csv"))
dat=data.frame(species=unique(trendsf$species))

for(i in 1:nrow(dat)){
  obj=name_backbone(name=gsub("[^0-9A-Za-z///' ]","",dat$species[i],ignore.case=T),class="Insecta") #look for the given species
  if(length(obj$rank)>0 & !is.na(dat$species[i])){  #if obj is not empty
    if(obj$status=="SYNONYM"){obj=name_backbone(name=name_usage(key=obj$acceptedUsageKey)$data$scientificName)}
    dat$rang[i]=obj$rank  #keep the rank of the match
    dat$confi[i]=obj$confidence #confidence in the match
    dat$kingdom[i]=obj$kingdom #kingdom
    dat$speciesKey[i]=ifelse(obj$rank %in% c("SPECIES","SUBSPECIES"),obj$speciesKey,NA)
    dat$usageKey[i]=obj$usageKey
    dat$canon[i]=ifelse(obj$rank %in% c("SPECIES"),obj$canonicalName,NA) #latin species name according to the GBIF
    dat$canon[i]=ifelse(obj$rank %in% c("SUBSPECIES"),obj$species,dat$canon[i]) #latin species name according to the GBIF
    dat$family[i]=ifelse(length(obj$family)>0,obj$family,NA) #family name according to the GBIF
    dat$genus[i]=ifelse(length(obj$genus)>0,obj$genus,NA) #genus name according to the GBIF
    if(length(obj$order)>0){dat$order[i]=obj$order}else{dat$order[i]=NA} #keep the order according to the GBIF
    if(length(obj$class)>0){dat$class[i]=obj$class}else{dat$class[i]=NA} #genus name according to the GBIF
  }
}

# Getting the taxon keys for all bombus
dat=subset(dat,rang=="SPECIES")
species_keys <- dat$speciesKey
predobj <- pred_and( pred_within("POLYGON((34 34, 34 75, -15 75, -15 34, 34 34))"),
                     pred_in("speciesKey",species_keys),
                     pred("hasCoordinate", TRUE),
                     pred("hasGeospatialIssue", FALSE),
                     pred_gte("year", 1921),
                     pred_lte("year", 2021))

cle <- occ_download(predobj,user = "fduchenn", pwd = "Afud2og*", email = "duchenne.f@gmail.com")

occ_download_wait(cle[1]) #wait that download is ready
#download when ready:
my_download <- occ_download_get(cle[1], overwrite=T, path = paste0("data/final_and_intermediate_outputs/."))

# let's keep useful information for later and for publication, the unique key and DOI of the download, that allow reproducing exactly the same extrcation, even in 10 years:
cle=occ_download_doi("10.15468/dl.wpvffa")
my_down=occ_download_meta(cle$key)

downloads=data.frame(key= cle$key,DOI = gbif_citation(my_down)$download)

fwrite(downloads,paste0("data/final_and_intermediate_outputs/key_and_DOI", Sys.Date(), ".csv"))

################# OPENING the occurrence data filter columns and re-export

b=occ_download_import(key = cle$key,select=c("gbifID","day","month","year",
                                             "decimalLatitude","decimalLongitude","speciesKey","usageKey",
                                             "coordinateUncertaintyInMeters","coordinatePrecision","occurrenceStatus",
                                             "countryCode", "individualCount","genus", "family", "taxonRank", "datasetName",
                                             "institutionCode","basisOfRecord","datasetID","occurrenceID",
                                             "institutionID","publisher","acceptedScientificName"),
                        path = paste0("data/final_and_intermediate_outputs/.")) #with such amount of data we need to reduce data frame dimensions to decrease the memory usage, we keep only useful columns, I tried to select as few as possible, but maybe you would need other ones


b$isPseudodata=FALSE
b$endYear=b$year
b$endMonth=b$month
b$endDay=b$day
b$datasetProvider="GBIF"
b$country=b$countryCode
b$scientificNameAuthorship=b$acceptedScientificName
b$speciesKey=as.character(b$speciesKey)
dat=dat[!(dat$speciesKey %in% dat$speciesKey[duplicated(dat$speciesKey)]),]
b=merge(b,unique(dat[,c("speciesKey","species")]),by=c("speciesKey"))
b$scientificName=b$species

b_filtered =  subset(b,is.na(coordinateUncertaintyInMeters) | coordinateUncertaintyInMeters < 15000)

b_filtered <- b_filtered %>%
  # Remove invalid coordintes (outside 90, 180, etc)
  cc_val() %>%
  # Test if coordinates = 0 
  cc_zero() %>%
  # Test if coordinates equal
  cc_equ() %>%
  # This removes if the coordinates are the same as the country centroid center 
  cc_cen(test = "country") %>%
  # This removes duplicates that have the same coordinates, day, month, year
  cc_dupl(additions = c("day", "month", "year")) %>%
  # This removes occurences that are less than 500 m from the capital coordinates
  cc_cap(buffer = 500) 


fwrite(b_filtered,paste0("data/final_and_intermediate_outputs/extraction_gbif.csv")) #write the light extraction
