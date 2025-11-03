## Assemble Traits data
#this script explores the database for descriptive analysis.
pkgs <- c("data.table", "dplyr","tidyverse","rgbif") 
inst <- pkgs %in% installed.packages()
if (any(inst)) install.packages(pkgs[!inst])
pkg_out <- lapply(pkgs, require, character.only = TRUE)


#defining working folder:
project_folder="C:/Users/Duchenne/Documents/safeguard/"

select <- dplyr::select

Traits <- fread(paste0(project_folder,"data/hoverfly_bee_traits_2025_04_01.csv"))
#correct a small mistake of names:
Traits$Species[Traits$Species=="Syrphus niditifrons"]="Syrphus nitidifrons"

#select useful columns
Traits2 <- Traits %>% select(Species,Order,Family,Genus, Sociality, STI_Species_temperature_index, Larval_diet_breadth, ITD_F_mm, 
								Adult_Body_size, Larval_nutrition,Flight_ability,Flight_height)
								
# ADDING DESCRIPTION DATE FOR EACH SPECIES AND KEEP ONLY SPECIES THAT WE HAVE IN OUR DATASET
datf=fread(paste0(project_folder,"data/database_clean_filtered.csv"))
spec_li=datf[,c("scientificNameAuthorship","scientificName")]
spec_li=unique(spec_li)
spec_li$year_description=as.numeric(unlist(regmatches(spec_li$scientificNameAuthorship, gregexpr("[[:digit:]]+", spec_li$scientificNameAuthorship))))
#check that all species match
spec_li$scientificName[!(spec_li$scientificName %in% Traits2$Species)]
#merge tables
Traits2=merge(Traits2,spec_li,by.x="Species",by.y="scientificName",all.x=FALSE,all.y=TRUE)								


## Include IUCN Status.
IUCN_status <- fread(paste0(project_folder,"data/Pollinators_IUCN.csv"))

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
New_IUCN_2025 <- fread(paste0(project_folder,"data/BEES_UICN_EUROPE_2025.csv"))

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

Traits4$taxo_group="bees"
Traits4$taxo_group[Traits4$Order=="Diptera"]="hoverflies"

Traits4$Larval_diet_breadth[Traits4$Larval_diet_breadth=="na"]=NA

fwrite(Traits4,paste0(project_folder,"data/traits_table.csv"))



