###########################################
#' Check for packages and if necessary install into library 
#+ message = FALSE
rm(list=ls())
pkgs <- c("data.table", "dplyr","sparta") 

inst <- pkgs %in% installed.packages()
if (any(inst)) install.packages(pkgs[!inst])
pkg_out <- lapply(pkgs, require, character.only = TRUE)

#defining working folder:
setwd(dir="C:/Users/Duchenne/Documents/safeguard/data")

# Loading data
dat=fread("database_clean_filtered.csv")
taxi=fread("species_family_table,.csv")
dat=merge(dat,taxi,by="TAXON")

dat=subset(dat,region_20 %in% c("alpine","boreal","atlantic","continental","mediterranean"))

#roughly define sites
dat$site=dat$gridID_20

#define what is a survey:
dat[,survey:=do.call(paste,.SD),.SDcols = which(names(dat) %in% c("site","YEAR_2","MONTH_2"))]

#calculate the number of species detected for each survey
dat=dat[,list_length:=length(unique(TAXON)),by=survey] 

#nb of records per period:
dat %>% group_by(time_period) %>% count()

#richness per period:
b=dat %>% group_by(time_period) %>% summarise(richness=length(unique(TAXON)),survey=length(unique(survey)))

#latitude of records per period:
b=dat %>% group_by(time_period) %>% summarise(latitude_avg=mean(LATITUDE),longitude_avg=mean(LONGITUDE))
boxplot(LATITUDE~time_period,data=dat)

#removing the sites that have been visited only in one period
count_table_sites=dat %>% group_by(site) %>% summarise(nperiods=length(unique(time_period))) #count number of records per species
dat=subset(dat,site %in% subset(count_table_sites,nperiods>1)$site)

#generating the detection/non-detections matrices, over sites and visits
length(unique(dat$TAXON)) #number of species (ncol of the matrix)
length(unique(dat$survey)) #number of survey (nrow of the matrix)

#matrix is way too big, needs to be splitted in two parts
#focusing on common species first
## to avoid to get a too huge matrix, we can put all the rare species (that we can not study) together
dat[,species:=TAXON] #new species column
count_table=dat[, .N,by=species] #count number of records per species
dat[dat$species %in% subset(count_table,N<1000)$species,species:="others"] #all species with less than 1000 records are classified as "others"
length(unique(dat$species))


#create the matrix
mat1=dcast(dat,survey+list_length+YEAR_2+time_period+site+COUNTRY+region_20~species)

#export matrix
fwrite(mat1,"det_nondet_matrix_species_common.csv")

#focusing on rare species
dat[,species:=TAXON] #new species column
count_table=dat[, .N,by=species] #count number of records per species
dat[dat$species %in% subset(count_table,N>=1000)$species,species:="others"] #all species with more than 999 records are classified as "others"

#create the matrix
mat2=dcast(dat,survey+list_length+YEAR_2+time_period+site+COUNTRY2~species)

#export second matrix
fwrite(mat2,"det_nondet_matrix_species_rare.csv")




