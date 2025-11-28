###########################################
#' Check for packages and if necessary install into library 
#+ message = FALSE
rm(list=ls())
pkgs <- c("data.table", "dplyr","sf","tidyverse") 

inst <- pkgs %in% installed.packages()
if (any(inst)) install.packages(pkgs[!inst])
pkg_out <- lapply(pkgs, require, character.only = TRUE)

#defining working folder:
#project_folder="C:/Users/Duchenne/Documents/safeguard/"
project_folder=""

# Loading data
dat=fread(paste0(project_folder,"data/final_and_intermediate_outputs/database_clean_filtered.csv"))
nb_records_initial=nrow(dat)
hex_grid=st_read(paste0(project_folder,"data/raw_data/grids_shapefiles/grid_",50,"KM.shp"),crs="+proj=utm +zone=32 +ellps=WGS84")
hex_grid_p=st_centroid(hex_grid)
gride=cbind(data.frame(gridID_50=hex_grid_p$gridID_50),st_coordinates(hex_grid_p))
names(gride)[2:3]=c("long_50","lat_50")
nrow(dat)
dat=merge(dat,gride,by=c("gridID_50"))
nrow(dat)


subset(dat,region_50 %in% c("alpine","boreal","atlantic","continental","mediterranean")) %>% group_by(taxo_group) %>% count()

n1=nrow(dat)
dat=subset(dat,region_50 %in% c("alpine","boreal","atlantic","continental","mediterranean"))
nr_regions=n1-nrow(dat)
nb_records=nrow(dat)
dat %>% group_by(taxo_group) %>% count()

#roughly define sites
dat$site=dat$gridID_50

b=dat %>% dplyr::group_by(gridID_50,region_50,taxo_group) %>% dplyr::summarise(nyear=length(unique(endYear)),
nyear_b70=sum(unique(endYear)<=1970),nyear_a70=sum(unique(endYear)>1970),nb_records=length(endYear))

fwrite(b,paste0(project_folder,"data/final_and_intermediate_outputs/sites_studied_tot.csv"))

nr_month=subset(dat,is.na(endMonth))
nr_month %>% group_by(taxo_group) %>% count()
dat=subset(dat,!is.na(endMonth))

#define what is a survey:
dat[,survey:=do.call(paste,.SD),.SDcols = which(names(dat) %in% c("site","endYear","endMonth","taxo_group"))]

#calculate the number of species detected for each survey
dat=dat[,list_length:=length(unique(scientificName)),by=survey]

#calculate the number of records detected for each survey
dat=dat[,record_number:=length(scientificName),by=survey] 

#define time period
dat$year_grouped=plyr::round_any(dat$endYear,1,f=ceiling)

#nb of records per period:
dat %>% group_by(year_grouped) %>% count()

#richness per period:
b=dat %>% group_by(year_grouped) %>% summarise(richness=length(unique(scientificName)),survey=length(unique(survey)))

#latitude of records per period:
b=dat %>% group_by(year_grouped) %>% summarise(latitude_avg=mean(decimalLatitude),longitude_avg=mean(decimalLongitude))
#boxplot(LATITUDE~year_grouped,data=dat)

#removing the sites that have been visited only in one period
count_table_sites=dat %>% group_by(site,taxo_group) %>% summarise(nperiods=length(unique(year_grouped))) 
nr_sites=nrow(subset(count_table_sites,nperiods<2))
subset(count_table_sites,nperiods<2) %>% group_by(taxo_group) %>% count()
dat=subset(dat,site %in% subset(count_table_sites,nperiods>1)$site)
nb_surveys=length(unique(dat$survey))

#generating the detection/non-detections matrices, over sites and visits
length(unique(dat$scientificName)) #number of species (ncol of the matrix)
length(unique(dat$survey)) #number of survey (nrow of the matrix)


###EXPORT NUMBERS AND TABLES THAT WILL BE USEFUL LATTER
count_table=dat[, .N,by=c("scientificName","family","taxo_group")] #count number of records per species

bb=dat %>% dplyr::group_by(gridID_50,region_50,scientificName,family,survey,taxo_group) %>% dplyr::summarise(occupied=length(scientificName))
bb=bb %>% dplyr::group_by(region_50,taxo_group) %>% dplyr::mutate(nsurv=length(unique(survey)))

b= bb %>% dplyr::group_by(region_50,scientificName,family,taxo_group) %>% dplyr::summarise(occupancy_obs=sum(occupied>0)/mean(nsurv),nb_records=sum(occupied),nb_detect=sum(occupied>0))
b=b %>% group_by(scientificName,taxo_group,family) %>% mutate(nb_records_tot=sum(nb_records))

fwrite(b,paste0(project_folder,"data/final_and_intermediate_outputs/species_nb_records.csv"))

nb_sp_common=subset(count_table,N>=1000) %>% group_by(taxo_group) %>% count() %>%  deframe()
nb_sp=subset(b,nb_records_tot>=10 & nb_detect>=5) %>% group_by(taxo_group) %>% summarise(n=length(unique(scientificName))) %>%  deframe()
nsp_tot=b %>% group_by(taxo_group) %>% summarise(n=length(unique(scientificName))) %>%  deframe()
length(unique(dat$scientificName))

list_filtering=list(nb_records_initial,nr_regions,nb_records,nr_month,nr_sites,nb_surveys,nb_sp_common,nb_sp,nsp_tot)
save(list_filtering,file=paste0(project_folder,"data/final_and_intermediate_outputs/list_filtering.RData"))

b=dat %>% dplyr::group_by(gridID_50,region_50,taxo_group) %>% dplyr::summarise(nyear=length(unique(endYear)),
nyear_b70=sum(unique(endYear)<=1970),nyear_a70=sum(unique(endYear)>1970))

#commented to avoid overwrite
#fwrite(b,paste0(project_folder,"data/final_and_intermediate_outputs/sites_studied.csv"))

bb=dat %>% dplyr::group_by(region_50,endYear,gridID_50,taxo_group) %>% dplyr::summarise(nsurveys=length(unique(survey)),nrec=length(survey))
b=bb %>% dplyr::group_by(region_50,endYear,taxo_group) %>% dplyr::summarise(nsites=length(unique(gridID_50)),mean_nsurv=mean(nsurveys),nrec=sum(nrec))

fwrite(b,paste0(project_folder,"data/final_and_intermediate_outputs/regions_sampling.csv"))

######################################### PREPARE MATRIX OF DETECTION AND NON DETECTION FOR EACH GROUP
taxo_group_vec=c("bees","hoverflies")
for(j in taxo_group_vec){
	dat2=subset(dat,taxo_group==j)
	
	if(j=="bees"){
		count_table=dat2[, .N,by=c("scientificName","family")] #count number of records per species
		#matrix is way too big, needs to be splitted in two parts
		#focusing on common species first
		## to avoid to get a too huge matrix, we can put all the rare species (that we can not study) together
		dat2[,species:=scientificName] #new species column
		dat2[dat2$species %in% subset(count_table,N<1000)$scientificName,species:="others"] #all species with less than 1000 records are classified as "others"
		# sp_to_test=c("Andrena agilissima","Andrena strohmella","Halictus scabiosae","Lasioglossum minutulum","Bombus terrestris")
		# dat[!(dat$species %in% sp_to_test),species:="others"]

		#create the matrix
		mat1=dcast(dat2,survey+list_length+record_number+year_grouped+endMonth+time_period+site+long_50+lat_50+region_50~species)

		# mat1$log.list.length=log(mat1$list_length)
		# mat1=mat1 %>% dplyr::group_by(region_50) %>% mutate(log.list.length.c=log.list.length-mean(log.list.length))

		# fwrite(mat1,"det_nondet_matrix_species_test.csv")

		#export matrix
		fwrite(mat1,paste0(project_folder,"data/final_and_intermediate_outputs/",j,"_det_nondet_matrix_common.csv"))

		#focusing on rare species
		count_table=dat2[, .N,by=c("scientificName","family")] #count number of records per species
		dat2[,species:=scientificName] #new species column
		dat2[dat2$species %in% subset(count_table,N>=1000)$scientificName,species:="others"] #all species with more than 999 records are classified as "others"
		dat2[dat2$species %in% subset(count_table,N<10)$scientificName,species:="others"] #all species with less than 10 records are classified as "others"

		#create the matrix
		mat2=dcast(dat2,survey+list_length+record_number+year_grouped+endMonth+time_period+site+long_50+lat_50+region_50~species)

		#export second matrix
		fwrite(mat2,paste0(project_folder,"data/final_and_intermediate_outputs/",j,"_det_nondet_matrix_rare.csv"))
	}else{
		dat2[,species:=scientificName] #new species column
		#create the matrix
		mat1=dcast(dat2,survey+list_length+record_number+year_grouped+endMonth+time_period+site+long_50+lat_50+region_50~species)

		# mat1$log.list.length=log(mat1$list_length)
		# mat1=mat1 %>% dplyr::group_by(region_50) %>% mutate(log.list.length.c=log.list.length-mean(log.list.length))

		# fwrite(mat1,"det_nondet_matrix_species_test.csv")

		#export matrix
		fwrite(mat1,paste0(project_folder,"data/final_and_intermediate_outputs/",j,"_det_nondet_matrix.csv"))
	}

}



dat1=fread(paste0(project_folder,"data/final_and_intermediate_outputs/bees_det_nondet_matrix_common.csv"))
dat2=fread(paste0(project_folder,"data/final_and_intermediate_outputs/bees_det_nondet_matrix_rare.csv"))
dat3=fread(paste0(project_folder,"data/final_and_intermediate_outputs/hoverflies_det_nondet_matrix.csv"))


length((which(names(dat1)=="region_50")+1):(which(names(dat1)=="others")-1))

liste_species=rbind(
data.frame(species=names(dat1)[(which(names(dat1)=="region_50")+1):(which(names(dat1)=="others")-1)],taxo_group="bees",matrix="common"),
data.frame(species=names(dat2)[(which(names(dat2)=="region_50")+1):(which(names(dat2)=="others")-1)],taxo_group="bees",matrix="rare"),
data.frame(species=names(dat3)[(which(names(dat3)=="region_50")+1):ncol(dat3)],taxo_group="hoverflies",matrix=NA))
fwrite(liste_species,paste0(project_folder,"data/final_and_intermediate_outputs/liste_total_species_occupancy.csv"))
