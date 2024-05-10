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
dat=readRDS("Safeguard_bee_df_final_v8_20240509.rds")

#convert data.frame to data.table for more faster subsequent manipulations
dat=as.data.table(dat)

#define time periods
dat[,time_period:=ifelse(.SD<=1980,"1980",ifelse(.SD<=1990,"1990",ifelse(.SD<=2000,"2000",ifelse(.SD<=2010,"2010","2020")))),.SDcols = "YEAR_2"]

#roughly define sites (first approximation, sites should not be defined in WGS84 coordinates but using a coordinate system that conserve distances over latitude)
dat$LATITUDE=as.numeric(dat$LATITUDE)
dat$LONGITUDE=as.numeric(dat$LONGITUDE)
dat[,site:=do.call(paste,do.call(round,list(x=.SD,digits=1))),.SDcols = which(names(dat) %in% c("LATITUDE","LONGITUDE"))]

#define what is a survey:
dat[,survey:=do.call(paste,.SD),.SDcols = which(names(dat) %in% c("site","YEAR","MONTH_2","LOCALITY"))]

# remove data without time_period (because year is missing)
dat=subset(dat,!is.na(time_period))

#use sparta tools to have a first look at the dataset
results <- dataDiagnostics(taxa = dat$TAXON,
                           site = dat$site,
                           time_period = dat$time_period,
                           progress_bar = TRUE)

#nb of records per period:
dat %>% group_by(time_period) %>% count()

#richness per period:
dat %>% group_by(time_period) %>% summarise(richness=length(unique(TAXON)))

#latitude of records per period:
dat %>% group_by(time_period) %>% summarise(latitude_avg=mean(LATITUDE),longitude_avg=mean(LONGITUDE))



#generating the detection/non-detections matrices, over sites and visits
## to avoid to get a too huge matrix, we can put all the rare species (that we can not study) together
dat[,species:=TAXON]
count_table=dat[, .N,by=species]
dat[dat$species %in% subset(count_table,N<50)$species,species="others"]
length(unique(dat$species))
length(unique(dat$survey))
mat=formatOccData(taxa=dat$TAXON,site=dat$site,survey=dat$survey,closure_period=dat$time_period)

