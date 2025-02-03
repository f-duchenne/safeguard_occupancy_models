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

dat=subset(dat,region_50 %in% c("alpine","boreal","atlantic","continental","mediterranean"))

#roughly define sites
dat$site=dat$gridID_50

dat=subset(dat,!is.na(MONTH_2))

#define what is a survey:
dat[,survey:=do.call(paste,.SD),.SDcols = which(names(dat) %in% c("site","YEAR_2","MONTH_2"))]

#calculate the number of species detected for each survey
dat=dat[,list_length:=length(unique(TAXON)),by=survey]

#calculate the number of records detected for each survey
dat=dat[,record_number:=length(TAXON),by=survey] 

#define time period
dat$year_grouped=plyr::round_any(dat$YEAR_2,1,f=ceiling)

#nb of records per period:
dat %>% group_by(year_grouped) %>% count()

#richness per period:
b=dat %>% group_by(year_grouped) %>% summarise(richness=length(unique(TAXON)),survey=length(unique(survey)))

#latitude of records per period:
b=dat %>% group_by(year_grouped) %>% summarise(latitude_avg=mean(LATITUDE),longitude_avg=mean(LONGITUDE))
boxplot(LATITUDE~year_grouped,data=dat)

#removing the sites that have been visited only in one period
count_table_sites=dat %>% group_by(site) %>% summarise(nperiods=length(unique(year_grouped))) #count number of records per species
dat=subset(dat,site %in% subset(count_table_sites,nperiods>1)$site)

#generating the detection/non-detections matrices, over sites and visits
length(unique(dat$TAXON)) #number of species (ncol of the matrix)
length(unique(dat$survey)) #number of survey (nrow of the matrix)

#matrix is way too big, needs to be splitted in two parts
#focusing on common species first
## to avoid to get a too huge matrix, we can put all the rare species (that we can not study) together
dat[,species:=TAXON] #new species column
count_table=dat[, .N,by=species] #count number of records per species

b=dat %>% group_by(species) %>% summarise(nsite=length(unique(site)))
# liste_denis=fread("WP2_3_Species.csv")
# names(liste_denis)[1]="species"
# liste_denis=subset(liste_denis,!is.na(species) & species!="")
# b=merge(liste_denis,b,by=c("species"),all.x=TRUE,all.y=FALSE)
# b$analyzed="not in dataset"
# b$analyzed[b$nsite<10]="no"
# b$analyzed[b$nsite>=10]="yes"
# dim(b)

dat[dat$species %in% subset(count_table,N<1000)$species,species:="others"] #all species with less than 1000 records are classified as "others"
# sp_to_test=c("Andrena agilissima","Andrena strohmella","Halictus scabiosae","Lasioglossum minutulum","Bombus terrestris")
# dat[!(dat$species %in% sp_to_test),species:="others"]

length(unique(dat$species))

#create the matrix
mat1=dcast(dat,survey+list_length+record_number+year_grouped+MONTH_2+time_period+site+region_50~species)

# mat1$log.list.length=log(mat1$list_length)
# mat1=mat1 %>% dplyr::group_by(region_50) %>% mutate(log.list.length.c=log.list.length-mean(log.list.length))

# fwrite(mat1,"det_nondet_matrix_species_test.csv")

#export matrix
fwrite(mat1,"det_nondet_matrix_species_common_50.csv")

#focusing on rare species
dat[,species:=TAXON] #new species column
count_table=dat[, .N,by=species] #count number of records per species
dat[dat$species %in% subset(count_table,N>=1000)$species,species:="others"] #all species with more than 999 records are classified as "others"
dat[dat$species %in% subset(count_table,N<10)$species,species:="others"] #all species with less than 10 records are classified as "others"

#create the matrix
mat2=dcast(dat,survey+list_length+record_number+year_grouped+MONTH_2+time_period+site+region_50~species)

#export second matrix
fwrite(mat2,"det_nondet_matrix_species_rare.csv")



obi=mat1 %>% group_by(year_grouped,region_50,site) %>% summarise(nsurv=length(unique(survey)))
obi2=obi %>% group_by(region_50,year_grouped) %>% summarise(moy=mean(nsurv),med=median(nsurv),quant_low=quantile(nsurv,prob=0.1),quant_high=quantile(nsurv,prob=0.9))
ggplot(obi2,aes(x=year_grouped,y=moy,color=region_50,fill=region_50))+
geom_ribbon(aes(ymax=quant_high,ymin=quant_low),alpha=0.2)+
geom_line()+
#geom_pointrange(aes(ymax=quant_high,ymin=quant_low))+
scale_y_continuous(breaks=c(2,4,6,8,10,12,14,16,18))+
theme_bw()+facet_wrap(~region_50,ncol=5)+theme(panel.grid=element_blank(),legend.position="none",strip.background=element_blank())+
ylab("Number of repetition in sampling per period (2 years)")+geom_hline(yintercept=2,linetype="dashed")

#some tests
mat1$Y=mat1$"Bombus terrestris"
mat1$time_period.num=1
mat1$time_period.num[mat1$time_period=="1941-1960"]=2
mat1$time_period.num[mat1$time_period=="1961-1980"]=3
mat1$time_period.num[mat1$time_period=="1981-2000"]=4
mat1$time_period.num[mat1$time_period=="2001-2020"]=5
mat1$Y[mat1$Y>1]=1 #if many dets, put one
count.table=mat1 %>% group_by(region_50) %>% summarise(nr=sum(Y)) #count number of records per country

mat1$log.list.length=log(mat1$list_length)
mat1=mat1 %>% group_by(region_50) %>% mutate(log.list.length.c=log.list.length-mean(log.list.length))
mat1$YEAR_2=(mat1$YEAR_2-mean(mat1$YEAR_2))/sd(mat1$YEAR_2)
bidon=subset(mat1,region_50=="atlantic")

model1=glmmTMB(Y ~ YEAR_2+(1|site)+(1|MONTH_2)+log.list.length.c,data=bidon,family="binomial")






