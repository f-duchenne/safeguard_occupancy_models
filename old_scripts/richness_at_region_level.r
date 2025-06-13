###########################################
#' Check for packages and if necessary install into library 
#+ message = FALSE
rm(list=ls())
pkgs <- c("data.table", "dplyr","sf","tidyverse","iNEXT") 

inst <- pkgs %in% installed.packages()
if (any(inst)) install.packages(pkgs[!inst])
pkg_out <- lapply(pkgs, require, character.only = TRUE)

#defining working folder:
setwd(dir="C:/Users/Duchenne/Documents/safeguard/data")

# Loading data
dat=fread("database_clean_filtered.csv")
nb_records_initial=nrow(dat)
hex_grid=st_read(paste0("grid_",50,"KM.shp"),crs="+proj=utm +zone=32 +ellps=WGS84")
hex_grid_p=st_centroid(hex_grid)
gride=cbind(data.frame(gridID_50=hex_grid_p$gridID_50),st_coordinates(hex_grid_p))
names(gride)[2:3]=c("long_50","lat_50")
nrow(dat)
dat=merge(dat,gride,by=c("gridID_50"))
nrow(dat)

n1=nrow(dat)
dat=subset(dat,region_50 %in% c("alpine","boreal","atlantic","continental","mediterranean"))
nr_regions=n1-nrow(dat)
nb_records=nrow(dat)

#roughly define sites
dat$site=dat$gridID_50

b=dat %>% dplyr::group_by(gridID_50,region_50,taxo_group) %>% dplyr::summarise(nyear=length(unique(YEAR_2)),
nyear_b70=sum(unique(YEAR_2)<=1970),nyear_a70=sum(unique(YEAR_2)>1970),nb_records=length(YEAR_2))

fwrite(b,"sites_studied_tot.csv")

nr_month=nrow(subset(dat,is.na(MONTH_2)))
dat=subset(dat,!is.na(MONTH_2))

#define what is a survey:
dat[,survey:=do.call(paste,.SD),.SDcols = which(names(dat) %in% c("site","YEAR_2","taxo_group"))]

#calculate the number of species detected for each survey
dat=dat[,list_length:=length(unique(TAXON)),by=survey]

#calculate the number of records detected for each survey
dat=dat[,record_number:=length(TAXON),by=survey] 

#define time period
dat$year_grouped=plyr::round_any(dat$YEAR_2,1,f=ceiling)
dat$abund=1

reg="mediterranean"
tax_group="bees"
for(reg in unique(dat$region_50)){
for(tax_group in c("bees")){
b=subset(dat,taxo_group==tax_group & region_50==reg)
comm=split(b, b$year_grouped)
comm=lapply(comm,function(x){dcast(x,year_grouped+TAXON~survey,fun.aggregate =max,value.var="abund",fill=0)}) #put each data frame in the right format for iNEXT function
comm=lapply(comm,function(x){x[,-c(1:2)]}) #remove id columns
comm=lapply(comm,as.data.frame)

rari=iNEXT(comm,q=0,datatype ="incidence_raw")
save(rari,file=paste0("iNEXT_",tax_group,"_",reg,".RData"))
}}


resf=NULL
for(reg in unique(dat$region_50)){
for(tax_group in c("bees")){
load(paste0("iNEXT_",tax_group,"_",reg,".RData"))

bidon=subset(rari$AsyEst,Diversity=="Species richness")
bidon$taxo_group=tax_group
bidon$region_50=reg
resf=rbind(resf,bidon)
}}

colo2=c("#44AA99","#117733","#332288","#CC6677","#DDCC77")

ggplot(data=resf,aes(y=Estimator,x=as.numeric(as.character(Assemblage)),color=region_50))+geom_pointrange(aes(ymin=LCL,ymax=UCL))+facet_wrap(~region_50,scales="free")+
theme_bw()+
theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.border=element_blank(),
panel.grid.minor = element_blank(),panel.background = element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),
strip.background=element_rect(fill=NA,color=NA),legend.position="none")+
scale_color_manual(values=colo2)