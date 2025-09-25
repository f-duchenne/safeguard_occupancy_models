###########################################
###########################################
#' Check for packages and if necessary install into library 
#+ message = FALSE
rm(list=ls())
pkgs <- c("data.table", "dplyr","glmmTMB","sf","ggplot2") 

inst <- pkgs %in% installed.packages()
if (any(inst)) install.packages(pkgs[!inst])
pkg.out <- lapply(pkgs, require, character.only = TRUE)


#defining working folder:
project_folder="C:/Users/Duchenne/Documents/safeguard/"

# Loading data
dat=fread(paste0(project_folder,"data/bees_det_nondet_matrix_common.csv"))

years=c(1960,2010)

species="Bombus lapidarius"

dat$Y=dat[,species,with=FALSE]
dat$Y[dat$Y>0]=1

bidon=subset(dat,year_grouped %in% years)

#SHAPE FILE of grid
# Define the UTM projection for a suitable UTM zone
utm_crs <- st_crs("+proj=utm +zone=32 +ellps=WGS84")
hex_grid=st_read(paste0(project_folder,"data/grid_",50,"KM.shp"),crs=utm_crs,quiet =TRUE)

obj=merge(hex_grid,bidon[,c("site","region_50","MONTH_2","Y","year_grouped")],by.x="gridID_50",by.y="site")


p1=ggplot(data=obj)+geom_sf(fill="orange")+facet_wrap(~year_grouped)

png(paste0(project_folder,"sampling.png"),width=1100,height=900,res=120)
p1
dev.off();


p2=ggplot(data=obj)+geom_sf(aes(fill=as.factor(Y)))+facet_wrap(~year_grouped)+scale_fill_manual(values=c("orange","blue"))+labs(fill="presence/absence")+theme(legend.position="none")

png(paste0(project_folder,"pres-abs.png"),width=1100,height=900,res=120)
p2
dev.off();


p3=ggplot(data=obj)+geom_sf(aes(fill=as.factor(Y)),col=NA)+facet_wrap(~year_grouped+MONTH_2,ncol=6)+scale_fill_manual(values=c("orange","blue"))+labs(fill="detection/non detection")

png(paste0(project_folder,"detect_non_detect.png"),width=1400,height=1100,res=120)
p3
dev.off();


geom_point(data)

