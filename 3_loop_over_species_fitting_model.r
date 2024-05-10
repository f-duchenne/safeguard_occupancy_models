###########################################
###########################################
#' Check for packages and if necessary install into library 
#+ message = FALSE
rm(list=ls())
pkgs <- c("data.table", "dplyr","R2jags") 

inst <- pkgs %in% installed.packages()
if (any(inst)) install.packages(pkgs[!inst])
pkg_out <- lapply(pkgs, require, character.only = TRUE)

#defining working folder:
setwd(dir="C:/Users/Duchenne/Documents/safeguard/data")

#species index
i=1
index=i+6

#import data:
dat=fread("det_nondet_matrix_species_common.csv")
dat$Y=dat[,index,with=F]
count_table=dat[,sum(Y),by=COUNTRY] #count number of records per country

#exclude country without record of the focal species
dat=subset(dat,COUNTRY %in% subset(count_table,V1>0)$COUNTRY)

#det/nondet of the focal species
Y=dat$Y
Y[Y>1]=1 #if many dets, put one

N=nrow(dat)

Nc=length(unique(dat$COUNTRY))

Nperiod=length(unique(dat$time_period))

period_num=as.numeric(as.factor(dat$time_period))

country_num=as.numeric(as.factor(dat$COUNTRY))

site_num=as.numeric(as.factor(dat$site))

dat.model=list(Y=Y,N=N,Nc=Nc,Nperiod=Nperiod,time_period=dat$time_period,period_num=period_num,COUNTRY=dat$COUNTRY,country_num=country_num,list_length=dat$list_length,site=dat$site,site_num=site_num)

ParsStage <- c("country_period_det","country_period_det","alpha","site_eff","sd.occ","sd.det","sd.site")

Inits <- function(){list()}

t1=Sys.time()
results1 <- jags(model.file="model.txt", parameters.to.save=ParsStage, n.chains=1, data=dat.model,n.burnin = 50,  n.iter = 100, n.thin = 2,inits =Inits,jags.seed =2)
t2=Sys.time()
t2-t1








save(results1,file=paste0("chain_model_ZI",j,".RData"))
#
