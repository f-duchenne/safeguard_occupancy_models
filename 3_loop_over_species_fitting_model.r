###########################################
###########################################
#' Check for packages and if necessary install into library 
#+ message = FALSE
rm(list=ls())
pkgs <- c("data.table", "dplyr","R2jags") 

inst <- pkgs %in% installed.packages()
if (any(inst)) install.packages(pkgs[!inst])
pkg.out <- lapply(pkgs, require, character.only = TRUE)


# Collect command arguments
args <- commandArgs(trailingOnly = TRUE)
args_contents <- strsplit(args, ' ')
# Get first argument
i <- as.numeric(args_contents[[1]])
i=1

#defining working folder:
setwd(dir="C:/Users/Duchenne/Documents/safeguard/data")

#species index
index=i+7

#import data:
dat=fread("det_nondet_matrix_species_common.csv")


dat$Y=dat[,index,with=F]
#det/nondet of the focal species
Y=dat$Y
Y[Y>1]=1 #if many dets, put one
count.table=dat[,sum(Y),by=region_20] #count number of records per country

nsurvey_tot=nrow(dat)
#exclude country without record of the focal species
dat=subset(dat,region_20 %in% subset(count.table,V1>0)$region_20)
nsurvey_used=nrow(dat)


N=nrow(dat)

Nc=length(unique(dat$COUNTRY2))

Nperiod=length(unique(dat$time_period))

Nr=length(unique(dat$region_20))

Nsite=length(unique(dat$site))

period.num=as.numeric(as.factor(dat$time_period))

country.num=as.numeric(as.factor(dat$COUNTRY))

region.num=as.numeric(as.factor(dat$region_20))

site.num=as.numeric(as.factor(dat$site))

region.site=unique(as.data.table(cbind(site.num,region.num)))
region.site=region.site[order(region.site$site.num),]

dat.model=list(Y=Y,N=N,Nc=Nc,Nperiod=Nperiod,Nsite=Nsite,Nr=Nr,time_period=dat$time_period,period.num=period.num,COUNTRY=dat$COUNTRY,region_20=dat$region_20,region.num=region.num,country.num=country.num,list.length=dat$list_length,site=dat$site,site.num=site.num,region.site_num=region.site$region.num)

country_tab=unique(data.frame(COUNTRY=dat.model$COUNTRY,country.num=country.num))
period_tab=unique(data.frame(time_period=dat.model$time_period,period.num=period.num))
region_tab=unique(data.frame(region_20=dat.model$region_20,region.num=region.num))

species=names(dat)[index]

rm(dat)

ParsStage <- c("country.period.det","country.period.occ","alpha","eu_eff","sd.occ","sd.det","sd.site")

Inits <- function(){list()}

t1=Sys.time()
results1 <- jags.parallel(model.file="model.txt", parameters.to.save=ParsStage, n.chains=1, data=dat.model,n.burnin = 2,  n.iter = 20, n.thin = 1,inits =Inits,jags.seed =2)
t2=Sys.time()
t2-t1

lili=list(results1,country_tab,period_tab,species,nsurvey_tot,nsurvey_used)

save(lili,file=paste0("model_",i,".RData"))