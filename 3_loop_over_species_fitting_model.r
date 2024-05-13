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

#defining working folder:
setwd(dir="C:/Users/Duchenne/Documents/safeguard/data")

#species index
index=i+6

#import data:
dat=fread("det_nondet_matrix_species_common.csv")
bidon=subset(dat,COUNTRY=="France")
bidon2=data.frame(long=sapply(strsplit(unique(bidon$site)," "),function(x){x[2]}),lat=sapply(strsplit(unique(bidon$site)," "),function(x){x[1]}))

dat$Y=dat[,index,with=F]
count.table=dat[,sum(Y),by=COUNTRY] #count number of records per country

#exclude country without record of the focal species
dat=subset(dat,COUNTRY %in% subset(count.table,V1>0)$COUNTRY)

#det/nondet of the focal species
Y=dat$Y
Y[Y>1]=1 #if many dets, put one

N=nrow(dat)

Nc=length(unique(dat$COUNTRY))

Nperiod=length(unique(dat$time_period))

Nsite=length(unique(dat$site))

period.num=as.numeric(as.factor(dat$time_period))

country.num=as.numeric(as.factor(dat$COUNTRY))

site.num=as.numeric(as.factor(dat$site))

dat.model=list(Y=Y,N=N,Nc=Nc,Nperiod=Nperiod,Nsite=Nsite,time_period=dat$time_period,period.num=period.num,COUNTRY=dat$COUNTRY,country.num=country.num,list.length=dat$list_length,site=dat$site,site.num=site.num)

country_tab=unique(data.frame(COUNTRY=dat.model$COUNTRY,country.num=country.num))
period_tab=unique(data.frame(time_period=dat.model$time_period,period.num=period.num))
species=names(dat)[index]

rm(dat)

ParsStage <- c("country.period.det","country.period.occ","alpha","site.eff","sd.occ","sd.det","sd.site")

Inits <- function(){list()}

t1=Sys.time()
results1 <- jags.parallel(model.file="model.txt", parameters.to.save=ParsStage, n.chains=3, data=dat.model,n.burnin = 20000,  n.iter = 30000, n.thin = 3,inits =Inits,jags.seed =2)
t2=Sys.time()
t2-t1

lili=list(results1,country_tab,period_tab,species)

save(lili,file=paste0("model_",i,".RData"))