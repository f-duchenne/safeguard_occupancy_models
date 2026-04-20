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
# Capture the arguments passed from the command line
args <- commandArgs(trailingOnly = TRUE)
i <- as.numeric(args[1])+1
print(i)

#defining working folder:
setwd(dir="/data/")


#import data:
dat=fread("det_nondet_matrix_species_common_50.csv")

#combinations
tab=expand.grid(chain=1:3,species=names(dat)[8:(which(names(dat)=="others")-1)])
index=names(dat)[names(dat)==tab$species[i]]

print(index)

dat$Y=dat[,index,with=F]
dat$Y[dat$Y>1]=1 #if many dets, put one
#det/nondet of the focal species
Y=dat$Y
count.table=dat[,sum(Y),by=region_50] #count number of records per country

nsurvey_tot=nrow(dat)
#exclude country without record of the focal species
dat=subset(dat,region_50 %in% subset(count.table,V1>0)$region_50)
nsurvey_used=nrow(dat)

dat$log.list.length=log(dat$list_length)
dat=dat %>% group_by(region_50) %>% mutate(log.list.length.c=log.list.length-mean(log.list.length))

N=nrow(dat)

print(N)

Nc=length(unique(dat$COUNTRY))

Nperiod=length(unique(dat$time_period))

Nr=length(unique(dat$region_50))

Nsite=length(unique(dat$site))

dat$time_period.num=1
dat$time_period.num[dat$time_period=="1941-1960"]=2
dat$time_period.num[dat$time_period=="1961-1980"]=3
dat$time_period.num[dat$time_period=="1981-2000"]=4
dat$time_period.num[dat$time_period=="2001-2020"]=5

period.num=dat$time_period.num

country.num=as.numeric(as.factor(dat$COUNTRY))

region.num=as.numeric(as.factor(dat$region_50))

site.num=as.numeric(as.factor(dat$site))

region.site=unique(as.data.table(cbind(site.num,region.num)))
region.site=region.site[order(region.site$site.num),]

dat.model=list(Y=Y,N=N,Nc=Nc,Nperiod=Nperiod,Nsite=Nsite,Nr=Nr,time_period=dat$time_period,period.num=period.num,COUNTRY=dat$COUNTRY,region_50=dat$region_50,region.num=region.num,country.num=country.num,log.list.length.c=dat$log.list.length.c,site=dat$site,site.num=site.num,region.site_num=region.site$region.num)

print(N)

country_tab=unique(data.frame(COUNTRY=dat.model$COUNTRY,country.num=country.num))
period_tab=unique(data.frame(time_period=dat.model$time_period,period.num=period.num))
region_tab=unique(data.frame(region_50=dat.model$region_50,region.num=region.num))

lili=list(country_tab,period_tab,region_tab,nsurvey_tot,nsurvey_used,region.site)
if(tab$chain[i]==1){
save(lili,file=paste0("data_",tab$species[i],".RData"))
}

species=names(dat)[index]

rm(dat)

ParsStage <- c("region.period.det","region.period.occ","alpha","eu_eff","sd.occ","sd.det","sd.site")

Inits <- function(){list()}

t1=Sys.time()
results1 <- jags(model.file="model.txt", parameters.to.save=ParsStage, n.chains=1, data=dat.model,n.burnin = 110000,  n.iter = 120000, n.thin = 3,inits =Inits,jags.seed=i)
t2=Sys.time()
t2-t1


save(results1,file=paste0("model_",tab$species[i],"_chain_",tab$chain[i],".RData"))