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

i=1
print(i)

#defining working folder:
setwd(dir="/data/")
setwd(dir="C:/Users/Duchenne/Documents/safeguard/data")

#import data:
dat=fread("det_nondet_matrix_species_common_50.csv")

#combinations
tab=expand.grid(chain=1:3,species=names(dat)[(which(names(dat)=="region_50")+1):(which(names(dat)=="others")-1)])
index=names(dat)[names(dat)==tab$species[i]]

print(index)

dat$Y=dat[,index,with=F]
dat$Y[dat$Y>1]=1 #if many dets, put one
#det/nondet of the focal species
Y=dat$Y
count.table=dat[,.(n_records=sum(Y)),by=region_50] #count number of records per region


nsurvey_tot=nrow(dat)
#exclude region with less than 5 record of the focal species
dat=subset(dat,region_50 %in% subset(count.table,n_records>=5)$region_50)
nsurvey_used=nrow(dat)

#list length standardization
dat$log.list.length=log(dat$list_length)
dat=dat %>% group_by(region_50) %>% mutate(log.list.length.c=log.list.length-mean(log.list.length))

#list count standardization
dat$log.list.count=log(dat$record_number)

N=nrow(dat)

print(N)

dat$period.num=as.numeric(as.factor(dat$year_grouped))

dat$month.num=dat$MONTH_2

dat$region.num=as.numeric(as.factor(dat$region_50))

dat$site.num=as.numeric(as.factor(dat$site))

Nperiod=length(unique(dat$period.num))

Nr=length(unique(dat$region_50))

Nsite=length(unique(dat$site))

Nmonth=length(unique(dat$MONTH_2))

#keep a table for region-site labels
region.site=unique(dat[,c("site.num","region.num")])
region.site=region.site[order(region.site$site.num),]
#keep a table for period labels
period_tab=dat %>% group_by(year_grouped,period.num) %>% summarise(log.list.length.c.moy=mean(log.list.length.c),log.list.count.moy=mean(log.list.count))
#keep a table for region labels
region_tab=unique(data.frame(region_50=dat$region_50,region.num=dat$region.num))
region_tab=merge(region_tab,count.table,by="region_50")

# create list of data for JAGS
dat.model=list(Y=Y,N=N,Nmonth=Nmonth,Nperiod=Nperiod,Nsite=Nsite,Nr=Nr,year_grouped=dat$year_grouped,period.num=dat$period.num,region_50=dat$region_50,region.num=dat$region.num,log.list.length.c=dat$log.list.length.c,month.num=dat$month.num,site=dat$site,site.num=dat$site.num,region.site_num=region.site$region.num,log.list.count=dat$log.list.count)

print(N)


#save tables for labels
lili=list(region.site,period_tab,region_tab,nsurvey_tot,nsurvey_used)
if(tab$chain[i]==1){
  save(lili,file=paste0("data_",tab$species[i],".RData"))
}


rm(dat)

ParsStage <- c("period.det","region.period.occ","alpha","beta","eu_eff","sd.occ","sd.det","sd.site","sd.month","reg.month")

Inits <- function(){list()}

t1=Sys.time()
results1 <- jags(model.file="model.txt", parameters.to.save=ParsStage, n.chains=1, data=dat.model,n.burnin = 110,  n.iter = 120, n.thin = 3,inits =Inits,jags.seed=i)
t2=Sys.time()
t2-t1


save(results1,file=paste0("model_",tab$species[i],"_chain_",tab$chain[i],".RData"))