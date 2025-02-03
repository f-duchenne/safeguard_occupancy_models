###########################################
###########################################
#' Check for packages and if necessary install into library 
#+ message = FALSE
rm(list=ls())
pkgs <- c("data.table", "dplyr","R2jags") 

inst <- pkgs %in% installed.packages()
if (any(inst)) install.packages(pkgs[!inst])
pkg.out <- lapply(pkgs, require, character.only = TRUE)

#defining working folder:
setwd(dir="C:/Users/Duchenne/Documents/safeguard/data")

# Collect command arguments
# Capture the arguments passed from the command line
#import data:
dat=fread("det_nondet_matrix_species_test.csv")
sp_to_test=c("Andrena agilissima","Andrena strohmella","Halictus scabiosae","Lasioglossum minutulum")

spi=1

i=which(names(dat)==sp_to_test[spi])
index=i

print(index)

dat=subset(dat,region_50=="atlantic")

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

dat$time_period.num=as.numeric(as.factor(dat$YEAR_2))

dat$month.num=dat$MONTH_2

country.num=as.numeric(as.factor(dat$COUNTRY))

region.num=as.numeric(as.factor(dat$region_50))

site.num=as.numeric(as.factor(dat$site))

Nc=length(unique(dat$COUNTRY))

Nperiod=length(unique(dat$time_period.num))

Nr=length(unique(dat$region_50))

Nsite=length(unique(dat$site))

Nmonth=length(unique(dat$MONTH_2))


model1=glmmTMB(Y ~ time_period.num+(1|site)+(1|MONTH_2)+log.list.length.c,data=dat,family="binomial")


#import data:
dat2=fread("det_nondet_matrix_species_common_50.csv")
index2=which(names(dat2)==sp_to_test[spi])
dat2$Y=dat2[,index2,with=F]
dat2$Y[dat2$Y>1]=1 #if many dets, put one
#det/nondet of the focal species
dat2=subset(dat2,region_50=="atlantic")

dat2$log.list.length=log(dat2$list_length)
dat2=dat2 %>% group_by(region_50) %>% mutate(log.list.length.c=log.list.length-mean(log.list.length))

dat2$time_period.num=1
dat2$time_period.num[dat2$time_period=="1941-1960"]=2
dat2$time_period.num[dat2$time_period=="1961-1980"]=3
dat2$time_period.num[dat2$time_period=="1981-2000"]=4
dat2$time_period.num[dat2$time_period=="2001-2020"]=5

model2=glmmTMB(Y ~ time_period.num+(1|site)+log.list.length.c,data=dat2,family="binomial")


region.site=unique(as.data.table(cbind(site.num,region.num)))
region.site=region.site[order(region.site$site.num),]

dat.model=list(Y=Y,N=N,Nc=Nc,Nmonth=Nmonth,Nperiod=Nperiod,Nsite=Nsite,Nr=Nr,time_period=dat$time_period,period.num=dat$time_period.num,COUNTRY=dat$COUNTRY,region_50=dat$region_50,region.num=region.num,country.num=country.num,log.list.length.c=dat$log.list.length.c,month.num=dat$month.num,site=dat$site,site.num=site.num,region.site_num=region.site$region.num)

print(N)

country_tab=unique(data.frame(COUNTRY=dat.model$COUNTRY,country.num=country.num))
period_tab=unique(data.frame(time_period=dat$YEAR_2,period.num=dat$time_period.num))
region_tab=unique(data.frame(region_50=dat.model$region_50,region.num=region.num))

lili=list(country_tab,period_tab,region_tab,nsurvey_tot,nsurvey_used,region.site)

save(lili,file=paste0("test_data_",sp_to_test[spi],".RData"))


species=names(dat)[index]

rm(dat)

ParsStage <- c("region.period.det","region.period.occ","alpha","eu_eff","sd.occ","sd.det","sd.site")

Inits <- function(){list()}

t1=Sys.time()
results1 <- jags(model.file="model.txt", parameters.to.save=ParsStage, n.chains=3, data=dat.model,n.burnin = 4000,  n.iter = 5000, n.thin = 3,inits =Inits,jags.seed=i)
t2=Sys.time()
t2-t1


save(results1,file=paste0("test_model_",tab$species[i],".RData"))



suma=as.data.frame(results1$BUGSoutput$summary)
suma$varia=rownames(suma)
suma$type=suma$vari
suma$type[grep("det",suma$vari)]="detection"
suma$type[grep("period.occ",suma$vari)]="occupancy"
suma$type[grep("eu_eff",suma$vari)]="eu_eff"
suma$region.num=sapply(strsplit(gsub("]","",sapply(strsplit(suma$vari,"[",fixed=T),function(x){x[2]}),fixed=T),","),function(x){x[1]})
suma$period.num=sapply(strsplit(gsub("]","",sapply(strsplit(suma$vari,"[",fixed=T),function(x){x[2]}),fixed=T),","),function(x){x[2]})

suma=merge(suma,lili[[2]],by="period.num")
suma=merge(suma,lili[[3]],by="region.num")
suma$species=sp_to_test[spi]

suma$period.num=as.numeric(suma$period.num)

plot(inv.logit(mean)~time_period,data=subset(suma,type=="detection" & region_50=="atlantic"),xlab="years",ylab="detection probability")