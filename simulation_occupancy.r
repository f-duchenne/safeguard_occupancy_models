###########################################
#' Check for packages and if necessary install into library 
#+ message = FALSE
rm(list=ls())
pkgs <- c("data.table", "dplyr","sparta","R2jags","ggplot2","glmmTMB") 
library(ggeffects)

inst <- pkgs %in% installed.packages()
if (any(inst)) install.packages(pkgs[!inst])
pkg_out <- lapply(pkgs, require, character.only = TRUE)

#defining working folder:
setwd(dir="C:/Users/Duchenne/Documents/safeguard/data")
inv.logit=function(x){exp(x)/(1+exp(x))}

Nperiod=30
Nsite=20
sites=data.frame(site=letters[1:Nsite],region=rep(c(1,2),each=Nsite/2),site_eff=rnorm(Nsite,0,0.5))

years1=data.frame(year=1:Nperiod,detect=rnorm(Nperiod,-2.5,0.6),occupancy=sapply(seq(-2.5,-0.5,length.out=Nperiod),function(x){rnorm(1,x,0.001)}),region=1)
years2=data.frame(year=1:Nperiod,detect=years1$detect,occupancy=sapply(seq(0.5,-2.5,length.out=Nperiod),function(x){rnorm(1,x,0.001)}),region=2)
years=rbind(years1,years2)
years$detect[years$year>15]=years$detect[years$year>15]+1.5

survey=expand.grid(site=unique(sites$site),month=1:12,year=1:Nperiod)
survey$prob=sqrt(survey$year)/sum(sqrt(survey$year))

survey2=survey[sample(1:nrow(survey),round(nrow(survey)/7),replace=FALSE,prob=survey$prob),]
survey2$list.length=1+rnbinom(nrow(survey2),size=1,prob=0.3)
survey2$nb.data=floor(runif(rep(1,nrow(survey2)),survey2$list.length,survey2$list.length*survey2$year))

nb.occ.effect=1.1
list.length.effect=0.9

occ_tab=merge(sites,years,by=c("region"))
occ_tab$state=sapply(inv.logit(occ_tab$occupancy+occ_tab$site_eff),function(x){rbinom(1,1,x)})

dat=merge(survey2,occ_tab,by=c("site","year"))
dat$log.list.length=log(dat$list.length)
dat=dat %>% group_by(region) %>% mutate(log.list.length.c=log.list.length-mean(log.list.length))

#list count standardization
dat$log.list.count=log(dat$nb.data)

prob_vec=dat$state*inv.logit(0.5+dat$detect-(dat$month-6)^2/10+list.length.effect*dat$log.list.length+nb.occ.effect*dat$log.list.count)
dat$Y=sapply(prob_vec,function(x){rbinom(1,1,x)})
dat$occ.proba=inv.logit(dat$occupancy+dat$site_eff)

dat$year_grouped=plyr::round_any(dat$year,5,f=ceiling)
dat$period.num=as.numeric(as.factor(dat$year_grouped))

dat$month.num=dat$month

dat$region.num=dat$region

dat$site.num=as.numeric(as.factor(dat$site))

Nperiod=length(unique(dat$period.num))

Nr=length(unique(dat$region.num))

Nsite=length(unique(dat$site))

Nmonth=length(unique(dat$month))
Y=dat$Y
N=nrow(dat)

#keep a table for region-site labels
region.site=unique(dat[,c("site.num","region.num")])
region.site=region.site[order(region.site$site.num),]

# create list of data for JAGS
dat.model=list(Y=Y,N=N,Nmonth=Nmonth,Nperiod=Nperiod,Nsite=Nsite,Nr=Nr,year_grouped=dat$year_grouped,period.num=dat$period.num,region=dat$region,region.num=dat$region.num,log.list.length.c=dat$log.list.length.c,month.num=dat$month.num,site=dat$site,site.num=dat$site.num,log.list.count=dat$log.list.count,region.site_num=region.site$region.num)


ggplot()+
stat_summary(data=years,aes(x=year,y=inv.logit(occupancy)),fun.y = "mean",geom="line",linewidth = 1, size = 1,color="blue")+
stat_summary(data=dat,aes(x=year,y=occ.proba),fun.y = "mean",geom="line",linewidth = 1, size = 1,color="orange")+
stat_summary(data=dat,aes(x=year,y=Y),fun.y = "mean",geom="line",linewidth = 1, size = 1,color="grey")+
stat_summary(data=occ_tab,aes(x=year,y=state),fun.y = "mean",geom="line",linewidth = 1, size = 1,color="red")+
theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.border=element_blank(),
panel.grid.minor = element_blank(),panel.background = element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),
strip.background=element_rect(fill=NA,color=NA),legend.position="none")+facet_wrap(~region,scales="free")


model_string="
model{

### State model
for (i in 1:Nsite){
for (t in 1:Nperiod){
  z[i,t] ~ dbern(muZ[i,t]) 
  logit(muZ[i,t]) <- region.period.occ[region.site_num[i],t]+site.eff[i]
}
}

### Observation Model
for(j in 1:N){
  Y[j] ~ dbern(min(Py[j]+0.00000000000000001,0.9999999999999999))
  Py[j]<- p[j]*z[site.num[j],period.num[j]]
  logit(p[j]) <- period.det[period.num[j]]+reg.month[region.num[j],month.num[j]]+alpha*log.list.length.c[j]+beta*log.list.count[j]
}

# PRIORS
## State models priors
for(j in 1:Nr){
  sd.occ[j] ~ dt(0, 1, 1)T(0,) 
  tau.occ[j] <- 1/(sd.occ[j] * sd.occ[j])
}

for(j in 1:Nr){
  region.period.occ[j,1] ~ dnorm(0,0.5) 
  for(jj in 2:Nperiod){
    region.period.occ[j,jj] ~ dnorm(region.period.occ[j,jj-1],tau.occ[j])
  }
}

### random site effect
for(j in 1:Nsite){
site.eff[j] ~ dnorm(0,tau.site)
}

### EU occupancy
for(jj in 1:Nperiod){
eu_eff[jj]<-mean(z[1:Nsite,jj])
}

## Observation model priors
### temporal variation of detection probability
period.det.mu ~ dnorm(0,0.5) 
for(jj in 1:Nperiod){
   period.det[jj] ~ dnorm(period.det.mu,tau.det)
}


### month:region random effect
for(j in 1:Nr){
  mu_region[j] ~ dnorm(0,tau.reg)
  for(jjj in 1:Nmonth){
    reg.month[j,jjj] ~ dnorm(mu_region[j],tau.month)
  }
}

### sampling effects
alpha ~ dnorm(0,0.01)
beta ~ dnorm(0,0.01)


#HYPERPRIORS
tau.det <- 1/(sd.det * sd.det)
sd.det ~ dt(0, 1, 1)T(0,)
tau.site <- 1/(sd.site * sd.site)
sd.site ~ dt(0, 1, 1)T(0,)
tau.month <- 1/(sd.month * sd.month)
sd.month ~ dt(0, 1, 1)T(0,)
tau.reg <- 1/(sd.reg * sd.reg)
sd.reg ~ dt(0, 1, 1)T(0,)

}

"

writeLines(model_string,con="model.txt")


ParsStage <- c("period.det","region.period.occ","alpha","beta","eu_eff","sd.occ","sd.det","sd.site","sd.month","reg.month")

Inits <- function(){list()}

t1=Sys.time()
results1 <- jags(model.file="model.txt", parameters.to.save=ParsStage, n.chains=3, data=dat.model,n.burnin = 20000,  n.iter = 30000, n.thin = 3,inits =Inits,jags.seed=i)
t2=Sys.time()
t2-t1


# fit1 <- stan(
  # file = "model.stan",  # Stan program
  # data = dat.model,    # named list of data
  # chains = 3,             # number of Markov chains
  # warmup = 20000,          # number of warmup iterations per chain
  # iter = 30000,            # total number of iterations per chain
  # cores = 1,              # number of cores (could use one per chain)
  # refresh = 0             # no progress shown
  # )


suma=as.data.frame(results1$BUGSoutput$summary)
suma$varia=rownames(suma)
suma$type=suma$vari
suma$type[grep("period.det",suma$vari)]="detection"
suma$type[grep("period.occ",suma$vari)]="occupancy"
suma$type[grep("eu_eff",suma$vari)]="eu_eff"
first=sapply(strsplit(gsub("]","",sapply(strsplit(suma$vari,"[",fixed=T),function(x){x[2]}),fixed=T),","),function(x){x[1]})
second=sapply(strsplit(gsub("]","",sapply(strsplit(suma$vari,"[",fixed=T),function(x){x[2]}),fixed=T),","),function(x){x[2]})
suma$region.num=NA
suma$region.num[grep("reg",suma$vari)]=first[grep("reg",suma$vari)]
suma$period.num=NA
suma$period.num[grep("reg",suma$vari)]=second[grep("reg",suma$vari)]
suma$period.num[suma$type %in% c("detection","eu_eff")]=first[suma$type %in% c("detection","eu_eff")]
suma$region=suma$region.num

subset(suma,varia=="beta")
subset(suma,varia=="alpha")


modelg=glmmTMB(Y~log.list.length.c+log.list.count+(1|year)+(1|region/month.num),family=binomial,data=dat,ziformula=~year*as.factor(region)+(1|site))

summary(modelg)
b=ggpredict(modelg,c("year","region"))
b$region=b$group
ggplot()+
geom_point(data=subset(suma,type=="occupancy"),aes(x=as.numeric(period.num)*5,y=inv.logit(mean)),alpha=0.7)+
geom_line(data=subset(suma,type=="occupancy"),aes(x=as.numeric(period.num)*5,y=inv.logit(mean)))+
stat_summary(data=dat,aes(x=year,y=Y),fun.y = "mean",geom="line",linewidth = 1, size = 1,color="grey")+
stat_summary(data=years,aes(x=year,y=inv.logit(occupancy)),fun.y = "mean",geom="line",linewidth = 1, size = 1,color="blue")+
stat_summary(data=dat,aes(x=year,y=occ.proba),fun.y = "mean",geom="line",linewidth = 1, size = 1,color="orange")+
geom_line(data=b,aes(x=x,y=predicted),fun.y = "mean",geom="line",linewidth = 1, size = 1,color="green")+
theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.border=element_blank(),
panel.grid.minor = element_blank(),panel.background = element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),
strip.background=element_rect(fill=NA,color=NA),legend.position="none")+facet_wrap(~region,scales="free")

ggplot()+
geom_point(data=subset(suma,type=="detection"),aes(x=as.numeric(period.num),y=inv.logit(mean)),alpha=0.7)+
geom_line(data=subset(suma,type=="detection"),aes(x=as.numeric(period.num),y=inv.logit(mean)))