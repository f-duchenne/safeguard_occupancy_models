pkgs <- c("randomForest","data.table", "dplyr", "lubridate","lme4","R2jags","mcmcOutput","mcmcplots","MCMCvis","ggplot2") 
inst <- pkgs %in% installed.packages()
if (any(inst)) install.packages(pkgs[!inst])
pkg_out <- lapply(pkgs, require, character.only = TRUE)
setwd(dir="C:/Users/Duchenne/Documents/safeguard/data")
#import data:
dat=fread("det_nondet_matrix_species_common_50.csv")

#combinations
tab=expand.grid(species=names(dat)[8:(ncol(dat)-1)])

setwd(dir="C:/Users/Duchenne/Documents/safeguard/results_models")

sumaf=NULL
for(i in 1:46){
#combinations
index=names(dat)[names(dat)==tab$species[i]]

print(index)
vec=c("Andrena aeneiventris")
load(paste0("model_",tab$species[i],"_chain_1.RData"))
chain1=results1
load(paste0("model_",tab$species[i],"_chain_2.RData"))
chain2=results1
load(paste0("model_",tab$species[i],"_chain_3.RData"))
chain3=results1

obj1=as.mcmc(chain1)
obj2=as.mcmc(chain2)
obj3=as.mcmc(chain3)
#combine chains and get a summary:
mco <- mcmcOutput(mcmc.list(obj1,obj2,obj3))
suma=summary(mco)
suma$varia=rownames(suma)
suma$species=tab$species[i]
sumaf=rbind(sumaf,suma)
}

sumaf$type=sumaf$vari
sumaf$type[grep("det",sumaf$vari)]="detection"
sumaf$type[grep("period.occ",sumaf$vari)]="occupancy"
sumaf$type[grep("eu_eff",sumaf$vari)]="eu_eff"
sumaf$country.num=sapply(strsplit(gsub("]","",sapply(strsplit(sumaf$vari,"[",fixed=T),function(x){x[2]}),fixed=T),","),function(x){x[1]})
sumaf$period.num=sapply(strsplit(gsub("]","",sapply(strsplit(sumaf$vari,"[",fixed=T),function(x){x[2]}),fixed=T),","),function(x){x[2]})
sumaf$period.num[sumaf$type=="eu_eff"]=sumaf$country.num[sumaf$type=="eu_eff"]
sumaf$country.num[sumaf$type=="eu_eff"]=NA

b= sumaf %>% group_by(species) %>% summarise(prop_bad_rhat=sum(Rhat>1.1)/length(Rhat),prop_bad_MCE=sum(MCEpc>5)/length(MCEpc))

ggplot(data=subset(sumaf,type=="eu_eff"),aes(x=period.num,y=mean,group=species,color=species,fill=species))+geom_line()+geom_ribbon(aes(ymin=l95,ymax=u95),alpha=0.2,color=NA)+xlab("Periode")+theme_bw()+ylab("occupancy")+
theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.border=element_blank(),
panel.grid.minor = element_blank(),panel.background = element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),legend.position="none",
strip.background=element_rect(fill=NA,color=NA))+facet_wrap(~species)

ggplot(data=subset(sumaf,type=="alpha"),aes(x=species,y=mean))+geom_pointrange(aes(ymin=l95,ymax=u95),alpha=0.2)+xlab("Periode")+theme_bw()+ylab("occupancy")+
theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.border=element_blank(),
panel.grid.minor = element_blank(),panel.background = element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),legend.position="none",
strip.background=element_rect(fill=NA,color=NA))



dat$Y=dat[,index,with=F]
#det/nondet of the focal species
dat$Y[dat$Y>1]=1 #if many dets, put one
Y=dat$Y
count.table=dat[,sum(Y),by=region_20] #count number of records per country

nsurvey_tot=nrow(dat)
#exclude country without record of the focal species
dat=subset(dat,region_20 %in% subset(count.table,V1>0)$region_20)
nsurvey_used=nrow(dat)

bidon=subset(dat,Y>0)

length(unique(bidon$site))

bidon %>% group_by(region_20) %>% summarise(nsite=length(unique(site)))

pkgs <- c("randomForest","data.table", "dplyr", "lubridate","lme4","R2jags","mcmcOutput","mcmcplots","MCMCvis") 
inst <- pkgs %in% installed.packages()
if (any(inst)) install.packages(pkgs[!inst])
pkg_out <- lapply(pkgs, require, character.only = TRUE)
setwd(dir="C:/Users/Duchenne/Documents/safeguard/results_models")
vec=c("Andrena aeneiventris")
load(paste0("model_",vec[1],"_chain_1.RData"))
chain1=results1
load(paste0("model_",vec[1],"_chain_2.RData"))
chain2=results1
load(paste0("model_",vec[1],"_chain_3.RData"))
chain3=results1

obj1=as.mcmc(chain1)
obj2=as.mcmc(chain2)
obj3=as.mcmc(chain3)
#combine chains and get a summary:
mco <- mcmcOutput(mcmc.list(obj1,obj2,obj3))
suma=summary(mco)
suma$varia=rownames(suma)