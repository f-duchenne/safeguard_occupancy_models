pkgs <- c("randomForest","data.table", "dplyr", "lubridate","lme4","R2jags","mcmcOutput","mcmcplots","MCMCvis","ggplot2") 
inst <- pkgs %in% installed.packages()
if (any(inst)) install.packages(pkgs[!inst])
pkg_out <- lapply(pkgs, require, character.only = TRUE)
path_data="C:/Users/Duchenne/Documents/safeguard/data/"
setwd(dir="C:/Users/Duchenne/Documents/safeguard/results_models")

#import data:
dat=fread(paste0(path_data,"det_nondet_matrix_species_common_50.csv"))
#combinations
tab1=expand.grid(species=names(dat)[9:(ncol(dat)-1)])

sumaf=NULL
for(i in 1:75){
#combinations
index=names(dat)[names(dat)==tab1$species[i]]

print(index)
load(paste0("model_",tab1$species[i],"_chain_1.RData"))
chain1=results1
load(paste0("model_",tab1$species[i],"_chain_2.RData"))
chain2=results1
load(paste0("model_",tab1$species[i],"_chain_3.RData"))
chain3=results1

obj1=as.mcmc(chain1)
obj2=as.mcmc(chain2)
obj3=as.mcmc(chain3)
#combine chains and get a summary:
mco <- mcmcOutput(mcmc.list(obj1,obj2,obj3))
suma=summary(mco)
suma$varia=rownames(suma)
suma$species=tab1$species[i]
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

load(paste0("data_",tab1$species[i],".RData"))
suma=merge(suma,lili[[2]],by="period.num",all.x=TRUE,all.y=FALSE)
suma=merge(suma,lili[[3]],by="region.num",all.x=TRUE,all.y=FALSE)

sumaf=rbind(sumaf,suma)

}

#import data:
dat2=fread(paste0(path_data,"det_nondet_matrix_species_rare_50.csv"))
#combinations
tab2=expand.grid(species=names(dat2)[8:(ncol(dat2)-1)])
for(i in 1:length(unique(tab2$species))){
#combinations
index=names(dat2)[names(dat2)==tab2$species[i]]

print(index)
load(paste0("model_",tab2$species[i],"_chain_1.RData"))
chain1=results1
load(paste0("model_",tab2$species[i],"_chain_2.RData"))
chain2=results1
load(paste0("model_",tab2$species[i],"_chain_3.RData"))
chain3=results1

obj1=as.mcmc(chain1)
obj2=as.mcmc(chain2)
obj3=as.mcmc(chain3)
#combine chains and get a summary:
mco <- mcmcOutput(mcmc.list(obj1,obj2,obj3))
suma=summary(mco)
suma$varia=rownames(suma)
suma$species=tab2$species[i]
suma$type=suma$vari
suma$type[grep("det",suma$vari)]="detection"
suma$type[grep("period.occ",suma$vari)]="occupancy"
suma$type[grep("eu_eff",suma$vari)]="eu_eff"
suma$region.num=sapply(strsplit(gsub("]","",sapply(strsplit(suma$vari,"[",fixed=T),function(x){x[2]}),fixed=T),","),function(x){x[1]})
suma$period.num=sapply(strsplit(gsub("]","",sapply(strsplit(suma$vari,"[",fixed=T),function(x){x[2]}),fixed=T),","),function(x){x[2]})

load(paste0("data_",tab2$species[i],".RData"))
suma=merge(suma,lili[[2]],by="period.num")
suma=merge(suma,lili[[3]],by="region.num")

sumaf=rbind(sumaf,suma)


}
sumaf$period.num[sumaf$type=="eu_eff"]=sumaf$region.num[sumaf$type=="eu_eff"]
sumaf$region.num[sumaf$type=="eu_eff"]=NA

fwrite(sumaf,paste0(path_data,"results_bayes.csv"))

#################################################################
pkgs <- c("randomForest","data.table", "dplyr", "lubridate","lme4","R2jags","mcmcOutput","mcmcplots","MCMCvis","ggplot2") 
inst <- pkgs %in% installed.packages()
if (any(inst)) install.packages(pkgs[!inst])
pkg_out <- lapply(pkgs, require, character.only = TRUE)
path_data="C:/Users/Duchenne/Documents/safeguard/data/"

inv.logit=function(x){exp(x)/(1+exp(x))}

sumaf=fread(paste0(path_data,"results_bayes.csv"))

#CONVERGENCE CHEKS:
b= sumaf %>% group_by(species,region_50) %>% summarise(prop_bad_rhat=sum(Rhat>1.1)/length(Rhat),prop_bad_MCE=sum(MCEpc>5)/length(MCEpc))

ggplot(data=subset(sumaf,type=="occupancy" & species=="Andrena agilissima"),aes(x=as.numeric(period.num),y=inv.logit(mean),color=species))+geom_point(alpha=0.7)+xlab("Periode")+theme_bw()+ylab("detection")+
geom_line()+
theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.border=element_blank(),
panel.grid.minor = element_blank(),panel.background = element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),
strip.background=element_rect(fill=NA,color=NA),legend.position="none")+facet_wrap(~region_50,scales="free")

#REGIONAL INDEXES:
sumaf$estimate=sumaf$mean
b2=subset(sumaf,type=="occupancy") %>% group_by(region_50,year_grouped,period.num) %>% summarise(global_index=mean(inv.logit(estimate)),global_dev=sd(estimate)/sqrt(length(estimate)))


ggplot(data=b2,aes(x=as.numeric(year_grouped),y=global_index,color=region_50))+geom_point(alpha=0.7)+xlab("Periode")+theme_bw()+ylab("occupancy")+
geom_line()+
theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.border=element_blank(),
panel.grid.minor = element_blank(),panel.background = element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),
strip.background=element_rect(fill=NA,color=NA))

b3=subset(sumaf,type=="detection") %>% group_by(region_50,year_grouped,period.num) %>% summarise(global_index=mean(inv.logit(estimate)),global_dev=sd(estimate)/sqrt(length(estimate)))

ggplot(data=b3,aes(x=as.numeric(year_grouped),y=global_index,color=region_50))+geom_point(alpha=0.7)+xlab("Periode")+theme_bw()+ylab("detection")+
geom_line()+
theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.border=element_blank(),
panel.grid.minor = element_blank(),panel.background = element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),
strip.background=element_rect(fill=NA,color=NA))



#### SPECIES trends
trends=NULL
for(i in unique(b$species)){
bidon=subset(sumaf2,type=="occupancy" & species==i)
for(z in unique(bidon$region_50)){
bidon2=subset(bidon,region_50==z)
model=lm(mean~as.numeric(as.character(period.num)),data=bidon2,weight=1/sd)
suma=summary(model)
trends=rbind(trends,data.frame(species=i,trend=suma$coefficients[2,1],trend_se=suma$coefficients[2,2],region_50=z))
}
}

ggplot(data=trends,aes(x=trend,color=region_50,fill=region_50))+geom_density(,alpha=0.4)+geom_vline(xintercept=0,linetype="dashed")+
theme_bw()+
theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.border=element_blank(),
panel.grid.minor = element_blank(),panel.background = element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),
strip.background=element_rect(fill=NA,color=NA))

ggplot(data=trends,aes(y=trend,x=region_50,color=region_50,fill=region_50))+ stat_summary(fun.data = "mean_cl_boot",linewidth = 1, size = 1)+geom_hline(yintercept=0,linetype="dashed")+
theme_bw()+
theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.border=element_blank(),
panel.grid.minor = element_blank(),panel.background = element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),
strip.background=element_rect(fill=NA,color=NA))

#### BELGIUM trendsliste2=read.table("yearly_estimates_of_occupancy.txt",sep="\t",header=T)
liste=fread(paste0(path_data,"yearly_estimates_of_occupancy_BELGIUM.txt"))
belg=NULL
for(i in unique(liste$species)){
bidon2=subset(liste,species==i)
bidon2$Annee=bidon2$Annee-1950
se=logit(bidon2$mean)-logit(bidon2$quant_025)
model=lm(logit(mean)~Annee,data=bidon2,weights=1/se)
belg=rbind(belg,data.frame(trend_effect=model$coeff["Annee"],
trend_err=summary(model)$coefficient["Annee","Std. Error"],
trend_pval=car::Anova(model)["Annee",4],
prob50=bidon2[bidon2$Annee==0,"mean"],species=i))
}

bibi=merge(subset(trends,region_50=="atlantic"),belg,by="species")

ggplot(data=bibi,aes(x=trend_effect,y=trend))+geom_point()+xlab("Trends in Belgium")+
ylab("Trend from safeguard in atlantic region")

subset(bibi,species=="Bombus terrestris")


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




sp_to_test=c("Andrena agilissima","Andrena strohmella","Halictus scabiosae","Lasioglossum minutulum","Bombus terrestris")
sumaf=NULL
for(gr in 1:2){
for(i in sp_to_test){
#combinations

load(paste0("model_",i,"_chain_1_group_year_",gr,".RData"))
chain1=results1
load(paste0("model_",i,"_chain_2_group_year_",gr,".RData"))
chain2=results1
load(paste0("model_",i,"_chain_3_group_year_",gr,".RData"))
chain3=results1

obj1=as.mcmc(chain1)
obj2=as.mcmc(chain2)
obj3=as.mcmc(chain3)
#combine chains and get a summary:
mco <- mcmcOutput(mcmc.list(obj1,obj2,obj3))
suma=summary(mco)
suma$varia=rownames(suma)
suma$species=i
suma$type=suma$vari
suma$type[grep("det",suma$vari)]="detection"
suma$type[grep("period.occ",suma$vari)]="occupancy"
suma$type[grep("eu_eff",suma$vari)]="eu_eff"
suma$region.num=sapply(strsplit(gsub("]","",sapply(strsplit(suma$vari,"[",fixed=T),function(x){x[2]}),fixed=T),","),function(x){x[1]})
suma$period.num=sapply(strsplit(gsub("]","",sapply(strsplit(suma$vari,"[",fixed=T),function(x){x[2]}),fixed=T),","),function(x){x[2]})
suma$gr=gr

load(paste0("data_",i,".RData"))
suma=merge(suma,lili[[3]],by="region.num")

sumaf=rbind(sumaf,suma)
}}

b= sumaf %>% group_by(species,region_50,gr) %>% summarise(prop_bad_rhat=sum(Rhat>1.1)/length(Rhat),prop_bad_MCE=sum(MCEpc>5)/length(MCEpc))

ggplot(data=subset(sumaf,type=="detection"),aes(x=as.numeric(period.num),y=inv.logit(mean),color=region_50))+geom_point(alpha=0.7)+xlab("Periode")+theme_bw()+ylab("detection")+
geom_line()+
theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.border=element_blank(),
panel.grid.minor = element_blank(),panel.background = element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),
strip.background=element_rect(fill=NA,color=NA))+facet_wrap(~region_50,scales="free")+facet_wrap(~species+gr)