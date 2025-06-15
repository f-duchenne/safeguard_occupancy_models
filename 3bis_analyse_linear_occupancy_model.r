###########################################
###########################################
#' Check for packages and if necessary install into library 
#+ message = FALSE
rm(list=ls())
pkgs <- c("data.table","dplyr","ggeffects","glmmTMB","emmeans") 

inst <- pkgs %in% installed.packages()
if (any(inst)) install.packages(pkgs[!inst])
pkg.out <- lapply(pkgs, require, character.only = TRUE)


# Collect command arguments
args <- commandArgs(trailingOnly = TRUE)
args_contents <- strsplit(args, ' ')
# Get first argument
i <- as.numeric(args_contents[[1]])

print(i)

#defining working folder:
setwd(dir="/home/duchenne/safeguard/")
#import data:
if(i<=1364){
	taxo_group="bees"
	if(i<=317){
		dat=fread(paste0("bees_det_nondet_matrix_common.csv"))
	}else{
		dat=fread(paste0("bees_det_nondet_matrix_rare.csv"))
		i=i-317
	}
}else{
  i=i-1364
	taxo_group="hoverflies"
  dat=fread(paste0("hoverflies_det_nondet_matrix.csv"))
  dat$others=NA
}

#combinations
tab1=expand.grid(species=names(dat)[(which(names(dat)=="region_50")+1):(which(names(dat)=="others")-1)])

logit=function(x){log(x/(1-x))}
func=function(x,y){zi_vcov[x,y]}
func2=function(x,y){vc[x,y]}

trendsf=NULL
nonlf=NULL

#combinations
index=names(dat)[names(dat)==tab1$species[i]]

#files=list.files("/home/duchenne/safeguard/results/",pattern="predicts_")

#if(!(paste0("predicts_",tab1$species[i],".RData") %in% files)){
dat=subset(dat,year_grouped>=1971)
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
dat$period.num2=numFactor(dat$period.num)
dat$group <- factor(rep(1,nrow(dat)))
dat$period.num_s=scale(dat$period.num)
dat$log.list.count_s=scale(dat$log.list.count)
dat$log.list.length.c_s=scale(dat$log.list.length.c)

#keep a table for period labels

print(index)
load(paste0("results/","model_non_linear",tab1$species[i],".RData"))

regions=lili2[[length(lili2)]]
trendsf=NULL

modelt=up2date(lili2[[1]])
zi_vcov=vcov(modelt)[[2]]
deg=which.min(lili2[[2]])
jj=deg  

if(length(regions)>1){
b=ggpredict(modelt,c("period.num2[all]","region_50"),type="zi_prob",condition=c(log.list.length.c_s=0,log.list.count_s=0,MONTH_2=NA,site=NA))
b2=as.data.frame(b)
}else{
b=ggpredict(modelt,c("period.num2[all]"),type="zi_prob",condition=c(log.list.length.c_s=0,log.list.count_s=0,MONTH_2=NA,site=NA))
b2=as.data.frame(b)
b2$group=regions
}


b2$year=b2$x*lili2[[6]][1,"std"]+lili2[[6]][1,"moy"]
b2$deg=deg
b2$species=tab1$species[i]
b2$acim=AIC(modelt)
b2$taxo_group=taxo_group

b2$pred=1-b2$predicted
if(length(regions)>1){
b2=b2 %>% group_by(group) %>% mutate(max.occ=max(pred),min.occ=min(pred))
}else{
b2=b2 %>% group_by(species) %>% mutate(max.occ=max(pred),min.occ=min(pred))
}

trendsf=rbind(trendsf,b2)

lili=list(trendsf)

save(lili,file=paste0("results/predicts_non_linear_",tab1$species[i],".RData"))
#}
