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

#defining working folder:
setwd(dir="/home/duchenne/safeguard/")
#import data:
if(i<=317){
dat=fread(paste0("det_nondet_matrix_species_common_50.csv"))
}else{
dat=fread(paste0("det_nondet_matrix_species_rare.csv"))
i=i-317
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

dat$period.num=as.numeric(as.factor(dat$year_grouped))
dat$period.num2=factor(dat$period.num,levels=1:100)
dat$group <- factor(rep(1,nrow(dat)))
dat$period.num_s=scale(dat$period.num)
dat$log.list.count_s=scale(dat$log.list.count)
dat$log.list.length.c_s=scale(dat$log.list.length.c)

print(index)
load(paste0("results/","model_",tab1$species[i],".RData"))

baselines_vec=lili2[[8]]
regions=lili2[[length(lili2)]]
trendsf=NULL

for(j in 1:length(baselines_vec)){
  dat2=subset(dat,year_grouped>=baselines_vec[j])
  modelt=up2date(lili2[[j]])
  zi_vcov=vcov(modelt)[[2]]
  trend=fixef(modelt)[[2]]["period.num_s"]+c(0,fixef(modelt)[[2]][grep("period.num_s:",names(fixef(modelt)[[2]]))])
  
  if(length(regions)>1){
  trend=fixef(modelt)[[2]]["period.num_s"]+c(0,fixef(modelt)[[2]][grep("period.num_s:",names(fixef(modelt)[[2]]))])
  errors=zi_vcov["zi~period.num_s","zi~period.num_s"]+
  c(0,mapply(func,grep("period.num_s:",colnames(zi_vcov)),grep("period.num_s:",colnames(zi_vcov)))+
  2*mapply(func,grep("period.num_s:",colnames(zi_vcov)),rep(which(colnames(zi_vcov)=="zi~period.num_s"),length(regions)-1)))
  }else{
  trend=fixef(modelt)[[2]]["period.num_s"]
  errors=zi_vcov["zi~period.num_s","zi~period.num_s"]
  }

  trends=data.frame(region_50=regions,trend=-1*trend,sde=sqrt(errors),species=tab1$species[i],convergence=modelt$fit$convergence,acim=AIC(modelt))
  trends$trend=trends$trend/lili2[[12]][1,"std"]
  trends$sde=trends$sde/lili2[[12]][1,"std"]
  trends$baseline=j
  
  trendsf=rbind(trendsf,trends)
}

lili=list(trendsf)

save(lili,file=paste0("results/predicts_",tab1$species[i],".RData"))
