###########################################
###########################################
#' Check for packages and if necessary install into library 
#+ message = FALSE
rm(list=ls())
pkgs <- c("data.table", "ggeffects","glmmTMB","emmeans") 

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
setwd(dir=path_data)
#import data:
if(i<=317){
dat=fread(paste0("det_nondet_matrix_species_common_50.csv"))
}else{
dat=fread(paste0("det_nondet_matrix_species_rare.csv"))
i=i-317
}
which(names(dat)=="Andrena vaga")

dat$period.num=as.numeric(as.factor(dat$year_grouped))-1
#combinations
tab1=expand.grid(species=names(dat)[9:(ncol(dat)-1)])

logit=function(x){log(x/(1-x))}
func=function(x,y){zi_vcov[x,y]}

trendsf=NULL
nonlf=NULL

#combinations
index=names(dat)[names(dat)==tab1$species[i]]
dat$Y=dat[,index,with=FALSE]
count.table=dat[,.(n_records=sum(Y)),by=region_50] #count number of records per region

print(index)
setwd(dir="C:/Users/Duchenne/Documents/safeguard/results_glm")
load(paste0("model_",tab1$species[i],".RData"))
modelt=up2date(lili[[3]])
regions=lili[[7]]
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

trends=data.frame(region_50=regions,trend=-1*trend,sde=sqrt(errors),species=tab1$species[i],convergence=modelt$fit$convergence,conv_message=modelt$fit$message)
trends$trend=trends$trend/lili[[6]][1,"std"]
trends$sde=trends$sde/lili[[6]][1,"std"]

model=up2date(lili[[2]])
if(length(regions)>1){
b=ggpredict(model,c("period.num2","region_50"),type="zi_prob")
}else{
b=ggpredict(model,c("period.num2"),type="zi_prob")
}
nonl=as.data.frame(b)
nonl$pred=1-nonl$predicted
nonl$conf.h=1-nonl$conf.low
nonl$conf.l=1-nonl$conf.high
nonl$occ_logit=logit(nonl$pred)
nonl$x=as.numeric(nonl$x)
nonl$convergence=model$fit$convergence

if(length(regions)>1){
modeltbis=lm(occ_logit~x*group,data=nonl,weight=1/std.error)
tt=as.data.frame(emtrends(modeltbis, var="x", specs="group"))
tt$conv2=model$fit$convergence
trends=cbind(trends,tt[,c("x.trend","SE")])
}else{
modeltbis=lm(occ_logit~x,data=nonl,weight=1/std.error)
tt=as.data.frame(emtrends(modeltbis, var="x"))
tt$conv2=model$fit$convergence
trends=cbind(trends,tt[,c("x.trend","SE")])
}

#trends_period=NULL
#periods=data.frame(starting=c(1,51,30,80),ending=c(50,100,80,100))
#for(j in 1:nrow(periods)){
#bidon=subset(nonl,x>=periods$starting[i] & x<=periods$ending[i])
#if(length(regions)>1){
#modeltbis=lm(occ_logit~x*group,data=bidon,weight=1/std.error)
#tt=as.data.frame(emtrends(modeltbis, var="x", specs="group"))
#tt$conv2=model$fit$convergence
#tt$period=paste0(periods$starting[i],"-",periods$ending[i])
#tt$region_50=regions
#tt$species=tab1$species[i]
#}else{
#modeltbis=lm(occ_logit~x,data=bidon,weight=1/std.error)
#tt=as.data.frame(emtrends(modeltbis, var="x"))
#tt$conv2=model$fit$convergence
#tt$period=paste0(periods$starting[i],"-",periods$ending[i])
#tt$region_50=regions
#tt$species=tab1$species[i]
#}
#trends_period=rbind(trends_period,tt)
#}


lili=list(trends,nonl)

save(lili,file=paste0("results/predicts_",tab1$species[i],".RData"))