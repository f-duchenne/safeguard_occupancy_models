###########################################
###########################################
#' Check for packages and if necessary install into library 
#+ message = FALSE
rm(list=ls())
pkgs <- c("data.table", "dplyr","glmmTMB") 

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

print(i)

#import data:
#dat=fread("det_nondet_matrix_species_common_50.csv")
if(i<=317){
dat=fread(paste0("det_nondet_matrix_species_common_50.csv"))
}else{
dat=fread(paste0("det_nondet_matrix_species_rare.csv"))
i=i-317
}

#todo=fread("species_to_repeat.csv", sep=",")
#combinations
tab=expand.grid(species=names(dat)[(which(names(dat)=="region_50")+1):(which(names(dat)=="others")-1)])
index=names(dat)[names(dat)==tab$species[i]]

print(index)


files=list.files("/home/duchenne/safeguard/results/")

#if(!(paste0("model_",tab$species[i],".RData") %in% files)){

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

Nperiod=length(unique(dat$period.num))

Nr=length(unique(dat$region_50))

Nsite=length(unique(dat$site))

Nmonth=length(unique(dat$MONTH_2))

regions=levels(factor(dat$region_50))

#keep a table for period labels
period_tab=dat %>% group_by(year_grouped,period.num) %>% summarise(log.list.length.c.moy=mean(log.list.length.c),log.list.count.moy=mean(log.list.count))

scaling=data.frame(varia=c("period.num","log.list.count","log.list.length.c"),moy=c(mean(dat$period.num),mean(dat$log.list.count),mean(dat$log.list.length.c)),std=c(sd(dat$period.num),sd(dat$log.list.count),sd(dat$log.list.length.c)))

optim_vec=c("Nelder-Mead", "BFGS", "CG")

######## LINEAR TREND
lili=list()
baselines_vec=c(1921,1951,1961,1971,1981,1991,2001)
for(j in 1:length(baselines_vec)){

dat2=subset(dat,year_grouped>=baselines_vec[j])

if(length(unique(dat$region_50))>1){
modelt=glmmTMB(Y~log.list.length.c_s+log.list.count_s+(1|region_50/MONTH_2),family=binomial,data=dat2,ziformula=~period.num_s*region_50+(1|site),control=glmmTMBControl(optCtrl =list(iter.max=1e5,eval.max=1e3)))
}else{
modelt=glmmTMB(Y~log.list.length.c_s+log.list.count_s+(1|MONTH_2),family=binomial,data=dat2,ziformula=~period.num_s+(1|site),control=glmmTMBControl(optCtrl =list(iter.max=1e5,eval.max=1e3)))
}

### If it did not converge try to change the solver
b=0
while((modelt$fit$convergence!=0 | is.na(AIC(modelt))) & b<3){
b=b+1
if(length(unique(dat$region_50))>1){
modelt=glmmTMB(Y~log.list.length.c_s+log.list.count_s+(1|region_50/MONTH_2),family=binomial,data=dat2,ziformula=~period.num_s*region_50+(1|site),control=glmmTMBControl(optimizer=optim,optArgs=list(method=optim_vec[b])))
}else{
modelt=glmmTMB(Y~log.list.length.c_s+log.list.count_s+(1|MONTH_2),family=binomial,data=dat2,ziformula=~period.num_s+(1|site),control=glmmTMBControl(optimizer=optim,optArgs=list(method=optim_vec[b])))
}
}
lili[[j]]=modelt
}

#b=ggpredict(model,c("period.num2","region_50"),type="re.zi")
#ggplot(data=b,aes(x=as.numeric(x),y=predicted,color=group,fill=group))+geom_ribbon(aes(ymin=conf.low,ymax=conf.high),alpha=0.2)+geom_line()


lili2=c(lili,list(baselines_vec,period_tab,nsurvey_tot,nsurvey_used,scaling,regions))

save(lili2,file=paste0("results/model_",tab$species[i],".RData"))
#}
