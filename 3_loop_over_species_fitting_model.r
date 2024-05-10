###########################################
###########################################
#' Check for packages and if necessary install into library 
#+ message = FALSE
rm(list=ls())
pkgs <- c("data.table", "dplyr","R2jags") 

inst <- pkgs %in% installed.packages()
if (any(inst)) install.packages(pkgs[!inst])
pkg_out <- lapply(pkgs, require, character.only = TRUE)

#defining working folder:
setwd(dir="C:/Users/Duchenne/Documents/safeguard/data")

#species index
i=1
index=i+6

#import data:
mat=fread("det_nondet_matrix_species_common.csv")

#det/nondet of the focal species
Y=dat[,index]
Y[Y>1]=1 #if many dets, put one


ParsStage <- c("country_period_det","country_period_det","alpha","site_effect",
"site_eff")

Inits <- function(){list()}

t1=Sys.time()
results1 <- jags(model.file="model.txt", parameters.to.save=ParsStage, n.chains=1, data=dat,n.burnin = 50,  n.iter = 100, n.thin = 2,inits =Inits,jags.seed =2)
t2=Sys.time()
t2-t1

save(results1,file=paste0("chain_model_ZI",j,".RData"))


ParsStage <- c("barrier_infer","match_infer","Interceptpz","traitBarrier","Intercept","traitMismatch","pheno","abond","plant_effect","site_effect",
"temp_effect","sitebird_effect","sd.plant","sd.bird","sd.site","sd.temp","r","edec","samp")

Inits <- function(){list()}

t1=Sys.time()
results1 <- jags(model.file="model.txt", parameters.to.save=ParsStage, n.chains=1, data=dat,n.burnin = 50,  n.iter = 100, n.thin = 2,inits =Inits,jags.seed =2)
t2=Sys.time()
t2-t1

save(results1,file=paste0("chain_model_ZI",j,".RData"))
#
