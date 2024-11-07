#defining working folder:
setwd(dir="C:/Users/Duchenne/Documents/safeguard/data")


#####################################################################dat#################################
#### JAGS model file
######################################################################################################


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
  Py[j]<- z[site.num[j],period.num[j]]*p[j]
  logit(p[j]) <- region.period.det[region.num[j],period.num[j]]+alpha*log.list.length.c[j]
}

# PRIORS
# Observation and State models priors
for(j in 1:Nr){
  region.period.det[j,1] ~ dnorm(0,0.5)
  region.period.occ[j,1] ~ dnorm(0,0.5)
  for(jj in 2:Nperiod){
    region.period.det[j,jj] ~ dnorm(region.period.det[j,jj-1],tau.det)
    region.period.occ[j,jj] ~ dnorm(region.period.occ[j,jj-1],tau.occ)
  }
}

# EU occupancy
for(jj in 1:Nperiod){
	eu_eff[jj]<-mean(z[1:Nsite,jj])
}

alpha ~ dnorm(0,0.5)

#random site effect, with one variance per region
for(j in 1:Nsite){
site.eff[j] ~ dnorm(0,tau.site)
}

#HYPERPRIORS
tau.det <- 1/(sd.det * sd.det)
sd.det ~ dt(0, 1, 1)T(0.00001,20)
tau.occ <- 1/(sd.occ * sd.occ)
sd.occ ~ dt(0, 1, 1)T(0.00001,20)
tau.site <- 1/(sd.site * sd.site)
sd.site ~ dt(0, 1, 1)T(0.00001,20)

}
"

writeLines(model_string,con="model.txt")

