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

