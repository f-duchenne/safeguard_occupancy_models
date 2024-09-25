#defining working folder:
setwd(dir="C:/Users/Duchenne/Documents/safeguard/data")


#####################################################################dat#################################
#### JAGS model file
######################################################################################################


model_string="
model{

# MODEL
for(i in 1:N){

	### Observation Model
	Y[i] ~ dbern(min(Py[i]+0.00000000000000001,0.9999999999999999))
	Py[i]<- z[i]*p[i]
	logit(p[i]) <- country.period.det[country.num[i],period.num[i]]+alpha*log(list.length[i]) 
	
	### State model
	z[i] ~ dbern(muZ[i]) 
	logit(muZ[i]) <- region.period.occ[region.num[i],period.num[i]]+site.eff[site.num[i],region.num[i]] 
}

# PRIORS
# Observation and State models priors
for(j in 1:Nc){
country.period.det[j,1] ~ dnorm(0,0.5)
for(jj in 2:Nperiod){
country.period.det[j,jj] ~ dnorm(country.period.det[j,jj-1],tau.det)
}
}

for(j in 1:Nr){
region.period.occ[j,1] ~ dnorm(0,0.5)
for(jj in 2:Nperiod){
region.period.occ[j,jj] ~ dnorm(region.period.occ[j,jj-1],tau.occ)
}
}



# EU occupancy
for(jj in 1:Nperiod){
eu_eff[jj]<-mean(region.period.occ[1:Nr,jj])
}

alpha ~ dnorm(0,0.5)

#random site effect, with one variance per region
for(j in 1:Nsite){
site.eff[j] ~ dnorm(0,tau.site[region.site_num[j]])
}


#HYPERPRIORS
tau.det <- 1/(sd.det * sd.det)
sd.det ~ dt(0, 1, 1)T(0,)
tau.occ <- 1/(sd.occ * sd.occ)
sd.occ ~ dt(0, 1, 1)T(0,)
for(j in 1:Nr){
sd.site[j] ~ dt(0, 1, 1)T(0,)
tau.site[j] <- 1/(sd.site[j] * sd.site[j])
}

}
"

writeLines(model_string,con="model.txt")
