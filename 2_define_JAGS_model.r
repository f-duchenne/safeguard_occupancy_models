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
	logit(muZ[i]) <- country.period.occ[country.num[i],period.num[i]]+site.eff[site.num[i]] 
}

# PRIORS
# Observation and State models priors
for(j in 1:Nc){
country.period.det[j,1] ~ dnorm(0,0.5)
country.period.occ[j,1] ~ dnorm(0,0.5)

for(jj in 2:Nperiod){
country.period.det[j,jj] ~ dnorm(country.period.det[j,jj-1],tau.det)
country.period.occ[j,jj] ~ dnorm(country.period.occ[j,jj-1],tau.occ)
}
}

alpha ~ dnorm(0,0.5)

for(j in 1:Nsite){
site.eff[j] ~ dnorm(0,tau.site)
}

#HYPERPRIORS
tau.det <- 1/(sd.det * sd.det)
sd.det ~ dt(0, 1, 1)T(0,)
tau.occ <- 1/(sd.occ * sd.occ)
sd.occ ~ dt(0, 1, 1)T(0,)
tau.site <- 1/(sd.site * sd.site)
sd.site ~ dt(0, 1, 1)T(0,)

}
"

writeLines(model_string,con="model.txt")
