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
	Y[i] ~ dbern(Py[i])
	Py[i]<- z[i]*p[i]
	logit(p[i]) <- country_period_det[country_num[i],period_num[i]]+alpha*log(list_length[i]) 
	
	### State model
	z[i] ~ dbern(muZ[i]) 
	logit(muZ[i]) <- country_period_occ[country_num[i],period_num[i]]+site_eff[site_num[i]] 
}

# PRIORS
# Observation and State models priors
for(j in 1:Nc){
country_period_det[j,1] <- dnorm(0,0.5)
country_period_occ[j,1] <- dnorm(0,0.5)

for(jj in 2:Nperiod){
country_period_det[j,jj] ~ dnorm(country_period_det[j,jj-1],tau.det)
country_period_occ[j,jj] ~ dnorm(country_period_occ[j,jj-1],tau.occ)
}
}

alpha ~ dnorm(0,0.5)

for(j in 1:Nsite){
site_eff[j] ~ dnorm(0,tau.site)
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
