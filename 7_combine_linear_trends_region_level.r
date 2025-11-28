############# CALCULATE AVERAGE TREND PER REGION:
pkgs <- c("data.table", "dplyr","lme4","ggplot2","ggridges","metafor","cowplot","emmeans","tidyverse", "metafor") 
inst <- pkgs %in% installed.packages()
if (any(inst)) install.packages(pkgs[!inst])
pkg_out <- lapply(pkgs, require, character.only = TRUE)

#project_folder="C:/Users/Duchenne/Documents/safeguard/"
project_folder <- ""

trendsf=fread(paste0(project_folder,"data/final_and_intermediate_outputs/all_trends.csv"))

bidon=subset(trendsf,baseline>=1921)
nb_weird=nrow(subset(bidon,abs(trend)>=1))
bidon=subset(bidon,abs(trend)<1)
baselines_vec=unique(trendsf$baseline)

for(jj in unique(bidon$taxo_group)){
	for(j in 1:length(baselines_vec)){
		bidon2=subset(bidon,convergence==0 & !is.na(acim) & baseline==baselines_vec[j] & !is.na(sde) & taxo_group==jj)

		model=rma(trend ~ region_50, sei=sde,data=bidon2, digits=10)
		model2=rma(trend ~ region_50, sei=sde,data=subset(bidon2,max.occ>1e-4), digits=10)
		lis_bas=list(model,model2)
		save(lis_bas,file=paste0(project_folder,"data/final_and_intermediate_outputs/models/model_",baselines_vec[j],"_",jj,".RData"))
	}
}


for(jj in unique(bidon$taxo_group)){
	for(j in 1:length(baselines_vec)){
		bidon2=subset(bidon,convergence==0 & !is.na(acim) & baseline==baselines_vec[j] & !is.na(sde) & taxo_group==jj)
		model=rma(trend ~ 1, sei=sde,data=bidon2, digits=10)
		model2=rma(trend ~ 1, sei=sde,data=subset(bidon2,max.occ>1e-4), digits=10)
		lis_bas=list(model,model2)
		save(lis_bas,file=paste0(project_folder,"data/final_and_intermediate_outputs/models/model_total_",baselines_vec[j],"_",jj,".RData"))
	}
}
