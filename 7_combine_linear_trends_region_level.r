
############# FIRST ASSEMBLE ALL SPECIES TOGETHER
pkgs <- c("data.table", "dplyr") 
inst <- pkgs %in% installed.packages()
if (any(inst)) install.packages(pkgs[!inst])
pkg_out <- lapply(pkgs, require, character.only = TRUE)

#project_folder="C:/Users/Duchenne/Documents/safeguard/"
project_folder <- ""

setwd(dir=paste0(project_folder,"/predicts_linear"))
lifile=list.files( pattern ="predicts_")
lifile=lifile[grep("non_linear",lifile,invert=TRUE)]

trendsf=NULL
for(i in lifile){
	index=gsub(".RData","",gsub("predicts_","",i))
	load(i)
	trends=lili[[1]]
	trendsf=rbind(trendsf,trends)
}

fwrite(trendsf,paste0(project_folder,"data/results_TMB_trends_with_counts.csv"))

inv.logit=function(x){exp(x)/(1+exp(x))}
logit=function(x){log(x/(1-x))}

tab_spec=fread(paste0(project_folder,"data/species_nb_records.csv"))
tab_spec=subset(tab_spec,nb_detect>=5 & nb_records_tot>=10) %>% group_by(scientificName) %>% mutate(nb_detect_tot=sum(nb_detect),occ_min=min(occupancy_obs))

nrow(trendsf)
trendsf=merge(trendsf,tab_spec,by.x=c("species","region_50","taxo_group"),by.y=c("scientificName","region_50","taxo_group"))
nrow(trendsf)

b=unique(trendsf[,c("species","nb_records_tot","nb_detect_tot","convergence","occ_min")])

trendsf$significant="no"
trendsf$significant[trendsf$trend<0 & (trendsf$trend+1.96*trendsf$sde)<0]="yes"
trendsf$significant[trendsf$trend>0 & (trendsf$trend-1.96*trendsf$sde)>0]="yes"

trendsf$growth.rate=100*(exp(trendsf$trend)-1)
trendsf=subset(trendsf,convergence==0 & !is.na(acim) & !is.na(sde))

fwrite(trendsf,paste0(project_folder,"data/all_trends.csv"))


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


######warning######### rma not working
for(jj in unique(bidon$taxo_group)){
	for(j in 1:length(baselines_vec)){
		bidon2=subset(bidon,convergence==0 & !is.na(acim) & baseline==baselines_vec[j] & !is.na(sde) & taxo_group==jj)

		model=rma(trend ~ region_50, sei=sde,data=bidon2, digits=10)
		model2=rma(trend ~ region_50, sei=sde,data=subset(bidon2,max.occ>1e-4), digits=10)
		lis_bas=list(model,model2)
		save(lis_bas,file=paste0(project_folder,"data/model_",baselines_vec[j],"_",jj,".RData"))
	}
}


for(jj in unique(bidon$taxo_group)){
	for(j in 1:length(baselines_vec)){
		bidon2=subset(bidon,convergence==0 & !is.na(acim) & baseline==baselines_vec[j] & !is.na(sde) & taxo_group==jj)
		model=rma(trend ~ 1, sei=sde,data=bidon2, digits=10)
		model2=rma(trend ~ 1, sei=sde,data=subset(bidon2,max.occ>1e-4), digits=10)
		lis_bas=list(model,model2)
		save(lis_bas,file=paste0(project_folder,"data/model_total_",baselines_vec[j],"_",jj,".RData"))
	}
}
