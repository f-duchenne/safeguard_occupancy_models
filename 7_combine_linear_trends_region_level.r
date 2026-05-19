############# CALCULATE AVERAGE TREND PER REGION:
pkgs <- c("data.table", "dplyr","lme4","ggplot2","ggridges","metafor","cowplot","emmeans","tidyverse", "metafor") 
pkg_out <- lapply(pkgs, require, character.only = TRUE)

#project_folder="C:/Users/Duchenne/Documents/safeguard/"
project_folder <- ""

trendsf=fread(paste0(project_folder,"data/final_and_intermediate_outputs/all_trends.csv"))

bidon=subset(trendsf,baseline>=1921)

nb_weird_1921=nrow(subset(bidon,abs(trend)>=1 & baseline==1921))
nb_weird_1971=nrow(subset(bidon,abs(trend)>=1 & baseline==1971))
nb_weird_1921
nb_weird_1971

bidon=subset(bidon,abs(trend)<1)
nrow(subset(bidon,baseline==1921))
nrow(subset(bidon,baseline==1971))
bidon %>% group_by(taxo_group,baseline) %>% summarise(length(unique(species)))

baselines_vec=unique(trendsf$baseline)
bidon$genus=sapply(strsplit(bidon$species, " "),function(x){x[[1]]})

length(unique(bidon$species[bidon$baseline==1921]))

ggplot(data=subset(bidon,taxo_group=="bees" & baseline==1921),aes(x=genus,y=det_prob))+geom_boxplot()+
  coord_flip()+ylab("Detection probability")


bidon$signi_des="non"
bidon$signi_des[bidon$described_pval<0.05]="yes"
sum(!is.na(bidon$described_eff) & bidon$baseline==1921)
sum(bidon$described_pval<0.05 & bidon$baseline==1921,na.rm=TRUE)
ggplot(data=subset(bidon,baseline==1921),
       aes(x=taxo_group,y=described_eff,color=signi_des))+geom_hline(yintercept=0,linetype="dashed")+geom_boxplot()+
  ylab("Effect of the publication of the species description on detection probability")+
  xlab("Group")+labs(color="Significant:")

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
