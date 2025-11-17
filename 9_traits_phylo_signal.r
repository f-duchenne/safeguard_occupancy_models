pkgs <- c("data.table", "dplyr","lme4","ggplot2","ggridges","metafor","cowplot","emmeans","tidyverse","ape","brms") 
inst <- pkgs %in% installed.packages()
if (any(inst)) install.packages(pkgs[!inst])
pkg_out <- lapply(pkgs, require, character.only = TRUE)
#project_folder="C:/Users/Duchenne/Documents/safeguard/"
project_folder=""

#MERGE TRENDS AND TRAITS
trends=fread(paste0(project_folder,"data/final_and_intermediate_outputs/all_trends.csv"))
traits=fread(paste0(project_folder,"data/final_and_intermediate_outputs/traits_table.csv"))
traits[traits==""]=NA
trends=merge(trends,traits,by.x=c("species","taxo_group"),by.y=c("Species","taxo_group"),all.x=TRUE,all.y=FALSE)
trends=subset(trends,abs(trend)<1)

#### RUNNING MODELS (MODELS TAKE FEW HOURS TO RUN, YOU CAN SKIP THAT SCRIPT AND LOAD THEM DIRECTLY TO DO FIGURE 4)

nit=20000 #set number of iterations for bayesian models

#### BEES
trends_bee=subset(trends,taxo_group=="bees" & baseline==1921)
trends_bee=subset(trends_bee,!is.na(Larval_diet_breadth) & !is.na(ITD_F_mm) & !is.na(STI_Species_temperature_index) & !is.na(Sociality_simplified))
trends_bee$phylo=gsub(" ","_",trends_bee$species)
trends_bee$genus=sapply(str_split(trends_bee$species," "),"[",1)

bp=read.tree(paste0(project_folder,"data/raw_data/phylo_bees_safeguard.tree"))
bp=makeNodeLabel(bp)

bp=drop.tip(bp,bp$tip.label[!(bp$tip.label %in% trends_bee$phylo)])

Ab <- ape::vcv.phylo(bp)

subset(trends_bee,taxo_group=="bees" & !(phylo %in% colnames(Ab)))
bidon=subset(trends_bee,taxo_group=="bees" & phylo %in% colnames(Ab))

model_bees <- brm(
  trend | mi(sde) ~ 1 + ITD_F_mm + Sociality_simplified + STI_Species_temperature_index +
  Larval_diet_breadth + (1|region_50)+(1|gr(phylo, cov = A)),
  data = bidon,
  family = gaussian(),
  data2 = list(A = Ab),
  chains = 3,cores=3,iter =nit)
 
 model_bees_null <- brm(
  trend | mi(sde) ~ 1+(1|region_50)+(1|gr(phylo, cov = A)),
  data = bidon,
  family = gaussian(),
  data2 = list(A = Ab),
  chains = 3,cores=3,iter =nit)

#### HOVERFLIES
trends_hov=subset(trends,taxo_group=="hoverflies" & baseline==1921)
trends_hov=subset(trends_hov,!is.na(Larval_nutrition) & !is.na(Adult_Body_size_num) & !is.na(STI_Species_temperature_index) & !is.na(Flight_height))
trends_hov$phylo=gsub(" ","_",trends_hov$species)
trends_hov$genus=sapply(str_split(trends_hov$species," "),"[",1)

hp=read.tree(paste0(project_folder,"data/raw_data/phylo_hovs_safeguard.tree"))

hp=drop.tip(hp,hp$tip.label[!(hp$tip.label %in% trends_hov$phylo)])

Ah <- ape::vcv.phylo(hp)

bidon2=subset(trends_hov,taxo_group=="hoverflies" & phylo %in% colnames(Ah))

model_hovs <- brm(
  trend | mi(sde) ~ 1+Adult_Body_size_num +STI_Species_temperature_index+ Flight_height+Larval_nutrition+(1|region_50)+(1|gr(phylo, cov = A)),
  data = bidon2,
  family = gaussian(),
  data2 = list(A = Ah),
  chains = 3,cores=3,iter =nit)
  
model_hovs_null <- brm(
  trend | mi(sde) ~ 1+(1|region_50)+(1|gr(phylo, cov = A)),
  data = bidon2,
  family = gaussian(),
  data2 = list(A = Ah),
  chains = 3,cores=3,iter =nit)
  
models=list(model_bees,model_hovs,model_bees_null,model_hovs_null)

save(models,file=paste0(project_folder,"data/final_and_intermediate_outputs/traits_models_brms.RData"))
#
