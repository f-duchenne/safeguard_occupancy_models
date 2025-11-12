pkgs <- c("data.table", "dplyr","lme4","ggplot2","ggridges","metafor","cowplot","emmeans","tidyverse","ape","brms") 
inst <- pkgs %in% installed.packages()
if (any(inst)) install.packages(pkgs[!inst])
pkg_out <- lapply(pkgs, require, character.only = TRUE)
project_folder="C:/Users/Duchenne/Documents/safeguard/"

#MERGE TRENDS AND TRAITS
trends=fread(paste0(project_folder,"data/all_trends.csv"))
traits=fread(paste0(project_folder,"data/traits_table.csv"))
traits[traits==""]=NA
trends=merge(trends,traits,by.x=c("species","taxo_group"),by.y=c("Species","taxo_group"),all.x=TRUE,all.y=FALSE)

#### RUNNING MODELS (MODELS TAKE FEW HOURS TO RUN, YOU CAN SKIP THAT SCRIPT AND LOAD THEM DIRECTLY TO DO FIGURE 4)

nit=20000 #set number of iterations for bayesian models

#### BEES
trends_bee=subset(trends,taxo_group=="bees" & baseline==1921)
trends_bee=subset(trends_bee,!is.na(Larval_diet_breadth) & !is.na(ITD_F_mm) & !is.na(STI_Species_temperature_index) & !is.na(Sociality_simplified))
trends_bee$phylo=gsub(" ","_",trends_bee$species)
trends_bee$genus=sapply(str_split(trends_bee$species," "),"[",1)

bp=read.tree(paste0(project_folder,"data/phylo_bees_safeguard.tree"))
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

hp=read.tree(paste0(project_folder,"data/phylo_hovs_safeguard.tree"))

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

save(models,file=paste0(project_folder,"final_and_intermediate_outputs/model_traits/model traits/traits_models_brms.RData"))

###############################################################################
pkgs <- c("data.table", "dplyr","lme4","ggplot2","ggridges","metafor","cowplot","emmeans","tidyverse","ape","brms") 
inst <- pkgs %in% installed.packages()
if (any(inst)) install.packages(pkgs[!inst])
pkg_out <- lapply(pkgs, require, character.only = TRUE)
project_folder="C:/Users/Duchenne/Documents/safeguard/"

load(paste0(project_folder,"data/model traits/traits_models_brms_taxo.RData"))
model_bees=models[[1]]
model_hovs=models[[2]]
model_bees_null=models[[3]]
model_hovs_null=models[[4]]

obj=as_draws_df(model_bees)
names(obj)[c(1,7:9)]=c("intercept","phylo","region","sigma")
obj$traits=as_draws_df(model_bees_null)$sigma-obj$sigma
variance_tab=obj[,c("phylo","region","sigma","traits")]^2
variance_tab$tot=apply(variance_tab,1,sum)
variance_tab=variance_tab/variance_tab$tot
variance_tab$tot_norm=apply(variance_tab,1,sum)

b=melt(variance_tab,id.vars=c("tot","tot_norm")) %>% group_by(variable) %>% summarise(moy=mean(value),lcl=quantile(value,probs=0.025),ucl=quantile(value,probs=0.975))
b$variable2=as.character(b$variable)
b$variable2[b$variable=="sigma"]="Residual"
b$variable2[b$variable=="phylo"]="Evolutionary history"
b$variable2[b$variable=="region"]="Biogeographic regions"
b$variable2[b$variable=="traits"]="Species traits"
b=b[order(b$moy,decreasing =TRUE),]
b$variable2=factor(b$variable2,levels=unique(b$variable2))

p1=ggplot(data=b,aes(x=variable2,y=moy))+geom_bar(stat="identity")+
geom_errorbar(aes(ymin=lcl,ymax=ucl),width=0.1)+
scale_y_continuous(labels=scales::percent,limits=c(0,1))+
theme_bw()+
theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.border=element_blank(),
panel.grid.minor = element_blank(),panel.background = element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),
strip.background=element_rect(fill=NA,color=NA),legend.position="none",axis.text.x=element_text(angle=45,hjust=1))+ylab("Percentage of the total variance")+xlab("")+ggtitle("a")

obj=as_draws_df(model_hovs)
names(obj)[1:5]=c("intercept","genus","phylo","region","sigma")
variance_tab=obj[,2:5]^2
variance_tab$tot=apply(variance_tab,1,sum)
variance_tab[,1:4]=variance_tab[,1:4]/variance_tab$tot
variance_tab$tot_norm=apply(variance_tab[,1:4],1,sum)

b=melt(variance_tab,id.vars=c("tot","tot_norm")) %>% group_by(variable) %>% summarise(moy=mean(value),lcl=quantile(value,probs=0.025),ucl=quantile(value,probs=0.975))
b$variable2=as.character(b$variable)
b$variable2[b$variable=="sigma"]="Residual"
b$variable2[b$variable=="phylo"]="Evolutionary history"
b$variable2[b$variable=="genus"]="Taxonomy (genus)"
b$variable2[b$variable=="region"]="Biogeographic regions"
b=b[order(b$moy,decreasing =TRUE),]
b$variable2=factor(b$variable2,levels=unique(b$variable2))

p2=ggplot(data=b,aes(x=variable2,y=moy))+geom_bar(stat="identity")+
geom_errorbar(aes(ymin=lcl,ymax=ucl),width=0.1)+
scale_y_continuous(labels=scales::percent)+
theme_bw()+
theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.border=element_blank(),
panel.grid.minor = element_blank(),panel.background = element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),
strip.background=element_rect(fill=NA,color=NA),legend.position="none",axis.text.x=element_text(angle=45,hjust=1))+ylab("Percentage of the total variance")+xlab("")+ggtitle("b")


plot_grid(p1,p2,align="hv")

pdf(paste0(project_folder,"Figure_3.pdf"),width=7,height=5)
plot_grid(p1,p2,align="hv")
dev.off();
