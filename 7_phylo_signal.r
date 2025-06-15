pkgs <- c("data.table", "dplyr","lme4","ggplot2","ggridges","metafor","cowplot","emmeans","tidyverse","ape","brms") 
inst <- pkgs %in% installed.packages()
if (any(inst)) install.packages(pkgs[!inst])
pkg_out <- lapply(pkgs, require, character.only = TRUE)
project_folder="C:/Users/Duchenne/Documents/safeguard/"

trendsf=fread(paste0(project_folder,"data/all_trends.csv"))
bidon=subset(trendsf,baseline==1921)
bidon=subset(bidon,abs(trend)<1)
bidon$phylo=gsub(" ","_",bidon$species)
bidon$genus=sapply(str_split(bidon$species," "),"[",1)


hp=read.tree(paste0(project_folder,"data/phylo_hovs_safeguard.tree"))
bp=read.tree(paste0(project_folder,"data/phylo_bees_safeguard.tree"))
bp=makeNodeLabel(bp)

hp=drop.tip(hp,hp$tip.label[!(hp$tip.label %in% bidon$phylo)])
bp=drop.tip(bp,bp$tip.label[!(bp$tip.label %in% bidon$phylo)])

Ah <- ape::vcv.phylo(hp)
Ab <- ape::vcv.phylo(bp)


bidon2=subset(bidon,taxo_group=="bees" & phylo %in% colnames(Ab))

model_bees <- brm(
  trend | mi(sde) ~ 1 +(1|region_50) +(1|FAMILY)+(1|genus)+(1|gr(phylo, cov = A)),
  data = bidon2,
  family = gaussian(),
  data2 = list(A = Ab),
  chains = 3,cores=3,iter =3000)
  
model_bees=model_simple

bidon2=subset(bidon,taxo_group=="hoverflies" & phylo %in% colnames(Ah))

model_hovs <- brm(
  trend | mi(sde) ~ 1 +(1|region_50)+(1|genus)+(1|gr(phylo, cov = A)),
  data = bidon2,
  family = gaussian(),
  data2 = list(A = Ah),
  chains = 3,cores=3,iter =3000)
  

models=list(model_bees,model_hovs)

save(models,file=paste0(project_folder,"data/variance_partition_brms.RData"))
###############################################################################
pkgs <- c("data.table", "dplyr","lme4","ggplot2","ggridges","metafor","cowplot","emmeans","tidyverse","ape","brms") 
inst <- pkgs %in% installed.packages()
if (any(inst)) install.packages(pkgs[!inst])
pkg_out <- lapply(pkgs, require, character.only = TRUE)
project_folder="C:/Users/Duchenne/Documents/safeguard/"

load(paste0(project_folder,"data/variance_partition_brms.RData"))
model_bees=models[[1]]
model_hovs=models[[2]]

obj=as_draws_df(model_bees)
names(obj)[1:6]=c("intercept","family","genus","phylo","region","sigma")
variance_tab=obj[,2:6]^2
variance_tab$tot=apply(variance_tab,1,sum)
variance_tab[,1:5]=variance_tab[,1:5]/variance_tab$tot
variance_tab$tot_norm=apply(variance_tab[,1:5],1,sum)

b=melt(variance_tab,id.vars=c("tot","tot_norm")) %>% group_by(variable) %>% summarise(moy=mean(value),lcl=quantile(value,probs=0.025),ucl=quantile(value,probs=0.975))
b$variable2=as.character(b$variable)
b$variable2[b$variable=="sigma"]="Residual"
b$variable2[b$variable=="phylo"]="Evolutionary history"
b$variable2[b$variable=="family"]="Taxonomy (fam.)"
b$variable2[b$variable=="genus"]="Taxonomy (genus)"
b$variable2[b$variable=="region"]="Biogeographic regions"
b=b[order(b$moy,decreasing =TRUE),]
b$variable2=factor(b$variable2,levels=unique(b$variable2))

p1=ggplot(data=b,aes(x=variable2,y=moy))+geom_bar(stat="identity")+
geom_errorbar(aes(ymin=lcl,ymax=ucl),width=0.1)+
scale_y_continuous(labels=scales::percent,limits=c(0,1))+
theme_bw()+
theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.border=element_blank(),
panel.grid.minor = element_blank(),panel.background = element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),
strip.background=element_rect(fill=NA,color=NA),legend.position="none",axis.text.x=element_text(angle=45,hjust=1))+ylab("Percentage of the total variance")+xlab("")

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
strip.background=element_rect(fill=NA,color=NA),legend.position="none",axis.text.x=element_text(angle=45,hjust=1))+ylab("Percentage of the total variance")+xlab("")


plot_grid(p1,p2,align="hv")