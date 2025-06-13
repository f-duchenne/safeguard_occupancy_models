
############# FIRST ASSEMBLE ALL SPECIES TOGETHER
pkgs <- c("data.table", "dplyr") 
inst <- pkgs %in% installed.packages()
if (any(inst)) install.packages(pkgs[!inst])
pkg_out <- lapply(pkgs, require, character.only = TRUE)

project_folder="C:/Users/Duchenne/Documents/safeguard/"

setwd(dir=paste0(project_folder,"results_glm_with_counts"))
lifile=list.files( pattern ="predicts_")

trendsf=NULL
for(i in lifile){
index=gsub(".RData","",gsub("predicts_","",i))

load(i)
trends=lili[[1]]
trendsf=rbind(trendsf,trends)
}

fwrite(trendsf,paste0(project_folder,"data/results_TMB_trends_with_counts.csv"))

############# SECOND PLOT RESULTS:
pkgs <- c("data.table", "dplyr","lme4","ggplot2","ggridges","metafor","cowplot","emmeans","tidyverse") 
inst <- pkgs %in% installed.packages()
if (any(inst)) install.packages(pkgs[!inst])
pkg_out <- lapply(pkgs, require, character.only = TRUE)
project_folder="C:/Users/Duchenne/Documents/safeguard/"

colo2=c("#44AA99","#117733","#332288","#CC6677","#DDCC77")

inv.logit=function(x){exp(x)/(1+exp(x))}
logit=function(x){log(x/(1-x))}

tab_spec=fread(paste0(project_folder,"data/species_nb_records.csv"))
tab_spec=subset(tab_spec,nb_detect>=5 & nb_records_tot>=10) %>% group_by(TAXON) %>% mutate(nb_detect_tot=sum(nb_detect),occ_min=min(occupancy_obs))

trendsf=fread(paste0(project_folder,"data/results_TMB_trends_with_counts.csv"))
trendsf=subset(trendsf,species!="Apis mellifera")
nrow(trendsf)
trendsf=merge(trendsf,tab_spec,by.x=c("species","region_50","taxo_group"),by.y=c("TAXON","region_50","taxo_group"))
nrow(trendsf)

b=unique(trendsf[,c("species","nb_records_tot","nb_detect_tot","convergence","occ_min")])

trendsf$significant="no"
trendsf$significant[trendsf$trend<0 & (trendsf$trend+1.96*trendsf$sde)<0]="yes"
trendsf$significant[trendsf$trend>0 & (trendsf$trend-1.96*trendsf$sde)>0]="yes"

trendsf$growth.rate=100*(exp(trendsf$trend)-1)
trendsf=subset(trendsf,convergence==0 & !is.na(acim) & !is.na(sde))

fwrite(trendsf,"all_trends.csv")

####### FOCUS ON THE OVERALL PERIOD; BASELINE = 1921
bidon=subset(trendsf,baseline==1921)
nb_weird=nrow(subset(bidon,abs(trend)>=1))
bidon=subset(bidon,abs(trend)<1)

tab_spec$trend_available="no"
tab_spec$trend_available[tab_spec$TAXON %in% bidon$species]="yes"

nb_not_converged=sum(!(unique(tab_spec$TAXON[tab_spec$taxo_group=="bees"]) %in% bidon$species))
nb_converged=length(unique(bidon$species))
nb_combi=nrow(bidon)

bidon_bee=subset(bidon,taxo_group=="bees")
bidon_hov=subset(bidon,taxo_group=="hoverflies")
model_bee=rma(trend ~ region_50+FAMILY, sei=sde,data=bidon_bee, digits=10)
model_hov=rma(trend ~ region_50+FAMILY, sei=sde,data=bidon_hov, digits=10)
# model_bee_2=rma(trend ~ region_50, sei=sde,data=subset(bidon_bee,max.occ>1e-4), digits=10)
# model_hov_2=rma(trend ~ region_50, sei=sde,data=subset(bidon_hov,max.occ>1e-4), digits=10)
# lis=list(model_bee,model_hov,model_bee_2,model_hov_2)
# save(lis,file=paste0(project_folder,"data/models_total.RData"))
sav <- emmprep(model_bee)
b=as.data.frame(emmeans(sav,specs="FAMILY"))

ggplot()+
geom_pointrange(data=b,aes(y=emmean,x= FAMILY ,fill= FAMILY ,color= FAMILY ,ymax=emmean+1.96*SE,ymin=emmean-1.96*SE),size=2)+
theme_bw()+geom_vline(xintercept=0,linetype="dashed")+
theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.border=element_blank(),
panel.grid.minor = element_blank(),panel.background = element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),
strip.background=element_rect(fill=NA,color=NA),legend.position="none")+
ylab("trend")


load(paste0(project_folder,"data/models_total.RData"))
bf=NULL
for(jj in unique(bidon$taxo_group)){
	if(jj=="bees"){
		model_1=lis[[1]]
		model_2=lis[[3]]
	}else{
		model_1=lis[[2]]
		model_2=lis[[4]]
	}
	sav <- emmprep(model_1)
	b=as.data.frame(emmeans(sav,specs="region_50"))
	b$baseline=1921
	b$basic_mean=subset(bidon,taxo_group==jj) %>% group_by(region_50) %>% summarise(mean(trend)) %>% deframe()
	b$basic_SE=subset(bidon,taxo_group==jj) %>% group_by(region_50) %>% summarise(sd(trend)/sqrt(length(trend))) %>% deframe()
	b$taxo_group=jj
	b$subset_species="all"
	bf=rbind(bf,b)
	
	sav <- emmprep(model_2)
	b=as.data.frame(emmeans(sav,specs="region_50"))
	b$baseline=1921
	b$basic_mean=subset(bidon,taxo_group==jj & max.occ>=1e-4) %>% group_by(region_50) %>% summarise(mean(trend)) %>% deframe()
	b$basic_SE=subset(bidon,taxo_group==jj & max.occ>=1e-4) %>% group_by(region_50) %>% summarise(sd(trend)/sqrt(length(trend))) %>% deframe()
	b$taxo_group=jj
	b$subset_species="excluding rare"
	bf=rbind(bf,b)
}

dodgi=0.1

pl1=ggplot()+
geom_density_ridges(data=bidon_bee,aes(x=trend,y=region_50,fill=region_50,color=region_50),alpha=0.3,scale = 0.9)+
xlim(c(-0.2,0.2))+
geom_point(data=subset(bf,taxo_group=="bees" & subset_species=="all"),aes(x=emmean,y=as.numeric(as.factor(region_50))+dodgi,fill=region_50,color=region_50),size=2)+
geom_errorbarh(data=subset(bf,taxo_group=="bees" & subset_species=="all"),aes(x=emmean,xmax=emmean+1.96*SE,xmin=emmean-1.96*SE,y=as.numeric(as.factor(region_50))+dodgi,fill=region_50,color=region_50),height = 0,size=1)+
theme_bw()+geom_vline(xintercept=0,linetype="dashed")+
theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.border=element_blank(),
panel.grid.minor = element_blank(),panel.background = element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),
strip.background=element_rect(fill=NA,color=NA),legend.position="none")+
scale_color_manual(values=colo2)+scale_fill_manual(values=colo2)+ylab("Bio-geographic region")+xlab(expression(paste("Species trend (log of growth rate)")))+coord_cartesian(expand=FALSE)+ggtitle("a")

pl1b=ggplot()+
geom_density_ridges(data=subset(bidon_bee,max.occ>1e-4),aes(x=trend,y=region_50,fill=region_50,color=region_50),alpha=0.3,scale = 0.9)+
xlim(c(-0.2,0.2))+
geom_point(data=subset(bf,taxo_group=="bees" & subset_species=="excluding rare"),aes(x=emmean,y=as.numeric(as.factor(region_50))+dodgi,fill=region_50,color=region_50),size=2)+
geom_errorbarh(data=subset(bf,taxo_group=="bees" & subset_species=="excluding rare"),aes(x=emmean,xmax=emmean+1.96*SE,xmin=emmean-1.96*SE,y=as.numeric(as.factor(region_50))+dodgi,fill=region_50,color=region_50),height = 0,size=1)+
theme_bw()+geom_vline(xintercept=0,linetype="dashed")+
theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.border=element_blank(),
panel.grid.minor = element_blank(),panel.background = element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),
strip.background=element_rect(fill=NA,color=NA),legend.position="none")+
scale_color_manual(values=colo2)+scale_fill_manual(values=colo2)+ylab("Bio-geographic region")+xlab(expression(paste("Species trend (log of growth rate)")))+coord_cartesian(expand=FALSE)+ggtitle("a")

pl2=ggplot()+
geom_density_ridges(data=bidon_hov,aes(x=trend,y=region_50,fill=region_50,color=region_50),alpha=0.3,scale = 0.9)+
xlim(c(-0.2,0.2))+
geom_point(data=subset(bf,taxo_group=="hoverflies" & subset_species=="all"),aes(x=emmean,y=as.numeric(as.factor(region_50))+dodgi,fill=region_50,color=region_50),size=2)+
geom_errorbarh(data=subset(bf,taxo_group=="hoverflies" & subset_species=="all"),aes(x=emmean,xmax=emmean+1.96*SE,xmin=emmean-1.96*SE,y=as.numeric(as.factor(region_50))+dodgi,fill=region_50,color=region_50),height = 0,size=1)+
theme_bw()+geom_vline(xintercept=0,linetype="dashed")+
theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.border=element_blank(),
panel.grid.minor = element_blank(),panel.background = element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),
strip.background=element_rect(fill=NA,color=NA),legend.position="none")+
scale_color_manual(values=colo2)+scale_fill_manual(values=colo2)+ylab("Bio-geographic region")+xlab(expression(paste("Species trend (log of growth rate)")))+coord_cartesian(expand=FALSE)+ggtitle("b")

pl2b=ggplot()+
geom_density_ridges(data=subset(bidon_hov,max.occ>1e-4),aes(x=trend,y=region_50,fill=region_50,color=region_50),alpha=0.3,scale = 0.9)+
xlim(c(-0.2,0.2))+
geom_point(data=subset(bf,taxo_group=="hoverflies" & subset_species=="all"),aes(x=emmean,y=as.numeric(as.factor(region_50))+dodgi,fill=region_50,color=region_50),size=2)+
geom_errorbarh(data=subset(bf,taxo_group=="hoverflies" & subset_species=="all"),aes(x=emmean,xmax=emmean+1.96*SE,xmin=emmean-1.96*SE,y=as.numeric(as.factor(region_50))+dodgi,fill=region_50,color=region_50),height = 0,size=1)+
theme_bw()+geom_vline(xintercept=0,linetype="dashed")+
theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.border=element_blank(),
panel.grid.minor = element_blank(),panel.background = element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),
strip.background=element_rect(fill=NA,color=NA),legend.position="none")+
scale_color_manual(values=colo2)+scale_fill_manual(values=colo2)+ylab("Bio-geographic region")+xlab(expression(paste("Species trend (log of growth rate)")))+coord_cartesian(expand=FALSE)+ggtitle("b")


plot_grid(pl1,pl2,ncol=1)

####### BASELINE shift
bidon=subset(trendsf,baseline>1921)
nb_weird=nrow(subset(bidon,abs(trend)>=1))
bidon=subset(bidon,abs(trend)<1)

# for(jj in unique(bidon$taxo_group)){
# for(j in 2:length(baselines_vec)){
# bidon2=subset(bidon,convergence==0 & !is.na(acim) & baseline==baselines_vec[j] & !is.na(sde) & taxo_group==jj)

# model=rma(trend ~ region_50, sei=sde,data=bidon2, digits=10)
# model2=rma(trend ~ region_50, sei=sde,data=subset(bidon2,max.occ>1e-4), digits=10)
# lis_bas=list(model,model2)
# save(lis_bas,file=paste0(project_folder,"data/model_",baselines_vec[j],"_",jj,".RData"))
# }
# }
baselines_vec=sort(unique(bidon$baseline))
for(jj in unique(bidon$taxo_group)){
for(j in 1:length(baselines_vec)){
bidon2=subset(bidon,convergence==0 & !is.na(acim) & baseline==baselines_vec[j] & !is.na(sde) & taxo_group==jj)
load(paste0(project_folder,"data/model_",baselines_vec[j],"_",jj,".RData"))
model=lis_bas[[1]]
boxplot(weights(model))
sav <- emmprep(model)
bb=as.data.frame(emmeans(sav,specs="region_50"))
bb$baseline=baselines_vec[j]
bb$basic_mean=bidon2 %>% group_by(region_50) %>% summarise(mean(trend)) %>% deframe()
bb$basic_SE=bidon2 %>% group_by(region_50) %>% summarise(sd(trend)/sqrt(length(trend))) %>% deframe()
bb$taxo_group=jj
bb$subset_species="all"
bf=rbind(bf,bb)

model=lis_bas[[2]]
boxplot(weights(model))
sav <- emmprep(model)
bb=as.data.frame(emmeans(sav,specs="region_50"))
bb$baseline=baselines_vec[j]
bb$basic_mean=bidon2 %>% group_by(region_50) %>% summarise(mean(trend)) %>% deframe()
bb$basic_SE=bidon2 %>% group_by(region_50) %>% summarise(sd(trend)/sqrt(length(trend))) %>% deframe()
bb$taxo_group=jj
bb$subset_species="excluding rare"
bf=rbind(bf,bb)
}
}

pl3=ggplot()+
geom_point(data=subset(bf,taxo_group=="bees" & subset_species=="excluding rare"),aes(y=emmean,x=baseline,fill=region_50,color=region_50),position=position_dodge(width=2),size=2)+
geom_errorbar(data=subset(bf,taxo_group=="bees" & subset_species=="excluding rare"),aes(y=emmean,ymax=emmean+1.96*SE,ymin=emmean-1.96*SE,x=baseline,fill=region_50,color=region_50),position=position_dodge(width=2),width=0,size=1)+
theme_bw()+geom_hline(yintercept=0,linetype="dashed")+
theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.border=element_blank(),
panel.grid.minor = element_blank(),panel.background = element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),
strip.background=element_rect(fill=NA,color=NA),legend.position="none")+
scale_color_manual(values=colo2)+scale_fill_manual(values=colo2)+xlab("Baseline year")+ylab(expression(paste("Species trend (log of growth rate)")))+ggtitle("c")


pl3b=ggplot()+
geom_point(data=subset(bf,taxo_group=="bees" & subset_species=="excluding rare"),aes(y=emmean,x=baseline,fill=region_50,color=region_50),position=position_dodge(width=2),size=2)+
geom_errorbar(data=subset(bf,taxo_group=="bees" & subset_species=="excluding rare"),aes(y=emmean,ymax=emmean+1.96*SE,ymin=emmean-1.96*SE,x=baseline,fill=region_50,color=region_50),position=position_dodge(width=2),width=0,size=1)+
theme_bw()+geom_hline(yintercept=0,linetype="dashed")+
theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.border=element_blank(),
panel.grid.minor = element_blank(),panel.background = element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),
strip.background=element_rect(fill=NA,color=NA),legend.position="none")+
scale_color_manual(values=colo2)+scale_fill_manual(values=colo2)+xlab("Baseline year")+ylab(expression(paste("Species trend (log of growth rate)")))+ggtitle("c")

pl4=ggplot()+
geom_point(data=subset(bf,taxo_group=="hoverflies" & subset_species=="all"),aes(y=emmean,x=baseline,fill=region_50,color=region_50),position=position_dodge(width=2),size=2)+
geom_errorbar(data=subset(bf,taxo_group=="hoverflies" & subset_species=="all"),aes(y=emmean,ymax=emmean+1.96*SE,ymin=emmean-1.96*SE,x=baseline,fill=region_50,color=region_50),position=position_dodge(width=2),width=0,size=1)+
theme_bw()+geom_hline(yintercept=0,linetype="dashed")+
theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.border=element_blank(),
panel.grid.minor = element_blank(),panel.background = element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),
strip.background=element_rect(fill=NA,color=NA),legend.position="none")+
scale_color_manual(values=colo2)+scale_fill_manual(values=colo2)+xlab("Baseline year")+ylab(expression(paste("Species trend (log of growth rate)")))+ggtitle("d")

pl4b=ggplot()+
geom_point(data=subset(bf,taxo_group=="hoverflies" & subset_species=="excluding rare"),aes(y=emmean,x=baseline,fill=region_50,color=region_50),position=position_dodge(width=2),size=2)+
geom_errorbar(data=subset(bf,taxo_group=="hoverflies" & subset_species=="excluding rare"),aes(y=emmean,ymax=emmean+1.96*SE,ymin=emmean-1.96*SE,x=baseline,fill=region_50,color=region_50),position=position_dodge(width=2),width=0,size=1)+
theme_bw()+geom_hline(yintercept=0,linetype="dashed")+
theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.border=element_blank(),
panel.grid.minor = element_blank(),panel.background = element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),
strip.background=element_rect(fill=NA,color=NA),legend.position="none")+
scale_color_manual(values=colo2)+scale_fill_manual(values=colo2)+xlab("Baseline year")+ylab(expression(paste("Species trend (log of growth rate)")))+ggtitle("d")

plot_grid(pl1,pl3,pl2,pl4,ncol=2,align="hv")

pdf(paste0(project_folder,"Figure_2.pdf"),width=8,height=6)
plot_grid(pl1,pl3,pl2,pl4,ncol=2,rel_widths=c(1.2,1))
dev.off();

pdf(paste0(project_folder,"Figure_S2.pdf"),width=8,height=6)
plot_grid(pl1b,pl3b,pl2b,pl4b,ncol=2,rel_widths=c(1.2,1))
dev.off();
############################ prop loosers and winners
trendsf$winner=0
trendsf$winner[trendsf$trend>0 & trendsf$significant=="yes"]=1

trendsf$looser=0
trendsf$looser[trendsf$trend<0 & trendsf$significant=="yes"]=1

tabi=subset(trendsf,baseline==1921) %>% group_by(region_50,taxo_group) %>% summarise(prop_looser=round(sum(looser)/length(looser),2),prop_winner=round(sum(winner)/length(winner),2))


subset(trendsf,baseline==2001 & taxo_group=="hoverflies") %>% group_by(region_50,taxo_group) %>% summarise(prop_looser=round(sum(looser)/length(looser),2),prop_winner=round(sum(winner)/length(winner),2))


################################### Variance same species
bidon=subset(trendsf,convergence==0 & !is.na(acim) & baseline2==1921 & !is.na(sde))
nb_weird=nrow(subset(bidon,abs(trend)>=1))
bidon=subset(bidon,abs(trend)<1)

bidon$cate=1
bidon$cate[bidon$region_50 %in% c("atlantic")]=2
bidon$cate[bidon$region_50 %in% c("continental")]=3
bidon$cate[bidon$region_50 %in% c("boreal","alpine")]=4
bidon=bidon %>% group_by(species) %>% mutate(nr=length(unique(cate)),min.cate=min(cate),max.cate=max(cate))
bidon$cate2="middle"
bidon$cate2[bidon$cate==bidon$min.cate]="south edge"
bidon$cate2[bidon$cate==bidon$max.cate]="north edge"

bidon=subset(bidon,nr>1)

model=lmer(trend~cate2+(1|FAMILY)+(1|species),data=bidon)

be=ggpredict(model,"cate2",type="fixed")
plot(be)+geom_hline(yintercept=0)+xlab("")