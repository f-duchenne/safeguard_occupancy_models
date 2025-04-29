pkgs <- c("data.table", "dplyr","lme4","ggplot2","ggridges","metafor","cowplot","emmeans","tidyverse") 
inst <- pkgs %in% installed.packages()
if (any(inst)) install.packages(pkgs[!inst])
pkg_out <- lapply(pkgs, require, character.only = TRUE)
project_folder="C:/Users/Duchenne/Documents/safeguard/"

colo2=c("#44AA99","#117733","#332288","#CC6677","#DDCC77")

inv.logit=function(x){exp(x)/(1+exp(x))}
logit=function(x){log(x/(1-x))}

#load trends
trendsf=fread("all_trends.csv")

####### FOCUS ON THE OVERALL PERIOD; BASELINE = 1921
bidon=subset(trendsf,baseline==1921)

#remove some weird trends
nb_weird=nrow(subset(bidon,abs(trend)>=1))
bidon=subset(bidon,abs(trend)<1)

#one table for each group
bidon_bee=subset(bidon,taxo_group=="bees")
bidon_hov=subset(bidon,taxo_group=="hoverflies")

#extract average estimate per region accoridng to meta-analyses model
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


plot_grid(pl1,pl2,ncol=1)
