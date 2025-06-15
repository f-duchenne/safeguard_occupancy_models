pkgs <- c("data.table", "dplyr","lme4","ggplot2","ggridges","metafor","cowplot","emmeans","tidyverse","gridExtra") 
inst <- pkgs %in% installed.packages()
if (any(inst)) install.packages(pkgs[!inst])
pkg_out <- lapply(pkgs, require, character.only = TRUE)
project_folder="C:/Users/Duchenne/Documents/safeguard/"

colo2=c("#44AA99","#117733","#332288","#CC6677","#DDCC77")

inv.logit=function(x){exp(x)/(1+exp(x))}
logit=function(x){log(x/(1-x))}

#load trends
trendsf=fread(paste0(project_folder,"data/all_trends.csv"))

####### FOCUS ON THE OVERALL PERIOD; BASELINE = 1921
bidon=subset(trendsf,baseline==1921)

#remove some weird trends
nb_weird=nrow(subset(bidon,abs(trend)>=1))
bidon=subset(bidon,abs(trend)<1)

#one table for each group
bidon_bee=subset(bidon,taxo_group=="bees")
bidon_hov=subset(bidon,taxo_group=="hoverflies")

#extract average estimate per region accoridng to meta-analyses model
bf=NULL
bf2=NULL
for(jj in unique(bidon$taxo_group)){
    load(paste0(project_folder,"data/model_",1921,"_",jj,".RData"))
	model_1=lis_bas[[1]]
	model_2=lis_bas[[2]]
	sav <- emmprep(model_1)
	b=as.data.frame(emmeans(sav,specs="region_50"))
	b$baseline=1921
	b$basic_mean=subset(bidon,taxo_group==jj) %>% group_by(region_50) %>% summarise(mean(trend)) %>% deframe()
	b$basic_SE=subset(bidon,taxo_group==jj) %>% group_by(region_50) %>% summarise(sd(trend)/sqrt(length(trend))) %>% deframe()
	b$taxo_group=jj
	b$n_trends=model_1$k
	b$subset_species="all"
	bf=rbind(bf,b)
	
	sav <- emmprep(model_2)
	b=as.data.frame(emmeans(sav,specs="region_50"))
	b$baseline=1921
	b$basic_mean=subset(bidon,taxo_group==jj & max.occ>=1e-4) %>% group_by(region_50) %>% summarise(mean(trend)) %>% deframe()
	b$basic_SE=subset(bidon,taxo_group==jj & max.occ>=1e-4) %>% group_by(region_50) %>% summarise(sd(trend)/sqrt(length(trend))) %>% deframe()
	b$taxo_group=jj
	b$n_trends=model_2$k
	b$subset_species="excluding rare"
	bf=rbind(bf,b)
	
	load(paste0(project_folder,"data/model_total_",1921,"_",jj,".RData"))
	model_1=lis_bas[[1]]
	model_2=lis_bas[[2]]
	sav <- emmprep(model_1)
	b=as.data.frame(emmeans(sav,specs="1"))
	b$baseline=1921
	b$taxo_group=jj
	b$n_trends=model_1$k
	b$subset_species="all"
	bf2=rbind(bf2,b)
	load(paste0(project_folder,"data/model_total_",1971,"_",jj,".RData"))
	model_1=lis_bas[[1]]
	model_2=lis_bas[[2]]
	sav <- emmprep(model_1)
	b=as.data.frame(emmeans(sav,specs="1"))
	b$baseline=1971
	b$taxo_group=jj
	b$n_trends=model_1$k
	b$subset_species="all"
	bf2=rbind(bf2,b)
}

dodgi=0.1

bidon_bee$region2=bidon_bee$region_50
bidon_bee$region2[bidon_bee$region2=="mediterranean"]="medit."
bidon_bee$region2[bidon_bee$region2=="continental"]="conti."
bidon_bee$rare="no"
bidon_bee$rare[bidon_bee$max.occ<=1e-4]="yes"
bidon_hov$region2=bidon_hov$region_50
bidon_hov$region2[bidon_hov$region2=="mediterranean"]="medit."
bidon_hov$region2[bidon_hov$region2=="continental"]="conti."
bidon_hov$rare="no"
bidon_hov$rare[bidon_hov$max.occ<=1e-4]="yes"

bf=as.data.frame(bf)
bf$region2=as.character(bf$region_50)
bf$region2[bf$region2=="mediterranean"]="medit."
bf$region2[bf$region2=="continental"]="conti."

table_S2=subset(bf, subset_species=="all")
n_tab=bidon %>% group_by(taxo_group,region_50) %>% count()
table_S2=merge(table_S2,n_tab,by=c("taxo_group","region_50"))
fwrite(table_S2,paste0(project_folder,"data/table_S2.csv"))

nrow(subset(bidon_bee,abs(trend)>0.2))/nrow(bidon_bee)
nrow(subset(bidon_hov,abs(trend)>0.2))/nrow(bidon_hov)

### PANEL A-B
pl1=ggplot()+
geom_density_ridges(data=bidon_bee,aes(x=trend,y=region2,fill=region2,color=region2),alpha=0.3,scale = 0.9)+
xlim(c(-0.2,0.2))+
geom_point(data=subset(bf,taxo_group=="bees" & subset_species=="all"),aes(x=emmean,y=as.numeric(as.factor(region2))+dodgi,fill=region2,color=region2),size=1.5)+
geom_errorbarh(data=subset(bf,taxo_group=="bees" & subset_species=="all"),aes(x=emmean,xmax=emmean+1.96*SE,xmin=emmean-1.96*SE,y=as.numeric(as.factor(region2))+dodgi,fill=region2,color=region2),height = 0,size=1)+
theme_bw()+geom_vline(xintercept=0,linetype="dashed")+
theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.border=element_blank(),
panel.grid.minor = element_blank(),panel.background = element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),
strip.background=element_rect(fill=NA,color=NA),legend.position="none")+
scale_color_manual(values=colo2)+scale_fill_manual(values=colo2)+ylab("Biogeographic region")+xlab(expression(paste("Species trend (log of growth rate)")))+coord_cartesian(expand=FALSE)+ggtitle("a")

pl2=ggplot()+
geom_density_ridges(data=bidon_hov,aes(x=trend,y=region2,fill=region2,color=region2),alpha=0.3,scale = 0.9)+
xlim(c(-0.2,0.2))+
geom_point(data=subset(bf,taxo_group=="hoverflies" & subset_species=="all"),aes(x=emmean,y=as.numeric(as.factor(region2))+dodgi,fill=region2,color=region2),size=1.5)+
geom_errorbarh(data=subset(bf,taxo_group=="hoverflies" & subset_species=="all"),aes(x=emmean,xmax=emmean+1.96*SE,xmin=emmean-1.96*SE,y=as.numeric(as.factor(region2))+dodgi,fill=region2,color=region2),height = 0,size=1)+
theme_bw()+geom_vline(xintercept=0,linetype="dashed")+
theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.border=element_blank(),
panel.grid.minor = element_blank(),panel.background = element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),
strip.background=element_rect(fill=NA,color=NA),legend.position="none")+
scale_color_manual(values=colo2)+scale_fill_manual(values=colo2)+ylab("Biogeographic region")+xlab(expression(paste("Species trend (log of growth rate)")))+coord_cartesian(expand=FALSE)+ggtitle("b")

### PANEL C-D
bidon_bee$winlo=NA
bidon_bee$winlo[bidon_bee$trend>0 & bidon_bee$significant=="yes"]="Winners"
bidon_bee$winlo[bidon_bee$trend<0 & bidon_bee$significant=="yes"]="Losers"

bidon_hov$winlo=NA
bidon_hov$winlo[bidon_hov$trend>0 & bidon_hov$significant=="yes"]="Winners"
bidon_hov$winlo[bidon_hov$trend<0 & bidon_hov$significant=="yes"]="Losers"

pl3=ggplot(data=subset(bidon_bee,!is.na(winlo)),aes(x=winlo,fill=region2),position="stack")+
geom_bar()+
theme_bw()+
theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.border=element_blank(),
panel.grid.minor = element_blank(),panel.background = element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),
strip.background=element_rect(fill=NA,color=NA),legend.position="none")+
scale_color_manual(values=colo2)+scale_fill_manual(values=colo2)+ylab("nb. of species")+xlab("")+coord_cartesian(expand=FALSE)+ggtitle("c")

pl4=ggplot(data=subset(bidon_hov,!is.na(winlo)),aes(x=winlo,fill=region2),position="stack")+
geom_bar()+
theme_bw()+
theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.border=element_blank(),
panel.grid.minor = element_blank(),panel.background = element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),
strip.background=element_rect(fill=NA,color=NA),legend.position="none")+
scale_color_manual(values=colo2)+scale_fill_manual(values=colo2)+ylab("nb. of species")+xlab("")+coord_cartesian(expand=FALSE)+ggtitle("d")

p1=plot_grid(pl1,pl3,pl2,pl4,ncol=2,rel_widths=c(2/4,1/4))


### PANEL E
pl5=ggplot(data=subset(bf2,subset_species=="all"),aes(x=taxo_group,y=exp(emmean)-1,color=as.factor(baseline)))+
geom_errorbar(aes(ymax=exp(emmean+1.96*SE)-1,ymin=exp(emmean-1.96*SE)-1),width = 0,size=1,position = position_dodge(width = 0.4))+
geom_point(size=2,position = position_dodge(width = 0.4))+
theme_bw()+geom_hline(yintercept=0,linetype="dashed")+
theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.border=element_blank(),
panel.grid.minor = element_blank(),panel.background = element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),
strip.background=element_rect(fill=NA,color=NA),legend.position="right")+
xlab("Group")+ylab("European average occupancy trend (%/year)")+ggtitle("e")+
scale_y_continuous(labels=scales::percent)+labs(color="baseline")+
scale_color_manual(values=c("black","grey"))

grid.arrange(p1,pl5,ncol=2,widths=c(2,0.75))


pdf(paste0(project_folder,"/Figure_2.pdf"),width=8,height=5)
grid.arrange(p1,pl5,ncol=2,widths=c(2,1))
dev.off();


########################################
bidon_bee$group_abund=cut(log(bidon_bee$nb_records),10)
ggplot(data=bidon_bee,aes(y=abs(trend),x=group_abund))+
geom_boxplot()+
theme_bw()+geom_vline(xintercept=0,linetype="dashed")+
theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.border=element_blank(),
panel.grid.minor = element_blank(),panel.background = element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),
strip.background=element_rect(fill=NA,color=NA),legend.position="right")+
scale_color_manual(values=colo2)+xlab("Biogeographic region")+ylab(expression(paste("Species trend (log of growth rate)")))+ggtitle("a")


ggplot()+
geom_point(data=subset(bf,taxo_group=="hoverflies" & subset_species=="all"),aes(x=emmean,y=as.numeric(as.factor(region2)),fill=region2,color=region2),size=2)+
geom_errorbar(data=subset(bf,taxo_group=="hoverflies" & subset_species=="all"),aes(x=emmean,xmax=emmean+1.96*SE,xmin=emmean-1.96*SE,y=as.numeric(as.factor(region2)),fill=region2,color=region2),height = 0,size=1)


