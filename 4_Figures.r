
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
nrow(trendsf)
trendsf=merge(trendsf,tab_spec,by.x=c("species","region_50"),by.y=c("TAXON","region_50"))
nrow(trendsf)

b=unique(trendsf[,c("species","nb_records_tot","nb_detect_tot","convergence","occ_min")])

trendsf$significant="no"
trendsf$significant[trendsf$trend<0 & (trendsf$trend+1.96*trendsf$sde)<0]="yes"
trendsf$significant[trendsf$trend>0 & (trendsf$trend-1.96*trendsf$sde)>0]="yes"

trendsf$growth.rate=100*(exp(trendsf$trend)-1)

baselines_vec=c(1921,1951,1961,1971,1981,1991,2001)
trendsf$baseline2=sapply(trendsf$baseline,function(x){baselines_vec[x]})

####### FOCUS ON THE OVERALL PERIOD; BASELINE = 1921
bidon=subset(trendsf,convergence==0 & !is.na(acim) & baseline2==1921 & !is.na(sde))
nb_weird=nrow(subset(bidon,abs(trend)>=1))
bidon=subset(bidon,abs(trend)<1)

tab_spec$trend_available="no"
tab_spec$trend_available[tab_spec$TAXON %in% bidon$species]="yes"

nb_not_converged=sum(!(unique(tab_spec$TAXON[tab_spec$taxo_group=="bees"]) %in% bidon$species))
nb_converged=length(unique(bidon$species))
nb_combi=nrow(bidon)


# model=rma(trend ~ region_50, sei=sde,data=bidon, digits=10)
# save(model,file=paste0(project_folder,"data/model_total.RData"))
load(paste0(project_folder,"data/model_total.RData"))
sav <- emmprep(model)
b=as.data.frame(emmeans(sav,specs="region_50"))
b$baseline=1921
b$basic_mean=bidon %>% group_by(region_50) %>% summarise(mean(trend)) %>% deframe()
b$basic_SE=bidon %>% group_by(region_50) %>% summarise(sd(trend)/sqrt(length(trend))) %>% deframe()

dodgi=0.1

pl1=ggplot()+
geom_density_ridges(data=bidon,aes(x=trend,y=region_50,fill=region_50,color=region_50),alpha=0.3,scale = 0.9)+
xlim(c(-0.2,0.2))+
geom_point(data=b,aes(x=emmean,y=as.numeric(as.factor(region_50))+dodgi,fill=region_50,color=region_50),size=2)+
geom_errorbarh(data=b,aes(x=emmean,xmax=emmean+1.96*SE,xmin=emmean-1.96*SE,y=as.numeric(as.factor(region_50))+dodgi,fill=region_50,color=region_50),height = 0,size=1)+
theme_bw()+geom_vline(xintercept=0,linetype="dashed")+
theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.border=element_blank(),
panel.grid.minor = element_blank(),panel.background = element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),
strip.background=element_rect(fill=NA,color=NA),legend.position="none")+
scale_color_manual(values=colo2)+scale_fill_manual(values=colo2)+ylab("Bio-geographic region")+xlab(expression(paste("Species trend (",log(tau[r]),")")))+coord_cartesian(expand=FALSE)+ggtitle("a")


pl2=ggplot()+
geom_density_ridges(data=bidon,aes(x=trend,y=region_50,fill=region_50,color=region_50),alpha=0.3,scale = 0.9)+
xlim(c(-0.2,0.2))+
geom_point(data=b,aes(x=emmean,y=as.numeric(as.factor(region_50))+dodgi,fill=region_50,color=region_50),size=2)+
geom_errorbarh(data=b,aes(x=emmean,xmax=emmean+1.96*SE,xmin=emmean-1.96*SE,y=as.numeric(as.factor(region_50))+dodgi,fill=region_50,color=region_50),height = 0,size=1)+
theme_bw()+geom_vline(xintercept=0,linetype="dashed")+
theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.border=element_blank(),
panel.grid.minor = element_blank(),panel.background = element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),
strip.background=element_rect(fill=NA,color=NA),legend.position="none")+
scale_color_manual(values=colo2)+scale_fill_manual(values=colo2)+ylab("Bio-geographic region")+xlab(expression(paste("Species trend (",log(tau[r]),")")))+coord_cartesian(expand=FALSE)+ggtitle("b")


####### BASELINE shift
bidon=subset(trendsf,convergence==0 & !is.na(acim) & baseline2>1921 & !is.na(sde))
nb_weird=nrow(subset(bidon,abs(trend)>=1))
bidon=subset(bidon,abs(trend)<1)

# for(j in 2:length(baselines_vec)){
# bidon2=subset(bidon,convergence==0 & !is.na(acim) & baseline2==baselines_vec[j] & !is.na(sde))

# model=rma(trend ~ region_50, sei=sde,data=bidon2, digits=10)
# save(model,file=paste0(project_folder,"data/model_",baselines_vec[j],".RData"))
# }

bf=b
for(j in 2:length(baselines_vec)){
bidon2=subset(bidon,convergence==0 & !is.na(acim) & baseline2==baselines_vec[j] & !is.na(sde))
load(paste0(project_folder,"data/model_",baselines_vec[j],".RData"))
boxplot(weights(model))
sav <- emmprep(model)
bb=as.data.frame(emmeans(sav,specs="region_50"))
bb$baseline=baselines_vec[j]
bb$basic_mean=bidon2 %>% group_by(region_50) %>% summarise(mean(trend)) %>% deframe()
bb$basic_SE=bidon2 %>% group_by(region_50) %>% summarise(sd(trend)/sqrt(length(trend))) %>% deframe()
bf=rbind(bf,bb)
}

pl3=ggplot()+
geom_point(data=bf,aes(y=emmean,x=baseline,fill=region_50,color=region_50),position=position_dodge(width=2),size=2)+
geom_errorbar(data=bf,aes(y=emmean,ymax=emmean+1.96*SE,ymin=emmean-1.96*SE,x=baseline,fill=region_50,color=region_50),position=position_dodge(width=2),width=0,size=1)+
theme_bw()+geom_hline(yintercept=0,linetype="dashed")+
theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.border=element_blank(),
panel.grid.minor = element_blank(),panel.background = element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),
strip.background=element_rect(fill=NA,color=NA),legend.position="none")+
scale_color_manual(values=colo2)+scale_fill_manual(values=colo2)+xlab("Baseline year")+ylab(expression(paste("Species trend (",log(tau[r]),")")))+ggtitle("c")

pl4=ggplot()+
geom_point(data=bf,aes(y=emmean,x=baseline,fill=region_50,color=region_50),position=position_dodge(width=2),size=2)+
geom_errorbar(data=bf,aes(y=emmean,ymax=emmean+1.96*SE,ymin=emmean-1.96*SE,x=baseline,fill=region_50,color=region_50),position=position_dodge(width=2),width=0,size=1)+
theme_bw()+geom_hline(yintercept=0,linetype="dashed")+
theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.border=element_blank(),
panel.grid.minor = element_blank(),panel.background = element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),
strip.background=element_rect(fill=NA,color=NA),legend.position="none")+
scale_color_manual(values=colo2)+scale_fill_manual(values=colo2)+xlab("Baseline year")+ylab(expression(paste("Species trend (",log(tau[r]),")")))+ggtitle("c")

plot_grid(pl1,pl3,pl2,pl4,ncol=2,align="hv")

pdf(paste0(project_folder,"Figure_2.pdf"),width=8,height=6)
plot_grid(pl1,pl3,pl2,pl4,ncol=2,align="hv",rel_widths=c(1.2,1))
dev.off();

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

model=lmer(trend~cate2+(1|species),data=bidon,weights=1/(sde^2))


