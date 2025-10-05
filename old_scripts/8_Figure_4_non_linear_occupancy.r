############# FIRST ASSEMBLE ALL SPECIES TOGETHER
pkgs <- c("data.table", "dplyr","car") 
inst <- pkgs %in% installed.packages()
if (any(inst)) install.packages(pkgs[!inst])
pkg_out <- lapply(pkgs, require, character.only = TRUE)

project_folder="C:/Users/Duchenne/Documents/safeguard/"


lifile=list.files(path =paste0(project_folder,"predicts_non_linear/"), pattern ="predicts_")

trendsf=NULL
for(i in lifile){
index=gsub(".RData","",gsub("predicts_","",i))
load(paste0(project_folder,"predicts_non_linear/",i))
trends=lili[[1]]
trendsf=rbind(trendsf,trends)
}

fwrite(trendsf,paste0(project_folder,"data/all_trends_non_linear.csv"))


############# SECOND PLOT RESULTS:
pkgs <- c("data.table", "dplyr","lme4","ggplot2","ggridges","metafor","cowplot","emmeans","tidyverse","car") 
inst <- pkgs %in% installed.packages()
if (any(inst)) install.packages(pkgs[!inst])
pkg_out <- lapply(pkgs, require, character.only = TRUE)
project_folder="C:/Users/Duchenne/Documents/safeguard/"

colo2=c("#44AA99","#117733","#332288","#CC6677","#DDCC77")

inv.logit=function(x){exp(x)/(1+exp(x))}
logit=function(x){log(x/(1-x))}

trendsf=fread(paste0(project_folder,"data/all_trends_non_linear.csv"))
trendsf$x=as.character(trendsf$x)
trendsf$x=gsub("(","",trendsf$x,fixed=TRUE)
trendsf$x=gsub(")","",trendsf$x,fixed=TRUE)
trendsf$year=1971+as.numeric(trendsf$x)-1
trendsf=trendsf %>% group_by(species,group) %>% mutate(ner=sum(is.na(std.error)))
trendsf=subset(trendsf,ner==0 & !is.na(acim) & std.error<5)

bidon=trendsf %>% dplyr::group_by(year,taxo_group,group,species) %>% dplyr::summarise(moy=mean(pred))

ggplot(data=bidon,aes(x=year,y=moy,color=group,group=species))+geom_line()+
facet_grid(cols=vars(taxo_group))

inv.logit=function(x){exp(x)/(1+exp(x))}
resf=NULL
for(jjj in unique(trendsf$group)){
	for(jj in unique(trendsf$taxo_group)){
		for(j in 1971:2020){
			bidon2=subset(trendsf,year==j & taxo_group==jj & group==jjj)
			func=function(x,y){rnorm(1000,x,y)}
			res=data.frame(mean_occupancy=apply(inv.logit(t(mapply(func,x=car::logit(bidon2$pred,adjust=1e-15),y=bidon2$std.error))),2,mean),year=j,taxo_group=jj,region=jjj)
			resf=rbind(resf,res)
		}
	}
}


bidon=resf %>% dplyr::group_by(year,taxo_group,region) %>% dplyr::summarise(moy=mean(mean_occupancy,na.rm=TRUE),lcl=quantile(mean_occupancy,probs=c(0.025),na.rm=TRUE),ucl=quantile(mean_occupancy,probs=c(0.975),na.rm=TRUE))

ggplot(data=bidon,aes(x=year,y=moy,color=region))+geom_line()+
facet_grid(cols=vars(taxo_group))