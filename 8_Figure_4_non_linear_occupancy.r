############# FIRST ASSEMBLE ALL SPECIES TOGETHER
pkgs <- c("data.table", "dplyr") 
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
pkgs <- c("data.table", "dplyr","lme4","ggplot2","ggridges","metafor","cowplot","emmeans","tidyverse") 
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
trendsf=subset(trendsf,ner==0)

bidon=trendsf %>% dplyr::group_by(group,year,taxo_group) %>% dplyr::summarise(moy=mean(pred))