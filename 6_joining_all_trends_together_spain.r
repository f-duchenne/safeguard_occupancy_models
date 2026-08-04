#DO NOT RUN#
#This script is prepared to be run in a cluster, as it contains computationally demanding models.

############# FIRST ASSEMBLE ALL SPECIES TOGETHER
pkgs <- c("data.table", "dplyr") 
pkg_out <- lapply(pkgs, require, character.only = TRUE)

#project_folder="C:/Users/Duchenne/Documents/safeguard/"
project_folder <- ""

lifile=list.files(path="results/predicts_spain/", pattern ="predicts_",full.names =TRUE)

trendsf=NULL
for(i in lifile){
	index=gsub(".RData","",gsub("predicts_","",i))
	load(i)
	trends=lili[[1]]
	trendsf=rbind(trendsf,trends)
}

inv.logit=function(x){exp(x)/(1+exp(x))}
logit=function(x){log(x/(1-x))}


trendsf$significant="no"
trendsf$significant[trendsf$trend<0 & (trendsf$trend+1.96*trendsf$sde)<0]="yes"
trendsf$significant[trendsf$trend>0 & (trendsf$trend-1.96*trendsf$sde)>0]="yes"

trendsf$growth.rate=100*(exp(trendsf$trend)-1)

trash=subset(trendsf,(convergence!=0 | is.na(acim) | is.na(sde)))
trash %>% group_by(taxo_group,baseline) %>% summarise(length(unique(species)))


length(unique(trendsf$species[trendsf$baseline==1921]))
trendsf=subset(trendsf,convergence==0 & !is.na(acim) & !is.na(sde))
nrow(subset(trendsf,baseline==1921))
nrow(subset(trendsf,baseline==1971))
length(unique(trendsf$species[trendsf$baseline==1921]))
length(unique(trendsf$species[trendsf$baseline==1971]))

fwrite(trendsf,paste0(project_folder,"data/final_and_intermediate_outputs/all_trends_spain.csv"))
