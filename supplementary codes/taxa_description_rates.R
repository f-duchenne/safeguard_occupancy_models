pkgs <- c("data.table", "dplyr","lme4","ggplot2","ggridges","metafor","cowplot","emmeans","tidyverse","gridExtra") 
inst <- pkgs %in% installed.packages()
if (any(inst)) install.packages(pkgs[!inst])
pkg_out <- lapply(pkgs, require, character.only = TRUE)
#project_folder="C:/Users/Duchenne/Documents/safeguard/"
project_folder <- ""

colo2=c("#44AA99","#117733","#332288","#CC6677","#DDCC77")

#Load Our Sup Mat with species per region
dat <- fread(file = paste0(project_folder,"data/final_and_intermediate_outputs/database_clean_filtered.csv"))
head(dat)

#### nb records species per bioregion
sm=dat %>% group_by(scientificName,region_50,taxo_group,scientificNameAuthorship) %>% count()

#extract year of description
sm$year_description=as.numeric(unlist(regmatches(sm$scientificNameAuthorship, gregexpr("[[:digit:]]+", sm$scientificNameAuthorship))))
sm$val=1

des_reg= sm %>% group_by(region_50,taxo_group,year_description) %>% summarise(nb=sum(val))
des_reg= des_reg %>% arrange(year_description) %>% group_by(region_50,taxo_group) %>% mutate(cumul=cumsum(nb))

###3 add final values for all region
des_reg= des_reg %>% group_by(region_50,taxo_group) %>% mutate(max_year=max(year_description))

final_state=subset(des_reg,year_description==max_year)
final_state$year_description=2023

des_reg=rbind(des_reg,final_state)

des_reg=subset(des_reg,region_50 %in% c("alpine","boreal","atlantic","continental","mediterranean"))


ggplot(data=des_reg,aes(x=year_description,y=cumul,color=region_50))+geom_line(size=1.3)+theme_bw()+
theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.border=element_blank(),
panel.grid.minor = element_blank(),panel.background = element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),
strip.background=element_rect(fill=NA,color=NA),legend.position="right")+
scale_color_manual(values=colo2)+scale_fill_manual(values=colo2)+ylab("Number of described species")+xlab("Years")+facet_wrap(~taxo_group)+
labs(color="Biogeographic region")

