pkgs <- c("data.table", "dplyr","lme4","ggplot2","ggridges","metafor","cowplot","emmeans","tidyverse","gridExtra") 
inst <- pkgs %in% installed.packages()
if (any(inst)) install.packages(pkgs[!inst])
pkg_out <- lapply(pkgs, require, character.only = TRUE)
project_folder="C:/Users/Duchenne/Documents/safeguard/"


colo2=c("#44AA99","#117733","#332288","#CC6677","#DDCC77")

#Load Our Sup Mat with species per region
dat <- fread(file = paste0(project_folder,"data/database_clean_filtered.csv"))
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



################# PART THAT WAS CODED BY NACHO INITIALLY:

#CARLOS: I think its better if we have a list of species present per biogeoregion, as here I am only using species with trends.

#select Med species (Trend_mediterranean)
med <- sm[which(!is.na(sm$Trend_mediterranean)),]
head(med)
#select Cont species (Trend_atlantic)
atl <- sm[which(!is.na(sm$Trend_atlantic)),]
#select Med species (Trend_boreal)
bor <- sm[which(!is.na(sm$Trend_boreal)),]
#select Med species (Trend_continental)
cont <- sm[which(!is.na(sm$Trend_continental)),]
#select Med species (Trend_alpine)
alp <- sm[which(!is.na(sm$Trend_alpine)),]

#subset from master, med species ----
d_med <- subset(master, friendly_name %in% med$TAXON) #Assuming taxonomy matches!!
head(d_med)

med_year <- c()
matches <- regmatches(d_med$taxonomic_authority, gregexpr("[[:digit:]]+", d_med$taxonomic_authority))
med_year <- as.numeric(unlist(matches))

#subset from master, continental species ----
d_cont <- subset(master, friendly_name %in% cont$TAXON) #Assuming taxonomy matches!!

cont_year <- c()
matches <- regmatches(d_cont$taxonomic_authority, gregexpr("[[:digit:]]+", d_cont$taxonomic_authority))
cont_year <- as.numeric(unlist(matches))

#subset from master, alp species----
d_alp <- subset(master, friendly_name %in% alp$TAXON) #Assuming taxonomy matches!!

alp_year <- c()
matches <- regmatches(d_alp$taxonomic_authority, gregexpr("[[:digit:]]+", d_alp$taxonomic_authority))
alp_year <- as.numeric(unlist(matches))

#subset from master, bor species----
d_bor <- subset(master, friendly_name %in% bor$TAXON) #Assuming taxonomy matches!!

bor_year <- c()
matches <- regmatches(d_bor$taxonomic_authority, gregexpr("[[:digit:]]+", d_bor$taxonomic_authority))
bor_year <- as.numeric(unlist(matches))

#subset from master, atl species----
d_atl <- subset(master, friendly_name %in% atl$TAXON) #Assuming taxonomy matches!!

atl_year <- c()
matches <- regmatches(d_atl$taxonomic_authority, gregexpr("[[:digit:]]+", d_atl$taxonomic_authority))
atl_year <- as.numeric(unlist(matches))

#Put all lines toguether

plot(as.numeric(cumsum(table(med_year))) ~ as.numeric(names(table(med_year))),
     las = 1, 
     type = "l", 
     ylab = "Species Described", 
     xlab = "Year", 
     col = 1, 
     main = "Cumulative Species Descriptions Over Time")
lines(as.numeric(cumsum(table(bor_year))) ~ as.numeric(names(table(bor_year))),
      col = 2)
lines(as.numeric(cumsum(table(alp_year))) ~ as.numeric(names(table(alp_year))),
      col = 3)
lines(as.numeric(cumsum(table(cont_year))) ~ as.numeric(names(table(cont_year))),
      col = 4)
lines(as.numeric(cumsum(table(atl_year))) ~ as.numeric(names(table(atl_year))),
      col = 5)

#CARLOS: Match colors with the rest of the ms.
