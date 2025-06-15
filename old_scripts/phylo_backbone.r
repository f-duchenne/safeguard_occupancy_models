#Calculate phylogenetic matrix for Europe species
library(phangorn)
library(ape)
library(phytools)
library(dplyr)
library(data.table)
library(stringr)
setwd(dir="C:/Users/Duchenne/Documents/safeguard/data")
#Load extracted data
pref <- fread("liste_total_species_occupancy.csv")

#Check levels of species and cover names
#First, Species
species = subset(pref,taxo_group=="bees")$species 
genus=sapply(str_split(species," "),"[",1)

#-----
#PHYLO
bee.tree=read.tree(file="BEE_mat7_fulltree_tplo35_sf20lp.nwk")
bee.tree=force.ultrametric(bee.tree)
bee.tree$tip.label=gsub("~","_",bee.tree$tip.label)
genus_of_tips=sapply(str_split(bee.tree$tip.label,"_"),"[",1)
bee.tree$tip.label=paste0(sapply(str_split(bee.tree$tip.label,"_"),"[",1),"_",sapply(str_split(bee.tree$tip.label,"_"),"[",2))
drop.tip(bee.tree,bee.tree$tip.label[!(genus_of_tips %in% genus)])

# bee.tree=read.tree(file="BEE_mat7_fulltree.nwk")
# bee.tree$tip.label=gsub("~","_",bee.tree$tip.label)
# bee.tree$tip.label=sapply(str_split(bee.tree$tip.label,"_"),"[",1)

#creates Seladonia genus:
nodi=getMRCA(bee.tree, bee.tree$tip.label[grep("Halictus",bee.tree$tip.label)])
pos=bee.tree$edge.length[nodi==bee.tree$edge[,2]]/2
bee.tree=bind.tip(bee.tree,"Seladonia_sp",where=Ancestors(bee.tree, nodi, type = c("parent")),posistion=pos)

species[!(species %in% str_replace_all(bee.tree$tip.label,"_"," "))]
# plot(bee.tree)
# nodelabels()
# tiplabels()

#add dummy species labels
bee.tree$tip.label<-paste(bee.tree$tip.label," dum",sep="")

#Add species tips
for(i in 1:length(species)){
    bee.tree<-add.species.to.genus(bee.tree,species[i],
                                      where="root")
}

## prune out dummy taxa
ii<-grep("dum",bee.tree$tip.label)
bee.tree<-drop.tip(bee.tree,bee.tree$tip.label[ii])
#Our tree
plot(bee.tree, cex = 0.6)

##Check for missing species
setdiff(species,str_replace_all(bee.tree$tip.label,"_"," "))
write.tree(bee.tree,"phylo_bees_safeguard.tree")


setwd(dir="C:/Users/Duchenne/Documents/safeguard/data")
#Check levels of species and cover names
#First, Species
species = subset(pref,taxo_group=="hoverflies")$species 
genus=sapply(str_split(species," "),"[",1)

#-----
#PHYLo substitutions:
hov.tree2=read.tree(file="syrphidae_phylo.tree")
hov.tree2$tip.label=gsub("'","",hov.tree2$tip.label,fixed=TRUE)
hov.tree2$tip.label=gsub("*","",hov.tree2$tip.label,fixed=TRUE)
ages=as.data.frame(tree.age(hov.tree2))
calib=ages[,c(2,1)]
calib$age.max=calib[,2]
calib$soft=FALSE
names(calib)=c("node","age.min","age.max","soft.bounds")
calib=calib[244:nrow(calib),]
calib$node=as.numeric(calib$node)
calib=calib[!duplicated(calib$age.min),]
calib=subset(calib,age.min>1)
calib=fread("calibration.csv")
names(calib)=c("node","age.min","age.max","soft.bounds")
hov.tree3=chronos(hov.tree2,calibration=calib)

hov.tree3$tip.label

#add dummy species labels
hov.tree3$tip.label=gsub(" ","_",hov.tree3$tip.label,fixed=TRUE)
hov.tree3$tip.label<-paste(hov.tree3$tip.label,"_dum",sep="")

#Add species tips
for(i in 1:length(species)){
    hov.tree3<-add.species.to.genus(hov.tree3,species[i],
                                      where="root")
}

## prune out dummy taxa
ii<-grep("_dum",hov.tree3$tip.label)
hov.tree3<-drop.tip(hov.tree3,hov.tree3$tip.label[ii])
#Our tree
plot(hov.tree3, cex = 0.6)

##Check for missing species
setdiff(species,str_replace_all(hov.tree3$tip.label, "_", " "))
write.tree(hov.tree3,"phylo_hovs_safeguard.tree")











hov.tree=read.nexus(file="syrphidae_phylo_time_calibrated.nexus")
hov.tree$tip.label=gsub("'","",hov.tree$tip.label,fixed=TRUE)
#rename tips with species names
access=fread("accession_syrphidae_phylogeny.csv")
access$label=gsub(" ","_",access$label,fixed=TRUE)
hov.tree$tip.label[match(access$label,hov.tree$tip.label)]=access$species
