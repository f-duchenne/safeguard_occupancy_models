#Calculate phylogenetic matrix for Europe species
library(visreg) 
library(phangorn)
library(ape)
library(phytools)
library(dplyr)
library(readr)
library(stringr)
setwd(dir="C:/Users/Duchenne/Documents/safeguard/data")
#Load extracted data
pref <- read_csv("liste_total_species_occupancy.csv")

#Check levels of species and cover names
#First, Species
species = subset(pref,taxo_group=="bees")$species 
genus=sapply(str_split(species," "),"[",1)

#-----
#PHYLO
bee.tree=read.tree(file="BEE_mat7gen_p8pmAa_fst.nwk")
bee.tree$tip.label=gsub("~","_",bee.tree$tip.label)
bee.tree$tip.label=sapply(str_split(bee.tree$tip.label,"_"),"[",1)
drop.tip(bee.tree,bee.tree$tip.label[!(bee.tree$tip.label %in%genus)])
bee.tree=chronos(bee.tree)
# bee.tree=read.tree(file="BEE_mat7_fulltree.nwk")
# bee.tree$tip.label=gsub("~","_",bee.tree$tip.label)
# bee.tree$tip.label=sapply(str_split(bee.tree$tip.label,"_"),"[",1)


#creates Seladonia genus:
nodi=which("Halictus"==bee.tree$tip.label)
pos=bee.tree$edge.length[nodi==bee.tree$edge[,2]]/2
bee.tree=bind.tip(bee.tree,"Seladonia",where=Ancestors(bee.tree, nodi, type = c("parent")),posistion=pos)




bee.tree$tip.label
species

plot(bee.tree)
nodelabels()
tiplabels()


#add dummy species labels
bee.tree$tip.label<-paste(bee.tree$tip.label,"_dum",sep="")

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
setdiff(species,str_replace_all(bee.tree$tip.label, "_", " "))
write.tree(bee.tree,"phylo_bees_safeguard.nxs")