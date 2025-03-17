###########################################
###########################################
#' Check for packages and if necessary install into library 
#+ message = FALSE
rm(list=ls())
pkgs <- c("data.table", "dplyr","rgbif") 

inst <- pkgs %in% installed.packages()
if (any(inst)) install.packages(pkgs[!inst])
pkg.out <- lapply(pkgs, require, character.only = TRUE)

setwd(dir="C:/Users/Duchenne/Documents/safeguard/data")
dat=fread("Species_absent_in_P5.csv")

dat$nbocc=NA
for(i in 1:nrow(dat)){
dat$nbocc[i]=occ_count(scientificName=dat$V2[i],year= '2001,2021')
}