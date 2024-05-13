###########################################
###########################################
#' Check for packages and if necessary install into library 
#+ message = FALSE
rm(list=ls())
pkgs <- c("data.table", "dplyr","R2jags","ggplot2") 

inst <- pkgs %in% installed.packages()
if (any(inst)) install.packages(pkgs[!inst])
pkg.out <- lapply(pkgs, require, character.only = TRUE)

#defining working folder:
setwd(dir="C:/Users/Duchenne/Documents/safeguard/data")
resf=NULL
for(i in 1:10){

load(paste0("model_",i,".RData"))

suma=as.data.frame(lili[[1]]$BUGSoutput$summary)
suma$vari=rownames(suma)
names(suma)[1:7]=c("est","sd","low","first_quartile","median","third_quartile","high")
suma=suma[grep("site.eff",suma$vari,invert=T),]
suma$type=suma$vari
suma$type[grep("det",suma$vari)]="detection"
suma$type[grep("occ",suma$vari)]="occupancy"
suma$country.num=sapply(strsplit(gsub("]","",sapply(strsplit(suma$vari,"[",fixed=T),function(x){x[2]}),fixed=T),","),function(x){x[1]})
suma$period.num=sapply(strsplit(gsub("]","",sapply(strsplit(suma$vari,"[",fixed=T),function(x){x[2]}),fixed=T),","),function(x){x[2]})

suma=merge(suma,lili[[2]],by="country.num")
suma=merge(suma,lili[[3]],by="period.num")
suma$species=lili[[4]]

resf=rbind(resf,suma)
}



ggplot(data=subset(resf,type=="occupancy"),aes(x=time_period,y=est,color=COUNTRY,linetype=COUNTRY))+geom_pointrange(aes(ymin=low,ymax=high))+geom_line()+
theme_bw()+theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
panel.border = element_blank(),panel.background = element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),
strip.background=element_rect(fill=NA,color=NA))+
facet_wrap(~species,scales="free")+xlab("Time period")+scale_x_continuous(breaks=c(1980,1990,2000,2010,2020),labels=c("<1980","1981-1990","1991-2000","2001-2010","2011-2020"))+
ylab("Occupancy (logit scale)")