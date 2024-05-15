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
setwd(dir="C:/Users/Duchenne/Documents/safeguard/models")
resf=NULL
for(i in c(1,2,3,5)){

load(paste0("model_",i,".RData"))

suma=as.data.frame(lili[[1]]$BUGSoutput$summary)
suma$vari=rownames(suma)
names(suma)[1:7]=c("est","sd","low","first_quartile","median","third_quartile","high")
suma=suma[grep("site.eff",suma$vari,invert=T),]
suma$type=suma$vari
suma$type[grep("det",suma$vari)]="detection"
suma$type[grep("period.occ",suma$vari)]="occupancy"
suma$type[grep("eu_eff",suma$vari)]="eu_eff"
suma$country.num=sapply(strsplit(gsub("]","",sapply(strsplit(suma$vari,"[",fixed=T),function(x){x[2]}),fixed=T),","),function(x){x[1]})
suma$period.num=sapply(strsplit(gsub("]","",sapply(strsplit(suma$vari,"[",fixed=T),function(x){x[2]}),fixed=T),","),function(x){x[2]})
suma$period.num[suma$type=="eu_eff"]=suma$country.num[suma$type=="eu_eff"]
suma$country.num[suma$type=="eu_eff"]=NA

suma=merge(suma,lili[[2]],by="country.num",all.x=TRUE,all.y=FALSE)
suma=merge(suma,lili[[3]],by="period.num",all.x=TRUE,all.y=FALSE)
suma$species=lili[[4]]

resf=rbind(resf,suma)
}



ggplot()+geom_pointrange(data=subset(resf,type=="occupancy"),aes(x=time_period,y=est,color=COUNTRY2,ymin=low,ymax=high))+
geom_line(data=subset(resf,type=="occupancy"),aes(x=time_period,y=est,color=COUNTRY2),alpha=0.5)+
geom_pointrange(data=subset(resf,type=="eu_eff"),aes(x=time_period,y=est,ymin=low,ymax=high),col="black")+
geom_line(data=subset(resf,type=="eu_eff"),aes(x=time_period,y=est),col="black",size=1.3)+
theme_bw()+theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
panel.border = element_blank(),panel.background = element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),
strip.background=element_rect(fill=NA,color=NA))+
facet_wrap(~species,scales="free")+xlab("Time period")+scale_x_continuous(breaks=c(1980,1990,2000,2010,2020),labels=c("<1980","1981-1990","1991-2000","2001-2010","2011-2020"))+
ylab("Occupancy (logit scale)")+labs(colour="Country")