########################### NON STRATIFIED SAMPLING to use iNEXT with "abundance"
library(dplyr)
library(iNEXT)
library(reshape2)
library(gridExtra)
library(ggplot2)
# a pool of species with a detection probability for each
sp_pool=rbind(data.frame(sp=c(letters[1:26],LETTERS[1:26]),prob=exp(-runif(52,0,10)),site="site1"),
data.frame(sp=c(letters[1:26],LETTERS[1:26]),prob=exp(-runif(52,0,3)),site="site2"),
data.frame(sp=c(letters[1:26],LETTERS[1:26]),prob=runif(52,0,1),site="site3"))
sp_pool=sp_pool %>% group_by(site) %>% mutate(prob_rel=prob/sum(prob))
sp_pool=sp_pool %>% group_by(site) %>% mutate(rich=sum(prob_rel>1e-4))
hist(sp_pool$prob)
##put manually two rares s

# a vector of different sampling effort
samp_vec = c(50,100,500)

# generate the records
nbtrial=10
resf=NULL
for(si in unique(sp_pool$site)){
	pool=subset(sp_pool,site==si)
	for(i in samp_vec){
		for(j in 1:nbtrial){
			resf=rbind(resf,data.frame(sp=sample(pool$sp,i,replace = TRUE, prob = pool$prob_rel),samp=i,site=si,trial=j,rich=unique(pool$rich)))
		}
	}	
}

b=resf %>% group_by(samp,site,trial,sp,rich) %>% count(name ="abund")
b$id=paste(b$samp,b$site,b$trial,sep=".")

comm=split(b$abund, b$id)


bidon=unique(b[,c("id","site","samp","trial","rich")])

rari=iNEXT(comm,q=0,size=c(50,100,500,1000))
str(rari)
str(rari$iNextEst)
str(rari$AsyEst)

rich_out=subset(rari$AsyEst,Diversity=="Species richness")
rich_out=merge(rich_out,bidon,by.x="Assemblage",by.y="id")

pl1=ggplot(data=sp_pool,aes(x=prob_rel))+geom_histogram()+xlab("detection probability (abundance)")+ylab("Number of species")+
theme_bw()+theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),
strip.background=element_rect(fill=NA,color=NA),panel.border=element_blank(),legend.position="none")+
facet_wrap(~site,scale="free_x")

pl2=ggplot(data=rich_out,aes(x=as.factor(samp),y=Estimator))+geom_boxplot()+xlab("sampling pressure (number of bees observed)")+ylab("Richness estimation")+
theme_bw()+theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),
strip.background=element_rect(fill=NA,color=NA),panel.border=element_blank(),legend.position="none")+
facet_wrap(~site)+geom_hline(aes(yintercept=rich))

grid.arrange(pl1,pl2)

rich_out=subset(rari$iNextEst$coverage_based,Method=="Extrapolation")
rich_out=merge(rich_out,bidon,by.x="Assemblage",by.y="id")

pl1=ggplot(data=sp_pool,aes(x=prob))+geom_histogram()+xlab("detection probability (abundance)")+ylab("Number of species")+
theme_bw()+theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),
strip.background=element_rect(fill=NA,color=NA),panel.border=element_blank(),legend.position="none")+
facet_wrap(~site)

pl2=ggplot(data=rich_out,aes(x=as.factor(samp),y=qD,color=as.factor(m)))+geom_boxplot()+xlab("sampling pressure (number of bees observed)")+ylab("Richness estimation")+
theme_bw()+theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),
strip.background=element_rect(fill=NA,color=NA),panel.border=element_blank(),legend.position="none")+
facet_wrap(~site)

grid.arrange(pl1,pl2)

########################### STRATIFIED SAMPLING to use iNEXT with "incidence_raw"
library(dplyr)
library(iNEXT)
library(reshape2)
library(gridExtra)
library(ggplot2)
nbsp=100
# a pool of species with a detection probability for each
sp_pool_site1=data.frame(sp=paste0("sp",1:nbsp),prob=exp(-runif(nbsp,0,10)),site="site1")
sp_pool_site2=data.frame(sp=paste0("sp",1:nbsp),prob=exp(-runif(nbsp,0,3)),site="site2")
sp_pool_site3=data.frame(sp=paste0("sp",1:nbsp),prob=runif(nbsp,0,1),site="site3")

#sp_pool=sp_pool %>% group_by(site) %>% mutate(prob=prob/sum(prob))

sp_pool=rbind(sp_pool_site1,sp_pool_site2,sp_pool_site3) #put all sites together

pl1=ggplot(data=sp_pool,aes(x=prob))+geom_histogram(fill="white",color="black")+xlab("detection probability (abundance)")+ylab("Number of species")+
theme_bw()+theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),
strip.background=element_rect(fill=NA,color=NA),panel.border=element_blank(),legend.position="none")+
facet_wrap(~site)

pl1

hist(sp_pool$prob)
##put manually two rares s

# a vector of different sampling effort
samp_vec =c(1,2,3) #sampling intensity
samp_events=c(5,10,20) #number of sampling events

# generate the records
nbtrial=10
resf=NULL
for(si in unique(sp_pool$site)){
	pool=subset(sp_pool,site==si)
	for(ii in samp_vec){
		for(i in samp_events){
			for(z in 1:i){
				for(j in 1:nbtrial){
					vec=c(sapply(pool$prob,function(x){rbinom(1,ii,x)}))
					resf=rbind(resf,data.frame(sp=rep(pool$sp,vec),samp=ii,nbsamp=i,samp_event=z,site=si,trial=j))
				}
			}
		}
	}
}


b=resf %>% group_by(samp,site,trial,sp,nbsamp,samp_event) %>% count(name ="abund") #aggregate data per species, for each site and simulation parameters

b$id=paste(b$site,b$trial,b$nbsamp,b$samp,sep=".") #create a unique id for each simulations

comm=split(b, b$id) #creat a list fo community data, one dataframe per simulations
comm=lapply(comm,function(x){dcast(x,id+sp~samp_event,value.var="abund",fill=0)}) #put each data frame in the right format for iNEXT function
comm=lapply(comm,function(x){x[,-c(1:2)]}) #remove id columns

rari=iNEXT(comm,q=0,size=c(50,100,500,1000),datatype ="incidence_raw")
str(rari)
str(rari$iNextEst)
str(rari$AsyEst)


rich_out=subset(rari$AsyEst,Diversity=="Species richness") #extract the predicted richness
rich_out=merge(rich_out,unique(b[,c("id","site","samp","trial","nbsamp")]),by.x="Assemblage",by.y="id") #adding information about parameters of simulations to that

subset(rich_out,samp==1 & nbsamp==5 & site==3)

pl1=ggplot(data=sp_pool,aes(x=prob))+geom_histogram()+xlab("detection probability (abundance)")+ylab("Number of species")+
theme_bw()+theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),
strip.background=element_rect(fill=NA,color=NA),panel.border=element_blank(),legend.position="none")+
facet_wrap(~site,scale="free_x")

pl2=ggplot(data=rich_out,aes(x=as.factor(nbsamp),y=Estimator,color=as.factor(samp)))+geom_boxplot()+xlab("Number of sampling events")+ylab("Richness estimation")+
theme_bw()+theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),
strip.background=element_rect(fill=NA,color=NA),panel.border=element_blank(),legend.position="none")+
facet_wrap(~site)+geom_hline(aes(yintercept=100))


grid.arrange(pl1,pl2)

ggplot(data=rich_out,aes(x=nbsamp,y=Estimator))+geom_point()+
stat_smooth()+xlab("Number of sampling events")+ylab("Richness estimation")+
scale_y_log10()+
theme_bw()+theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),
strip.background=element_rect(fill=NA,color=NA),panel.border=element_blank(),legend.position="none")+
facet_wrap(~site,scales="free")+geom_hline(aes(yintercept=52))


