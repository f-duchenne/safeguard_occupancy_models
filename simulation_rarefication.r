########################### NON STRATIFIED SAMPLING to use iNEXT with "abundance"
library(dplyr)
library(iNEXT)
library(reshape2)
library(gridExtra)
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
# a pool of species with a detection probability for each
sp_pool=rbind(data.frame(sp=c(letters[1:26],LETTERS[1:26]),prob=exp(-runif(52,0,10)),site="site1"),
data.frame(sp=c(letters[1:26],LETTERS[1:26]),prob=exp(-runif(52,0,3)),site="site2"),
data.frame(sp=c(letters[1:26],LETTERS[1:26]),prob=runif(52,0,1),site="site3"))
sp_pool=sp_pool %>% group_by(site) %>% mutate(prob_rel=prob/sum(prob))
sp_pool=sp_pool %>% group_by(site) %>% mutate(rich=sum(prob_rel>1e-2))
hist(sp_pool$prob)
##put manually two rares s

# a vector of different sampling effort
samp_vec =c(1,2,3)
samp_events=c(5,15,50)

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
					resf=rbind(resf,data.frame(sp=rep(pool$sp,vec),samp=ii,nbsamp=i,samp_event=z,site=si,trial=j,rich=unique(pool$rich)))
				}
			}
		}
	}
}

b=resf %>% group_by(samp,site,trial,sp,rich,nbsamp,samp_event) %>% count(name ="abund")
b$id=paste(b$site,b$trial,b$nbsamp,b$samp,sep=".")
comm=split(b, b$id)

comm=lapply(comm,function(x){dcast(x,id+sp~samp_event,value.var="abund",fill=0)})
comm=lapply(comm,function(x){x[,-c(1:2)]})

bidon=unique(b[,c("id","site","trial","rich","nbsamp")])

rari=iNEXT(comm,q=0,size=c(50,100,500,1000),datatype ="incidence_raw")
str(rari)
str(rari$iNextEst)
str(rari$AsyEst)

rich_out=subset(rari$AsyEst,Diversity=="Species richness")
rich_out=merge(rich_out,bidon,by.x="Assemblage",by.y="id")

pl1=ggplot(data=sp_pool,aes(x=prob))+geom_histogram()+xlab("detection probability (abundance)")+ylab("Number of species")+
theme_bw()+theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),
strip.background=element_rect(fill=NA,color=NA),panel.border=element_blank(),legend.position="none")+
facet_wrap(~site,scale="free_x")

pl2=ggplot(data=rich_out,aes(x=as.factor(nbsamp),y=Estimator,color=samp))+geom_boxplot()+xlab("sampling pressure (number of bees observed)")+ylab("Richness estimation")+
theme_bw()+theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),
strip.background=element_rect(fill=NA,color=NA),panel.border=element_blank(),legend.position="none")+
facet_wrap(~site)+geom_hline(aes(yintercept=52))

grid.arrange(pl1,pl2)