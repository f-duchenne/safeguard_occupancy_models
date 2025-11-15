#this script explores the database for descriptive analysis.
pkgs <- c("data.table", "dplyr","ggplot2","cowplot","gridExtra", "sf") 
inst <- pkgs %in% installed.packages()
if (any(inst)) install.packages(pkgs[!inst])
ypkg_out <- lapply(pkgs, require, character.only = TRUE)

#colors for regions:
colo2=c("#44AA99","#117733","#332288","#CC6677","#DDCC77")

#defining working folder:
#project_folder="C:/Users/Duchenne/Documents/safeguard/"
project_folder=""

# Loading data
datf=fread(paste0(project_folder,"data/final_and_intermediate_outputs/database_clean_filtered.csv"))
datf %>% group_by(taxo_group) %>% count()
datf %>% group_by(taxo_group) %>% summarise(length(unique(scientificName)))

#define the total number of species expected, accoridng to Reverte:
datf$ntot_spec_europe_tax=2138
datf$ntot_spec_europe_tax[datf$taxo_group=="hoverflies"]=913
datf$ntot_spec_europe=3051

###
datf %>% group_by(taxo_group) %>% summarise(length(unique(scientificName))/max(ntot_spec_europe_tax))


#proportion of species detected in each period
dat=datf %>% group_by(time_period,taxo_group) %>% summarise(prop_sp=length(unique(scientificName))/unique(ntot_spec_europe_tax))

p1a=ggplot(data=dat,aes(x=time_period,y=prop_sp,fill=taxo_group))+geom_bar(stat="identity",position=position_dodge())+theme_bw()+
scale_fill_manual(values=c("gold3","dodgerblue3"))+
  theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.border=element_blank(),
panel.grid.minor = element_blank(),panel.background = element_blank(),axis.text.x=element_text(angle=45,hjust=1),legend.position="none",plot.title=element_text(size=16,face="bold",hjust = 0))+ggtitle("a")+
  scale_y_continuous(labels=scales::percent)+
  ylab("Percentage of species observed")+xlab("Time period")



#nb records per group per year
dat2=datf %>% group_by(endYear,taxo_group) %>% count()

p1b=ggplot(data=dat2,aes(x=endYear,y=n,color=taxo_group))+geom_line(size=1.2)+
scale_color_manual(values=c("gold3","dodgerblue3"))+
theme_bw()+ theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.border=element_blank(),
panel.grid.minor = element_blank(),panel.background = element_blank(),plot.title=element_text(size=16,face="bold",hjust = 0),plot.subtitle=element_text(size=12),legend.position="none")+xlab("Years")+ylab("Number of records")+labs(color="region")+ggtitle("b")+
 scale_y_log10() 

top=plot_grid(p1a,p1b,ncol=1,align="v")

############# MAP:
# Define the UTM projection for a suitable UTM zone
utm_crs <- st_crs("+proj=utm +zone=32 +ellps=WGS84")
hex_grid=st_read(paste0(project_folder,"data/raw_data/grids_shapefiles/grid_",50,"KM.shp"),crs=utm_crs,quiet =TRUE)

sites=datf %>% group_by(gridID_50,region_50) %>% count()
sites$region_50[!(sites$region_50 %in% c("alpine","boreal","atlantic","continental","mediterranean"))]="other regions"

obj=merge(hex_grid,sites,by="gridID_50")

p1c=ggplot()+theme_bw()+
  geom_sf(data=obj,aes(fill=region_50.y,alpha=n),color=NA)+xlab("Longitude") + ylab("Latitude")+
  scale_color_manual(values=c(colo2,"grey"),na.value=NA)+
  scale_alpha_continuous(trans = "log10")+
  scale_fill_manual(values=c(colo2,"grey"),na.value=NA)+labs(fill="Bioregions",color="Bioregions",alpha="nb. records")+
  theme(plot.title=element_text(size=16,face="bold",hjust = 0),plot.subtitle=element_text(size=14),panel.border = element_blank(),axis.title=element_blank(),axis.ticks=element_blank(),axis.line=element_blank(),axis.text=element_blank(),legend.position="bottom",legend.box="vertical")+scale_x_continuous(n.breaks=3)+ coord_sf(ylim = c( 3210000,10008220),xlim = c(NA,4700000),clip = "on",expand = F)+ggtitle("c")


pdf(paste0(project_folder,"figures/Figure_1.pdf"),width=8,height=6)
grid.arrange(top,p1c,widths=c(1,1.5))
dev.off()
