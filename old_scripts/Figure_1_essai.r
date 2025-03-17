pkgs <- c("data.table", "dplyr","ggplot2","sf","viridis","cowplot","ggeffects","lme4","metafor","emmeans","gridExtra") 

inst <- pkgs %in% installed.packages()
if (any(inst)) install.packages(pkgs[!inst])
pkg_out <- lapply(pkgs, require, character.only = TRUE)

#defining working folder:
path_data="C:/Users/Duchenne/Documents/safeguard/data/"

# Loading data numbers
load(paste0(path_data,"list_filtering.RData"))

nb_records_initial=scales::comma(list_filtering[[1]])

nr_regions=scales::comma(list_filtering[[2]])

nb_records=scales::comma(list_filtering[[3]])

nr_month=scales::comma(list_filtering[[4]])

nr_sites=list_filtering[[5]]

nb_surveys=scales::comma(list_filtering[[6]])

nb_sp_common=list_filtering[[7]]

nb_sp=scales::comma(list_filtering[[8]])

nb_sp_tot=scales::comma(list_filtering[[9]])


# Define the UTM projection for a suitable UTM zone
utm_crs <- st_crs("+proj=utm +zone=32 +ellps=WGS84")

# Define bounding box coordinates for Spain to northern Finland in the chosen UTM zone
bbox <- st_bbox(c(xmin = -1200000, xmax = 3500000, ymin = 3500000, ymax = 8000000), crs = utm_crs)
bboxsf=st_as_sf(st_as_sfc(st_bbox(bbox)))

#shapefile of biogeogrpahic regions:
bioregions <- bioregions <- st_read ("D:/land use change/biogeographic_regions/BiogeoRegions2016.shp")
bioregions <- st_transform(bioregions, crs = utm_crs)
bioregions2=st_crop(bioregions,bbox)

colo=c("dodgerblue3","orange","lightcyan","chartreuse3","black","darkslategray2","olivedrab1","gold3",NA,"purple","salmon")

pl1=ggplot()+theme_bw()+
  geom_sf(data=bioregions2,aes(fill=short_name),color=NA)+xlab("Longitude") + ylab("Latitude")+
  geom_sf(data=bioregions,fill=NA)+
  scale_fill_manual(values=colo,na.value =NA)+labs(fill="Bioregions")+scale_x_continuous(n.breaks=3)+coord_sf(ylim = c( 3210000,10008220),xlim = c(NA,4800000),clip = "on",expand = F)#+geom_sf(data=bboxsf,color="black",fill=NA)



hex_grid=st_read(paste0(path_data,"grid_",50,"KM.shp"),crs=utm_crs,quiet =TRUE)

sites=fread(paste0(path_data,"sites_studied_tot.csv"))

obj=merge(hex_grid,sites,by="gridID_50")


colo2=c("#44AA99","#117733","#332288","#CC6677","#DDCC77")


pl2a=ggplot()+theme_bw()+
  geom_sf(data=bioregions,fill=NA)+
  geom_sf(data=obj,aes(fill=region_50.y,alpha=nb_records),color=NA)+xlab("Longitude") + ylab("Latitude")+
  scale_color_manual(values=colo2,na.value=NA)+
  scale_alpha_continuous(trans = "log10")+
  scale_fill_manual(values=colo2,na.value=NA)+labs(fill="Bioregions",color="Bioregions",alpha="nb. records")+
  theme(plot.title=element_text(size=14,face="bold",hjust = 0),plot.subtitle=element_text(size=14),panel.border = element_blank(),axis.title=element_blank(),axis.ticks=element_blank(),axis.line=element_blank(),axis.text=element_blank(),legend.position="bottom",legend.box="vertical")+scale_x_continuous(n.breaks=3)+ coord_sf(ylim = c( 3210000,10008220),xlim = c(NA,4800000),clip = "on",expand = F)+ggtitle("a")
  
blank=ggplot() + theme(plot.title=element_text(size=14,face="bold",hjust = 0),plot.subtitle=element_text(size=14),panel.border = element_blank(),axis.title=element_blank(),axis.ticks=element_blank(),axis.line=element_blank(),panel.background=element_blank())+ggtitle("b")



pdf("map.pdf",width=7,height=5)
grid.arrange(pl2a,blank,blank,layout_matrix=rbind(c(1,2),c(1,3)),widths=c(2,1))
dev.off();

