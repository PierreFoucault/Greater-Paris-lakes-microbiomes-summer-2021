
#### Import packages ####

library(ggplot2)
library(ggpubr)
library(tidyverse)
library(dplyr)
library(readr)
library(stringr)
library(tidytext)
library(AggregateR)
library(ggh4x)
library(patchwork)
library(ggdendro)
library(ggbreak)
library(gg.gap)
library(MASS)

library(grImport2)
library(rsvg)
library(cowplot)
library(gridSVG)
library(magick)
library(factoextra)

#If not already downloaded
#
BiocManager::install("ggfortify")

####. ####

#### Import data ####

Data_param<- read_delim("Metadata_COM2LIFE_eng.csv",";",
                        escape_double = FALSE, trim_ws = TRUE)


Data_param<- Data_param %>% 
  replace(is.na(.), 0) %>%
  dplyr::group_by(Month_Site,.drop = FALSE) %>%
  dplyr::summarise_at(vars(7:(ncol(Data_param)-1)),
                      list(median=median, sd=sd)) %>%
  cbind(Site_letter=(Data_param$Site_letter[match(.$Month_Site,Data_param$Month_Site)]),
        Month_letter=(Data_param$Month_letter[match(.$Month_Site,Data_param$Month_Site)]),
        .) %>%
  mutate(across(Site_letter, factor,
                levels=c("CHA","VER","GDP","BOI","CTL","CRJ1","CRJ2","VAI","JAB")))
#View(Data_param)

# transform SVG to cairo format
rsvg_svg("/Users/pierre/Desktop/PhD/COM2LIFE/QGIS/map_COM2LIFE.svg", "map_COM2LIFE-cairo.svg")
# pre-read the cairo svg into an object
SVG_map_cairo <- readPicture("map_COM2LIFE-cairo.svg")


library(magick)

gridsvg("/Users/pierre/Desktop/PhD/COM2LIFE/QGIS/map_COM2LIFE.svg")
grid.rect(gp=gpar(col=NA, fill="grey80"))
SVG_map_gg <- pictureGrob(SVG_map_cairo, ext="gridSVG", delayContent=FALSE)
plot_grid(SVG_map_gg)
dev.off()

SVG_map_gg<-image_read("/Users/pierre/Desktop/PhD/COM2LIFE/QGIS/map_COM2LIFE.svg") |> image_ggplot()

####. ####

#### Palettes & strip ####

#set color order, indianred2 for full month, gray38 for usual.

palette_month<-c("gray38","indianred2","gray38","gray38","indianred2","gray38","indianred2")
palette_site<-c("#84b231","#84b231","#84b231",
                "#D8C5B2","#D8C5B2","#D8C5B2",
                "#377EB8","#377EB8","#377EB8")

# To color the horizontal strips (=site name) according to the 3 levels of eutrophication 

strip_color_lake<- strip_themed(
  background_x = elem_list_rect(fill = palette_site,
                                color = palette_site),
  text_x = elem_list_text(colour = "white",
                          face = "bold",size=10))

####. ####

#### Chla ####

pdf(file="Fig_Chla.pdf",
    width = 8.5, # The width of the plot in inches
    height = 6) # The height of the plot in inches

Fig_Chla <- ggplot(Data_param,aes(x=Month_letter,y=Chla_median))+theme_bw()+
  geom_bar(stat="identity",width=0.3,fill="darkgreen",alpha=0.8)+
  geom_linerange(aes(x=Month_letter,
                      ymin=if_else(Chla_median-Chla_sd>0,
                                   Chla_median-Chla_sd,0),
                      ymax=Chla_median+Chla_sd), colour="black", alpha=0.8, size=0.65)+
  theme(plot.title = element_text(hjust=0.5,size=12,face = "bold"),
        axis.title = element_blank(),
        axis.text.x=element_text(size=7.5,
                                 colour=palette_month,face="bold"),
        axis.ticks.x = element_blank())+
  labs(title="Chla (µg/L)")+ylim(c(0,NA))+
  facet_wrap2(~Site_letter, scales ="free_y", strip = strip_color_lake)



Fig_Chla
dev.off()

View(Data_param)
Fig_Chla_seuil<-
  ggplot(Data_param, aes(x=as.factor(Month_letter),y=Chla_median,group = 1 ))+theme_bw()+
  geom_hline(yintercept=40,color="green")+
  geom_hline(yintercept=10,color="blue")+
  geom_point(color="black",alpha=0.8)+
  geom_line(color="darkgreen",alpha=0.8)+
  geom_linerange(aes(x=Month_letter,
                     ymin=if_else(Chla_median-Chla_sd>0,
                                  Chla_median-Chla_sd,0),
                     ymax=Chla_median+Chla_sd), color="black", alpha=0.8, size=0.65)+
  theme(plot.title = element_blank(),
        panel.grid = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=14),
        axis.text.x=element_text(size=7.5,
                                 colour="black",face="bold"),
        axis.ticks.x = element_blank())+
  labs(y="Chl a (µg/L) log10 scaled\n")+
  scale_y_continuous(trans="log10",breaks = c(0,10,40,100,200))+
  facet_wrap2(~Site_letter, scales ="free_x", strip = strip_color_lake)

Fig_Chla_seuil
 

####.####

#### pH ####

pdf(file="Fig_pH.pdf",
    width = 8.5, # The width of the plot in inches
    height = 6) # The height of the plot in inches

Fig_pH <- ggplot(Data_param,aes(x=Month_letter,y=pH_median))+theme_bw()+
  geom_bar(stat="identity",width=0.3,fill="indianred2",alpha=0.8)+
  geom_linerange( aes(x=Month_letter,
                      ymin=if_else(pH_median-pH_sd>0,
                                   pH_median-pH_sd,0),
                      ymax=pH_median+pH_sd), colour="black", alpha=0.8, size=0.65)+
  theme(plot.title = element_text(hjust=0.5,size=12,face = "bold"),
        axis.title = element_blank(),
        axis.text.x=element_text(size=7.5,
                                 colour=palette_month,face="bold"),
        axis.ticks.x = element_blank())+
  labs(title="pH")+ylim(c(0,(0.2+max(Data_param$pH_median))))+
  facet_wrap2(~Site_letter, scales ="free", strip = strip_color_lake)

Fig_pH$data[,c(1,2,6)] %>% subset(., Month_letter %in% c("A","B","C","D")) %>%
  dplyr::summarise_at(vars("pH_median"),list(median=median, sd=sd))

Fig_Salinity$data[,c(1,2,4)] %>% subset(., Month_letter %in% c("A","B","C","D")) %>%
  dplyr::summarise_at(vars("Salinity_median"),list(median=median, sd=sd))

View(Fig_pH$data[,c(1,2,6)] %>% subset(., Month_letter %in% c("A","B","C","D")))
Fig_pH
dev.off()

param_data[,c(1,2,7)] %>%
  #subset(., month %in% c("A","B","C","D")) %>%
  #group_by(month) %>%
  dplyr::summarise_at(vars("pH_median"),list(mean=mean, sd=sd))

####.####

#### Salinity ####

pdf(file="Fig_Salinity.pdf",
    width = 8.5, # The width of the plot in inches
    height = 6) # The height of the plot in inches

Fig_Salinity <- ggplot(Data_param,aes(x=Month_letter,y=Salinity_median))+theme_bw()+
  geom_bar(stat="identity",width=0.3,fill="lightblue",alpha=0.8)+
  geom_linerange( aes(x=Month_letter,
                      ymin=if_else(Salinity_median-Salinity_sd>0,
                                   Salinity_median-Salinity_sd,0),
                      ymax=Salinity_median+Salinity_sd), colour="black", alpha=0.8, size=0.65)+
  theme(plot.title = element_text(hjust=0.5,size=12,face = "bold"),
        axis.title = element_blank(),
        axis.text.x=element_text(size=7.5,
                                 colour=palette_month,face="bold"),
        axis.ticks.x = element_blank())+
  labs(title="Salinity")+ylim(c(0,NA))+
  facet_wrap2(~Site_letter, scales ="free", strip = strip_color_lake)

Fig_Salinity
dev.off()

####.####

#### P ####

pdf(file="Fig_Phosphate.pdf",
    width = 8.5, # The width of the plot in inches
    height = 6) # The height of the plot in inches


Fig_Phosphate <- ggplot(Data_param,aes(x=Month_letter,y=Phosphate_median))+theme_bw()+
  geom_bar(stat="identity",width=0.3,fill="#66C2A5",alpha=1)+
  geom_linerange( aes(x=Month_letter,
                      ymin=if_else(Phosphate_median-Phosphate_sd>0,
                                   Phosphate_median-Phosphate_sd,0),
                      ymax=Phosphate_median+Phosphate_sd), colour="black", alpha=0.8, size=0.65)+
  theme(plot.title = element_text(hjust=0.5,size=12,face = "bold"),
        axis.title = element_blank(),
        axis.text.x=element_text(size=7.5,
                                 colour=palette_month,face="bold"),
        axis.ticks.x = element_blank())+
  labs(title="Phosphate"+ylim(c(0,NA)))+
  facet_wrap2(~Site_letter, scales ="free", strip = strip_color_lake)

Fig_Phosphate
dev.off()

####.####

#### MIO ####

N_P_summer<-read_delim("N_P_summer.csv",";", escape_double = FALSE,
                       trim_ws = TRUE,show_col_types = FALSE) %>%
  dplyr::mutate(
    lake=factor(lake,
                levels=c("CHA","VER","GDP","BOI","CTL","CRJ1","CRJ2","VAI","JAB")),
    lake_type=if_else(lake %in% c("CHA","VER","GDP"),"Hypereutrophic",
                      if_else(lake %in% c("BOI","CTL"),"Eutrophic",
                              if_else(lake %in% c("CRJ1","CRJ2"),"Mesotrophic",
                                      "Oligotrophic"))),
    lake_type=factor(lake_type,levels=c("Hypereutrophic","Eutrophic","Mesotrophic","Oligotrophic")),
    lake_month=paste0(lake,"_",month))

head(N_P_summer)

NH4_summer<-ggplot(N_P_summer,aes(x=month,y=NH4,fill=lake))+
  geom_point(show.legend = F)+geom_boxplot(show.legend = F)+theme_bw()+
  theme(axis.title = element_blank(),panel.grid = element_blank(),
        plot.title = element_text(hjust=0,size=15,face="bold"))+
  scale_fill_manual(values=palette_lake)+
  facet_wrap2(~ lake, scales ="free_x",ncol = 1,strip.position = "left",strip = strip_color_lake_y)+
  labs(title = expression(paste(NH[4]^"+")))

NO_summer<-ggplot(N_P_summer,aes(x=month,y=NO3_NO2,fill=lake))+
  geom_point(show.legend = F)+geom_boxplot(show.legend = F)+theme_bw()+
  theme(axis.title = element_blank(),panel.grid = element_blank(),
        plot.title = element_text(hjust=0,size=15,face="bold"))+
  scale_fill_manual(values=palette_lake)+
  facet_wrap2(~ lake, scales ="free_x",ncol = 1,strip.position = "left",strip = strip_color_lake_y)+
  labs(title = expression(paste(NO[3]^"-"," & ",NO[2]^"-")))

P_summer<-ggplot(N_P_summer,aes(x=month,y=PO4,fill=lake))+
  geom_point(show.legend = F)+geom_boxplot(show.legend = F)+theme_bw()+
  theme(axis.title = element_blank(),panel.grid = element_blank(),
        plot.title = element_text(hjust=0,size=15,face="bold"))+
  scale_fill_manual(values=palette_lake)+
  facet_wrap2(~ lake, scales ="free_x",ncol = 1,strip.position = "left",strip = strip_color_lake_y)+
  labs(title = expression(paste(PO[4]^"3-")))

NH4_summer+NO_summer+P_summer

####.####

#### POC ####

TPCN_df<-read_delim("POC.csv",";",
                 escape_double = FALSE,trim_ws = TRUE,show_col_types = FALSE) %>%
  dplyr::mutate(.,lake=factor(lake,
                levels=c("CHA","GDP","VER","CTL","BOI","CRJ1","CRJ2","JAB","VAI")))

Fig_TPC <- ggplot(TPCN_df,aes(x=month,y=TPC,color=lake,group=lake))+theme_bw()+
  geom_line(show.legend = F)+
  theme(axis.title.y = element_text(hjust=0.5,size=15,face = "bold"),
        panel.grid = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10,face="bold"),
        axis.ticks = element_blank())+
  scale_color_manual(values =palette_site )+ 
  scale_y_continuous(trans="log10",limits = c(NA,max(POC_df$POC_C)),breaks = c(0.2,0.5,1,2.5,5))+
  labs(y = "Total Particular Carbon in ug/mL (log10)\n")+
  facet_wrap2(~lake, scales ="fixed", strip = strip_color_lake)

Fig_TPC

Fig_TPN <- ggplot(TPCN_df,aes(x=month,y=TPN,color=lake,group=lake))+theme_bw()+
  geom_line(show.legend = F)+
  theme(axis.title.y = element_text(hjust=0.5,size=15,face = "bold"),
        panel.grid = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10,face="bold"),
        axis.ticks = element_blank())+
  scale_color_manual(values = palette_site)+ 
  #scale_y_continuous(trans="log10",limits = c(NA,max(POC_df$POC_N)),breaks = c(0.2,0.5,1,2.5,5))+
  labs(y = "Total Particular Nitrogen in ug/mL\n")+
  facet_wrap2(~lake, scales ="fixed", strip = strip_color_lake)

Fig_TPN

Fig_TPC_N <- ggplot(TPCN_df,aes(x=month,y=`TPC/TPN`,color=lake,group=lake))+theme_bw()+
  geom_line(show.legend = F)+
  theme(axis.title.y = element_text(hjust=0.5,size=15,face = "bold"),
        panel.grid = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10,face="bold"),
        axis.ticks = element_blank())+
  scale_color_manual(values = palette_site)+ 
  #scale_y_continuous(trans="log10",limits = c(NA,max(POC_df$POC_N)),breaks = c(0.2,0.5,1,2.5,5))+
  labs(y = "Total Particular C/N in ug/mL\n")+
  facet_wrap2(~lake, scales ="fixed", strip = strip_color_lake)

Fig_TPC_N

####.####

#### geo dist ####

# a phyloseq object with LON LAT as samdata
A_D_rarW

# import into a new dataframe with sample names
Havers_data<- data.frame(sample=A_D_rarW@sam_data$sample,LON=A_D_rarW@sam_data$LON,
                       LAT=A_D_rarW@sam_data$LAT) %>%
  column_to_rownames(.,var="sample")

# Create a distance matrix and change into a dataframe
Havers_matrix<-distm(Havers_data,fun = distHaversine) %>% as.data.frame()
?distHaversine
colnames(Havers_matrix)=A_D_rarW@sam_data$sample
Havers_matrix=cbind(sample=A_D_rarW@sam_data$sample,Havers_matrix)

# Transform the data in a long table
Havers_df<-Havers_matrix %>% pivot_longer(-sample,names_to = "sample2",values_to = "Haver_dist") %>%
  filter(as.character(sample) != as.character(sample2))

View(Havers_df)
View(Havers_data)

####.####

#### Land Cover ####

CLC_codex<-read_delim("CLC_codex.csv",delim=";")

CLC.df<-read_delim("land_cover_1km.csv",delim=";") %>%
  dplyr::mutate(.,land_type=CLC_codex$land_type[match(.$CLC_code,CLC_codex$CLC_code)]) %>%
  group_by(lake,land_type) %>% dplyr::summarise(.,sum_aera=sum(Area)) %>%
  group_by(lake) %>% dplyr::mutate(.,total_area=sum(sum_aera)) %>% ungroup() %>%
  dplyr::mutate(.,arear_percent=(sum_aera/total_area)*100)

Fig_CLC <- CLC.df %>%
  ggplot(.,aes(lake,arear_percent,fill=land_type))+
  geom_col()+theme_bw()+
  scale_fill_manual(values = palette_Charlotte)

# %>%
#   group_by(lake,land_type) %>% dplyr::summarise(land_percent=(sum(Area)/sum_aera)*100)
#   


####.####

#### Lake status ####

####__ parameters ####
# get data
all_A_E<-subset_samples(all_rarW,month %in% c("A","B","C","D"))

param_data <- all_A_E@sam_data %>% as.data.frame() %>%
  dplyr::group_by(lake_month,.drop = FALSE) %>%
  dplyr::summarise_at(vars(14:(ncol(all_A_E@sam_data)-1)),
                      list(median=median)) %>%
  cbind(lake=(all_A_E@sam_data$lake[match(.$lake_month,all_A_E@sam_data$lake_month)]),
        month=(all_A_E@sam_data$month[match(.$lake_month,all_A_E@sam_data$lake_month)]),
        .) %>% .[,c(-8)]

head(param_data)
colnames(param_data)

param_table<-param_data[,c(3,4:ncol(param_data))] %>% as.data.frame() %>% column_to_rownames(var = "lake_month" )
#View(param_table)
write.csv(param_table,"param_table_AD.csv")

ggplot(param_table,aes(y=Oxygen_median,x=Chla_median))+geom_point()

N_P_summer_median <-N_P_summer %>%
  dplyr::group_by(lake_month,.drop = FALSE) %>%
  dplyr::summarise_at(vars(6:8),list(median=median))

param_table_AD<-read_delim("param_table_AD.csv",";",show_col_types = F) %>% column_to_rownames(var = "...1") %>%
  dplyr::mutate(.,
                lake=factor(lake,levels = c("CHA","VER","GDP","BOI","CTL","CRJ1","CRJ2","VAI","JAB")),
                month_full=factor(month_full,levels = c("June","July","August","September"))) %>%
  dplyr::mutate(.,
    NH4=(N_P_summer_median$median_NH4_median[match(rownames(.),N_P_summer_median$lake_month)]),
    PO4=(N_P_summer_median$median_PO4_median[match(rownames(.),N_P_summer_median$lake_month)]),
    NO3_N02=(N_P_summer_median$median_NO3_NO2_median[match(rownames(.),N_P_summer_median$lake_month)]))
        
colnames(param_table_AD)
View(param_table_AD)
write.csv(param_table_AD[,-6],"param_table_2021.csv")

####__ PCA ####

param_table_AD %>% .[,c(12,13,19,20,21)] 
MDS_param.df<- param_table_AD %>% .[,c(12,13,17,19,20,21)] %>% group_by(lake) %>%
  dplyr::summarise_at(vars(1,2,3,5),list(mean=mean, sd=sd))

param_table_AD %>% .[,c(5,12,13,17,19,20,21)] %>% dplyr::summarise_at(vars(1,2,3,5,6,7),list(min=min, max=max))
View(param_table_AD)
tmp<-read_delim("param_table_AD.csv",delim = ";") %>% rename("lake_month"="...1") %>% column_to_rownames(.,var="lake_month") %>%
  dplyr::mutate(.,lake=param_table_AD$lake[match(rownames(.),rownames(param_table_AD))],.before=1)

MDS_param.df<- tmp %>%.[,c(1,2,4,5,7,8,10,12,13,19,20,21,22,23)] %>%
  scale(., center = T, scale = T)


MDS_param.df<- param_table_AD %>% .[,c(12,13,19,20,21)] %>%
  scale(., center = T, scale = T)

colnames(MDS_param.df)
MDS_param.eucli<-MDS_param.df %>% vegdist(., method="euclidean")
MDS_param.eucli.mds<- cmdscale(MDS_param.eucli,eig=TRUE, k=2)

pca_res <- prcomp(tmp %>%.[,c(2,4,5,7,8,10,12,13,19,20,21,22,23)], scale. = T,center = T)
var <- get_pca_var(pca_res)
var_dim1<-var$cos2 %>% .[,c(1,2)]
melted_corr <- melt(var_dim1)

corr_Fig<- melted_corr %>%  subset(.,Var2=="Dim.1") %>%
  ggplot(.,aes(x = Var2, y = Var1, fill = value, size=value))+
  geom_point(color='black',shape=21) + theme_bw()+
  geom_text(aes(x = Var2, y = Var1, label = round(value, 2)), color = "white",size = 3.8)+
  labs(fill = expression(Cos^"2"),size="",x="PC1 - [44%]")+
  theme(panel.grid=element_blank(),
        axis.title.x = element_text(face="bold",size = 15,hjust=0.5),
        axis.text.y = element_text(face="bold",size = 12),
        axis.text.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks = element_blank(),
        legend.title = element_text(face="bold",size = 10))+
  scale_fill_gradient2(midpoint = 0.5, mid ="grey70", low = "blue", high = "red",
                       limits = c(0, +1)) +
  scale_size_continuous(range = c(9,18),breaks = c(0.2,0.4,0.8),labels = c("0.2","0.4","0.6"))
corr_Fig

View(pca)
rownames(param_MDS_eucli.df)
View(MDS_param.eucli.mds)
param_MDS_eucli.df <- MDS_param.eucli.mds$points %>% as.data.frame(.) %>%
  rename(.,"dim1"="V1","dim2"="V2") %>%
  dplyr::mutate(.,
                month_full=param_table_AD$month_full[match(rownames(.),rownames(param_table_AD))],
                lake=param_table_AD$lake[match(rownames(.),rownames(param_table_AD))],
                lake=
                  if_else(lake=="CHA","Champs",
                          if_else(lake=="VER","Verneuil",
                                  if_else(lake=="GDP","Grande-Paroisse",
                                          if_else(lake=="BOI","Bois",
                                                  if_else(lake=="CTL","Créteil",
                                                          if_else(lake=="CRJ1","Cergy-small",
                                                                  if_else(lake=="CRJ2","Cergy-large",
                                                                          if_else(lake=="VAI","Vaires","Jablines")))))))),
                lake=factor(lake,levels=c("Verneuil","Champs",
                                          "Grande-Paroisse","Bois","Créteil","Cergy-small","Cergy-large",
                                          "Vaires","Jablines")))
                              
param_env <- envfit(MDS_param.eucli.mds$points ~ NH4+TPN+TPC+NO3_N02+PO4, data = as.data.frame(MDS_param.df), perm = 999)
param_env.df <- as.data.frame(scores(param_env, display = "vectors")) 
param_env.df <- cbind(param_env.df, param = rownames(param_env.df)) 

biplot()
# GGplot !
Fig_mds_param <- param_MDS_eucli.df %>%
  ggplot(.,aes(dim1,dim2,color = lake)) +
  geom_convexhull(aes(color = lake,fill=factor(lake),group = lake),alpha=0.4,show.legend = F,size=0.4)+
  geom_segment(data = param_env.df,
               mapping=aes(x = 0, xend = Dim1*3, y = 0, yend = Dim2*3),color = "black",inherit.aes = F,arrow = arrow(length = unit(0.25,"cm"),type = "closed")) +
  geom_text(data = param_env.df, mapping=aes(x = Dim1*3.12, y = Dim2*3.12, label = param),size = 3.5,fill="white",fontface = "bold",inherit.aes = F)+
  geom_point(size=2,shape=16,show.legend = F) + theme_bw() +
  theme(plot.title=element_blank(),
        axis.title = element_text(size=12,face = "bold"),
        axis.text = element_text(size=10),
        legend.title = element_text(size=15,face = "bold"),
        legend.text = element_text(size=12),
        plot.caption = element_text(size=12,face = "italic",color="black"),
        panel.grid = element_blank(),
        axis.ticks = element_blank())+
  scale_color_manual(values=palette_lake_3T)+
  scale_fill_manual(values=palette_lake_3T)+
  #scale_x_reverse()+
  scale_y_reverse()+
  #scale_shape_manual(values = Summer_2021)+
  guides(color = guide_legend(override.aes=list(size=6)),
         fill = guide_legend(override.aes=list(size=6)))+
  labs(color="Lake",shape="Month",
       x=paste0("PC1 [",round(MDS_param.eucli.mds$eig[1]*100/sum(MDS_param.eucli.mds$eig),1),"%]"),
       y=paste0("PC2 [",round(MDS_param.eucli.mds$eig[2]*100/sum(MDS_param.eucli.mds$eig),1),"%]"))
       #title=expression(paste("TPN, TPC, ",NH[4]^"-",", ",NO[3]^"-"," & ",NO[2]^"-",", ",PO[4]^"3-")),
       #title="Enviromental parameters and lake morphometry PCA")
Fig_mds_param

####__ PCA 1 ####

param_centroids_MDS_eucli.df <- param_MDS_eucli.df %>% dplyr::group_by(lake) %>% dplyr::summarise_at(vars(1:2),median)

Fig_lake_dim1<-ggplot(param_MDS_eucli.df,aes(x=0,y=dim1,color=lake,label = lake, fill = lake))+
  geom_point(size=1,show.legend = F)+theme_bw()+
  geom_point(param_centroids_MDS_eucli.df,inherit.aes = F,size=2.75,mapping=aes(x=0,y=dim1,color=lake),show.legend = F)+
  geom_label(param_centroids_MDS_eucli.df,color = "white",size = 4,fontface = "bold",show.legend = F,inherit.aes = F,
             mapping=aes(label = lake, fill = lake,
                         x=if_else(lake %in% c("Bois","Vaires","Verneuil","Créteil"),-0.05,
                                   if_else(lake %in% c("Cergy-large"),0.025,
                                           if_else(lake %in% c("Cergy-small"),-0.025,
                                           0.05))),
                         y=if_else(lake %in% c("Grande-Paroisse"),dim1-0.25,
                           if_else(lake %in% c("Bois"),dim1-0.2,
                           if_else(lake %in% c("Cergy-large","Jablines","Vaires"),dim1-0.2,
                           if_else(lake %in% c("Cergy-small"),dim1-0.12,dim1))))))+
  theme(panel.grid = element_blank(),
        axis.ticks = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_text(size=14),
        legend.title = element_text(size=15,face = "bold"),
        legend.text = element_text(size=12),
        plot.caption = element_text(size=12,face = "italic",color="black"),
        plot.title = element_text(hjust=0.5,size=15,face = "bold"),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size=20,face = "bold"))+
  scale_color_manual(values=palette_lake_3T)+
  scale_fill_manual(values=palette_lake_3T)+
  scale_x_continuous(limits = c(-0.2,0.2))+
  scale_y_continuous(breaks = c(-2,0,2,4))+
  #scale_y_reverse()+
  scale_x_reverse()+
  labs(y="Nutrients, TPC and TPN PC1 [44%]")+
       #caption=paste0("Abiotic parameters PCA 1 [",round(MDS_param.eucli.mds$eig[1]*100/sum(MDS_param.eucli.mds$eig),1),"%]"))+
  coord_flip()


####__ Chl-a ####
Chla_mean_MDS_eucli.df <- param_table_AD %>%.[,c(5,17)] %>% dplyr::group_by(lake) %>%
  dplyr::summarise_at(vars(1),list(mean=mean)) #

Chla_mean.df <- param_table_AD %>%.[,c(5,17)] %>% dplyr::group_by(lake) %>%
  dplyr::summarise_at(vars(1),list(mean=mean)) #

write.csv(param_table_AD,"param_table_AD.csv")

Fig_lake_median_Chla_hor<-param_table_AD %>% dplyr::group_by(lake) %>%
  dplyr::mutate(.,Chla_median_lake=mean(Chla_median)) %>%
  ggplot(.,aes(x=0,y=Chla_median,color=lake,label = lake, fill = lake))+
  geom_hline(yintercept = c(2.6), color = "grey", linetype = 2,size=0.6)+
  geom_hline(yintercept = c(7.3), color = "grey", linetype = 2,size=0.6)+
  geom_hline(yintercept = c(56), color = "grey", linetype = 2,size=0.6)+
  #geom_label(aes(x=-0.14,y=2),label="Oligotrophic",size=6,color="white",fill="#8CB8E2",fontface = "bold")+
  geom_label(aes(x=-0.14,y=4.1),label="Mesotrophic",size=6,color="white",fill="#E19B8E",fontface = "bold")+
  geom_label(aes(x=-0.14,y=19),label="Eutrophic",size=6,color="white",fill="#b19774",fontface = "bold")+
  geom_label(aes(x=-0.14,y=100),label="Hypereutrophic",size=6,color="white",fill="#71992A",fontface = "bold")+
  geom_point(size=1,show.legend = F)+theme_bw()+
  geom_point(size=2.75,mapping=aes(x=0,y=Chla_median_lake,color=lake),show.legend = F)+
  geom_label(color = "white",size = 5,fontface = "bold",show.legend = F,
             mapping=aes(label = lake, fill = lake,
                         x=if_else(lake %in% c("BOI","CRJ1","VAI","CHA"),-0.05,0.05),
                         y=Chla_median_lake))+
  theme(panel.grid = element_blank(),
        axis.ticks = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_text(size=14),
        legend.title = element_text(size=15,face = "bold"),
        legend.text = element_text(size=12),
        plot.title = element_text(hjust=0.5,size=15,face = "bold"),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size=20,face = "bold"))+
  scale_color_manual(values=palette_lake)+
  scale_fill_manual(values=palette_lake)+
  scale_x_continuous(limits = c(-0.2,0.2))+
  scale_y_continuous(trans="log10",breaks = c(0,2.6,7.3,56,200))+
  labs(y=expression(paste(bold("Chl"),bold(italic("a")))))+coord_flip()
Fig_lake_median_Chla_hor

Fig_lake_median_Chla_vert<-param_table_AD %>% dplyr::group_by(lake) %>%
  dplyr::mutate(.,Chla_median_lake=mean(Chla_median)) %>%
  ggplot(.,aes(x=0,y=Chla_median,color=lake,label = lake, fill = lake))+
  geom_hline(yintercept = c(2.6), color = "#8CB8E2", linetype = 2,size=0.6)+
  geom_hline(yintercept = c(7.3), color = "#b19774", linetype = 2,size=0.6)+
  geom_hline(yintercept = c(56), color = "#71992A", linetype = 2,size=0.6)+
  #geom_label(aes(x=-0.14,y=2),label="Oligotrophic",size=6,color="white",fill="#8CB8E2",fontface = "bold")+
  geom_label(aes(x=-0.14,y=2.6),label="Mesotrophic",size=5,color="white",fill="#8CB8E2",fontface = "bold")+
  geom_label(aes(x=-0.14,y=7.3),label="Eutrophic",size=5,color="white",fill="#b19774",fontface = "bold")+
  geom_label(aes(x=-0.14,y=56),label="Hypereutrophic",size=5,color="white",fill="#71992A",fontface = "bold")+
  geom_point(size=1,show.legend = F)+theme_bw()+
  geom_point(size=2.75,mapping=aes(x=0,y=Chla_median_lake,color=lake),show.legend = F)+
  geom_label(color = "white",size = 5,fontface = "bold",show.legend = F,
             mapping=aes(label = lake, fill = lake,
                         x=if_else(lake %in% c("BOI","CRJ1","VAI","CHA"),-0.05,0.05),
                         y=Chla_median_lake))+
  theme(panel.grid = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size=14),
        legend.title = element_text(size=15,face = "bold"),
        legend.text = element_text(size=12),
        plot.title = element_text(hjust=0.5,size=15,face = "bold"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=20,face = "bold"))+
  scale_color_manual(values=palette_lake)+
  scale_fill_manual(values=palette_lake)+
  scale_x_continuous(limits = c(-0.2,0.2))+
  scale_y_continuous(trans="log10",breaks = c(0,2.6,7.3,56,200))+
  labs(y=expression(paste(bold("Chl"),bold(italic("a")))))
Fig_lake_median_Chla_vert

####__ PO42- ####
View(param_table_AD)
P_mean_MDS_eucli.df <- param_table_AD %>%.[,c(6,17)] %>% dplyr::group_by(lake) %>% dplyr::summarise_at(vars(1),median)

Fig_lake_P_median<- ggplot(param_table_AD,aes(x=0,y=Phosphate_median,color=lake,label = lake, fill = lake))+
  geom_point(size=1,show.legend = F)+theme_bw()+
  geom_point(P_mean_MDS_eucli.df,inherit.aes = F,size=2.75,mapping=aes(x=0,y=Phosphate_median,color=lake),show.legend = F)+
  geom_label(P_mean_MDS_eucli.df,color = "white",size = 4,fontface = "bold",show.legend = F,inherit.aes = F,
             mapping=aes(label = lake, fill = lake,
                         x=if_else(lake %in% c("BOI","CRJ2","CHA","JAB"),-0.075,0.075),
                         y=if_else(lake %in% c("JAB"),Phosphate_median-0.02,Phosphate_median)))+
  theme(panel.grid = element_blank(),
        axis.ticks = element_blank(),
        axis.text.y = element_blank(),
        legend.title = element_text(size=15,face = "bold"),
        legend.text = element_text(size=12),
        plot.title = element_text(hjust=0.5,size=15,face = "bold"),
        axis.title = element_blank())+
  scale_color_manual(values=palette_lake)+
  scale_fill_manual(values=palette_lake)+
  scale_x_continuous(limits = c(-0.2,0.2))+
  scale_y_continuous(trans="log10")+
  labs(title="\nPO42-")+coord_flip()
Fig_lake_P_median

####__ TPC ####
TPC_mean_MDS_eucli.df <- param_table_AD %>%.[,c(13,17)] %>% dplyr::group_by(lake) %>% dplyr::summarise_at(vars(1),median)

Fig_lake_TPC_median<- ggplot(param_table_AD,aes(x=0,y=TPC,color=lake,label = lake, fill = lake))+
  geom_point(size=1,show.legend = F)+theme_bw()+
  geom_point(TPC_mean_MDS_eucli.df,inherit.aes = F,size=2.75,mapping=aes(x=0,y=TPC,color=lake),show.legend = F)+
  geom_label(TPC_mean_MDS_eucli.df,color = "white",size = 4,fontface = "bold",show.legend = F,inherit.aes = F,
             mapping=aes(label = lake, fill = lake,
                         x=if_else(lake %in% c("BOI","CRJ2","CHA","JAB"),-0.075,0.075),
                         y=TPC))+
  theme(panel.grid = element_blank(),
        axis.ticks = element_blank(),
        axis.text.y = element_blank(),
        legend.title = element_text(size=15,face = "bold"),
        legend.text = element_text(size=12),
        plot.title = element_text(hjust=0.5,size=15,face = "bold"),
        axis.title = element_blank())+
  scale_color_manual(values=palette_lake)+
  scale_fill_manual(values=palette_lake)+
  scale_x_continuous(limits = c(-0.2,0.2))+
  scale_y_continuous(trans="log10")+
  labs(title="\nTPC")+coord_flip()
Fig_lake_TPC_median

####__ TPN ####
TPN_mean_MDS_eucli.df <- param_table_AD %>%.[,c(12,17)] %>% dplyr::group_by(lake) %>% dplyr::summarise_at(vars(1),median)

Fig_lake_TPN_median<- ggplot(param_table_AD,aes(x=0,y=TPN,color=lake,label = lake, fill = lake))+
  geom_point(size=1,show.legend = F)+theme_bw()+
  geom_point(TPN_mean_MDS_eucli.df,inherit.aes = F,size=2.75,mapping=aes(x=0,y=TPN,color=lake),show.legend = F)+
  geom_label(TPN_mean_MDS_eucli.df,color = "white",size = 4,fontface = "bold",show.legend = F,inherit.aes = F,
             mapping=aes(label = lake, fill = lake,
                         x=if_else(lake %in% c("BOI","CRJ2","CHA","JAB"),-0.075,0.075),
                         y=TPN))+
  theme(panel.grid = element_blank(),
        axis.ticks = element_blank(),
        axis.text.y = element_blank(),
        legend.title = element_text(size=15,face = "bold"),
        legend.text = element_text(size=12),
        plot.title = element_text(hjust=0.5,size=15,face = "bold"),
        axis.title = element_blank())+
  scale_color_manual(values=palette_lake)+
  scale_fill_manual(values=palette_lake)+
  scale_x_continuous(limits = c(-0.2,0.2))+
  scale_y_continuous(trans="log10")+
  labs(title="\nTPN")+coord_flip()
Fig_lake_TPN_median

####__ Chla vs MDS-dim1 ####
Chla_dim1_centroids.df<-data.frame(sample=rownames(param_table_AD),Chla=param_table_AD$Chla_median,lake=param_table_AD$lake,month_full=param_table_AD$month_full) %>%
  remove_rownames() %>% column_to_rownames(var = "sample") %>%
  dplyr::mutate(.,MDS_dim1=param_MDS_eucli.df$dim1[match(rownames(.),rownames(param_MDS_eucli.df))],.after=1) %>%
  group_by(lake) %>% 
  dplyr::mutate(.,MDS_dim1_median=median(MDS_dim1),Chla_mean=mean(Chla),.before=1) %>% group_by(lake) %>% 
  dplyr::summarise_at(vars(1:2),median)
View(Chla_dim1_centroids.df)

spearman_gradient<-cor(Fig_Chla_dim1$data$MDS_dim1,Fig_Chla_dim1$data$Chla,method = "spearman")
lm_gradient_model <- lm(MDS_dim1 ~ Chla,data=Fig_Chla_dim1$data)
var_graident=paste("Spearman",round(spearman_gradient,2),", ","R2",round(summary(lm_gradient_model)[["adj.r.squared"]],2),", p<0.01")

write.csv(as.data.frame(Fig_Chla_dim1$data),
           "/Users/Pierre/Desktop/PhD/COM2LIFE/COM2LIFE_ADN/Chla_vs_dim1.csv")

####__ Fig panel ####

Panel_lake_status<-Fig_mds_param+Fig_Chla_dim1+
  plot_layout(guides = 'collect')+plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size=15,face="bold"))
Panel_lake_status

design_lake<-"
AAC
BBD
BBD"

Panel_lake<- SVG_map_gg+SVG_map_gg+Fig_lake_median_Chla_hor+Fig_mds_param+
  plot_layout(guides = 'collect',design = design_lake)+plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size=15,face="bold"))
Panel_lake


####____________________####
