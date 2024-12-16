
#### Import packages ####

library(BiocManager)
library(phyloseq) 
library(vegan) 
library(tidyverse)
library(tidytext)
library(ggConvexHull)
library(ggrepel)
library(ggh4x)
library(patchwork)
library(rstatix)
library(rcompanion)
library(reshape2)
library(codyn)
library(geosphere)
library(foreign)
library(pairwiseAdonis)
library(scales) 
library(lme4)
library(emmeans)
library(lmerTest)
#BiocManager::install("lmerTest")

#### Rdata ####

save.image(file=".Rdata")

#### Palettes & Strips ####

palette_lake_type_afem<-c("#377EB8","#b19774","#2CB550","#C9DAF8","#2E6D74")

palette_lake_chla_afem<-c("#71992A","#71992A","#E19B8E","#E19B8E","#E19B8E","#E19B8E","#E19B8E","#71992A","#71992A","#377EB8","#377EB8","#71992A")
palette_lake_afem<-c("#71992A","#71992A","#E19B8E","#E19B8E","#E19B8E","#E19B8E","#E19B8E","#377EB8","#377EB8")

# palette_lake_3T<-c("#1E7A36","#28A448",
#                    "#7B6A51","#998465","#b19774","#D1C1AD","#E3D9CD",
#                    "#8CB8E2","#377EB8")


#"#722300",
#palette_lake_3T<-c("#28A448","#1E7A36","#B23600","#CC3E00","#E44600","#F2A380","#F8D1BF","#5C9BD6","#2F6B9D")

tmp<-c("#28A448","#1E7A36","#B23600","#CC3E00",
       "#E44600","#EB7947","#F2AB8C","#5C9BD6","#2F6B9D")

palette_phyto<-c("#ffa07a",
                 "#235B27",
                 "#ffbf00",
                 "#CCEBC5",
                 "#008080",
                 #"#67B16B",
                 "#308238",
                 "#ae773f",
                 #"#ECF8ED",
                 "#67B16B")
pie(rep(1,length(tmp)), col=tmp)

pie(rep(1,length(palette_phyto)), col=palette_phyto)

palette_lake_3T<-c("#28A448","#1E7A36","#980000","#E40000","#E45000","#EB7947","#F2AB8C","#5C9BD6","#2F6B9D")

pie(rep(1,length(wes_gradient)), col=wes_gradient)
pie(rep(1,length(tmp)), col=tmp)
#"#7F300D","#BF4712","#FF5F19","#FFB08D","#FED7C6",
# "#B23600","#E44600","#FF5F19","#FF7E47","#FFA67F",
#"#295D89")
#"#377EB8")



wes_gradient<-c("#3B9AB2","#6FB2C2","#9ABD95","#CEC650","#E8C621","#E49200","#F98400","#EA5C00","#F21A00")

strip_color_lake_type_afem<- strip_themed(
  background_x = elem_list_rect(fill = palette_lake_type_afem,
                                color="black"),
  text_x = elem_list_text(colour = "white",
                          face = "bold",size=12))

strip_color_lake_afem<- strip_themed(
  background_x = elem_list_rect(fill = palette_lake_afem,
                                color="black"),
  text_x = elem_list_text(colour = "white",
                          face = "bold",size=12))

####|####

#### Figures ####


####__ PCA ####

Fig_mds_param <- param_MDS_eucli.df %>%
  ggplot(.,aes(dim1,dim2,color = lake)) +
  geom_convexhull(aes(color = lake,fill=factor(lake),group = lake),alpha=0.4,show.legend = F,size=0.4)+
  geom_segment(data = param_env.df,
               mapping=aes(x = 0, xend = Dim1*3, y = 0, yend = Dim2*3),color = "black",inherit.aes = F,arrow = arrow(length = unit(0.25,"cm"),type = "closed")) +
  geom_text(data = param_env.df, mapping=aes(x = Dim1*3.12, y = Dim2*3.12, label = param),size = 3.5,fill="white",fontface = "bold",inherit.aes = F)+
  geom_point(size=2,shape=16,show.legend = F) + theme_bw() +
  theme(axis.title = element_text(size=12,face = "bold"),
        axis.text = element_blank(),
        legend.title = element_text(size=15,face = "bold"),
        legend.text = element_text(size=12),
        plot.caption = element_blank(),
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
Fig_mds_param

####__ Chla ####

data_chla %>%
  dplyr::group_by(lake) %>%
  dplyr::summarise_at(vars("Chla"),list(mean=mean, sd=sd))


data_chla<-read_delim("chla_2021.csv",";", escape_double = FALSE,
                         trim_ws = TRUE,show_col_types = FALSE)

Fig_lake_median_Chla_hor<-data_chla %>% dplyr::group_by(lake) %>%
  dplyr::mutate(.,
                Chla_median=mean(Chla),
                Status=if_else(lake %in% c("CHA","VER"),"Hypereutrophic",if_else(lake %in% c("JAB","VAI"),"Mesotrophic","Eutrophic")),
                Status=factor(Status,levels=c("Hypereutrophic","Eutrophic","Mesotrophic")),
                lake_name=
                  if_else(lake=="CHA","CSM",
                          if_else(lake=="VER","VSS",
                                  if_else(lake=="GDP","LGP",
                                          if_else(lake=="BOI","BLR",
                                                  if_else(lake=="CTL","CRE",
                                                          if_else(lake=="CRJ1","CER-S",
                                                                  if_else(lake=="CRJ2","CER-L",
                                                                          if_else(lake=="VAI","VSM","JAB")))))))),
                lake_name=factor(lake_name,levels=rev(c("VSS","CSM",
                                                          "LGP","BLR","CRE",
                                                          "CER-S","CER-L",
                                                          "VSM","JAB")))) %>%
  ggplot(.,aes(x=0,y=Chla_median,color=lake_name,label = lake_name, fill = lake_name))+
  geom_point(size=1,show.legend=F,color="white")+theme_bw()+
  geom_hline(yintercept = c(2.6), color = "#2F6B9D", linetype = 2,size=0.8,alpha=0.9)+#E44600
  geom_hline(yintercept = c(7.3), color = "#E45000", linetype = 2,size=0.8,alpha=0.9)+
  geom_hline(yintercept = c(56), color = "#28A448", linetype = 2,size=0.8,alpha=0.9)+
  geom_label(aes(x=-0.14,y=4.1),label="Mesotrophic",size=4.5,color="white",fill="#2F6B9D",fontface = "bold")+
  geom_label(aes(x=-0.14,y=19),label="Eutrophic",size=4.5,color="white",fill="#E45000",fontface = "bold")+
  geom_label(aes(x=-0.14,y=100),label="Hypereutrophic",size=4.5,color="white",fill="#28A448",fontface = "bold")+
  geom_point(size=2.75,mapping=aes(x=0,y=Chla_median,color=lake_name),show.legend = T)+
  geom_label(color = "white",size = 4,fontface = "bold",show.legend = F,
             mapping=aes(label = lake_name, fill = lake_name,
                         x=if_else(lake_name %in% c("BLR","CER-S","VSM","CSM"),-0.05,0.05),
                         y=Chla_median))+
  theme(panel.grid = element_blank(),
        axis.ticks = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_text(size=10),
        legend.title = element_text(size=15,face = "bold"),
        legend.position = "below",
        legend.text = element_text(size=12),#legend.key.size = unit(0.8,"cm"),
        plot.title = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size=15,face = "bold"))+
  scale_color_manual(values=rev(palette_lake_3T),
                     labels = rev(c("Verneuil-sur-Seine (VSS)","Champs-sur-Marne (CSM)",
                                "La Grande-Paroisse (LGP)","Bois-le-Roi (BLR)","Créteil (CRE)",
                                "Cergy - small lake (CER-S)","Cergy - large lake (CER-L)",
                                "Vaires-sur-Marne (VSM)","Jablines (JAB)")))+
  scale_fill_manual(values=rev(palette_lake_3T),
                    labels = rev(c("Verneuil-sur-Seine (VSS)","Champs-sur-Marne (CSM)",
                                   "La Grande-Paroisse (LGP)","Bois-le-Roi (BLR)","Créteil (CRE)",
                                   "Cergy - small lake (CER-S)","Cergy - large lake (CER-L)",
                                   "Vaires-sur-Marne (VSM)","Jablines (JAB)")))+
  scale_x_continuous(limits = c(-0.2,0.2))+
  scale_y_continuous(trans="log10",breaks = c(0,2.6,7.3,56,100),limits = c(NA,150))+
  guides(fill=guide_legend(nrow=3,override.aes=list(shape=22,size=5)),color=F)+
  labs(y=expression(paste(bold("Chlorophyll "),bold(italic("a")),bold(" ["),bold(µg.L^"-1"),bold("]"))),
       fill="Lake")+coord_flip()

Fig_lake_median_Chla_hor

####__ Chla dim1 ####

Fig_Chla_dim1<- Fig_lake_median_Chla_hor$data %>%
  dplyr::mutate(.,samples=paste0(lake,"_",Month_letter)) %>%
  group_by(samples,lake,Month_letter) %>%
  dplyr::summarise(chla_mean=mean(Chla)) %>%
  column_to_rownames(var="samples") %>%
  dplyr::mutate(.,
                Status=if_else(lake %in% c("CHA","VER"),"Hypereutrophic",if_else(lake %in% c("JAB","VAI"),"Mesotrophic","Eutrophic")),
                Status=factor(Status,levels=c("Hypereutrophic","Eutrophic","Mesotrophic")),
                lake_name=
                  if_else(lake=="CHA","Champs-sur-Marne (CSM)",
                          if_else(lake=="VER","Verneuil-sur-Seine (VSS)",
                                  if_else(lake=="GDP","La Grande-Paroisse (LGP)",
                                          if_else(lake=="BOI","Bois-le-Roi (BLR)",
                                                  if_else(lake=="CTL","Créteil (CRE)",
                                                          if_else(lake=="CRJ1","Cergy - small lake (CER-S)",
                                                                  if_else(lake=="CRJ2","Cergy - large lake (CER-L)",
                                                                          if_else(lake=="VAI","Vaires-sur-Marne (VSM)","Jablines (JAB)")))))))),
                lake_name=factor(lake_name,levels=rev(c("Verneuil-sur-Seine (VSS)","Champs-sur-Marne (CSM)",
                                              "La Grande-Paroisse (LGP)","Bois-le-Roi (BLR)","Créteil (CRE)",
                                              "Cergy - small lake (CER-S)","Cergy - large lake (CER-L)",
                                              "Vaires-sur-Marne (VSM)","Jablines (JAB)")))) %>%
  cbind(.,MDS_dim1=Fig_lake_dim1$data$dim1[match(rownames(.),rownames(Fig_lake_dim1$data))]) %>%
  ggplot(.,aes(chla_mean,MDS_dim1,color=lake_name))+theme_bw()+
  #stat_smooth(method = "lm",show.legend = F,mapping=aes(group=1),color="black",alpha=0.15)+
  geom_point(size=3,show.legend = T)+
  #geom_label(mapping=aes(x=1,y=4,fill="white",color="black"),size=5,hjust=0,inherit.aes = F,
          # label=expression(paste0(italic("p"),"<0.01"),R^2,round(summary(lm_gradient_model)[["adj.r.squared"]],2)))+
  annotate(geom="text", x=1, y=4,color="black",size=5,hjust=0,label="paste(italic(p), \"<0.01 \", italic(R) ^ 2, \" 0.74\")", parse = TRUE)+
  theme(panel.grid = element_blank(),
        axis.title.y = element_text(size=12,face = "bold"),
        axis.title.x = element_text(size=14,face = "bold"),
        axis.text = element_text(size=10),
        axis.ticks = element_blank(),
        legend.title = element_text(size=15,face="bold"),
        legend.text = element_text(size=12),
        #legend.key.size = unit(0.5,"cm"),
        legend.position = "bottom",legend.direction = "vertical")+
  scale_color_manual(values=rev(palette_lake_3T))+
  scale_x_continuous(trans="log10",breaks = c(0,5,10,40,100,200))+
  guides(color=guide_legend(nrow=3,override.aes=list(shape=15,size=5)))+
  #scale_y_reverse()+
  labs(x=expression(paste(bold("Chlorophyll "),bold(italic("a")),bold(" ["),bold(µg.L^"-1"),bold("]"))),color="Lake",
       y=paste0("PC1 [44%]"))
Fig_Chla_dim1

spearman_chla_dim1<-cor.test(Fig_Chla_dim1$data$chla_mean,Fig_Chla_dim1$data$MDS_dim1,method = "spearman")
spearman_chla_dim1

lm_gradient_model <- lm(MDS_dim1 ~ Chla_median,data=Fig_Chla_dim1$data)
var_gradient=round(summary(lm_gradient_model)[["adj.r.squared"]],2)

Fig_Chla_dim1

####__ barplot Phyto ####

Phyto_summer<-read_delim("Phyto_summer.csv",";", escape_double = FALSE,
                         trim_ws = TRUE,show_col_types = FALSE) 

Biovol_summer<-read_delim("Biovol_phyto.csv",";", escape_double = FALSE,
                          trim_ws = TRUE,show_col_types = FALSE) 

Fig_Phyto.df <- Phyto_summer %>%
  dplyr::mutate(.,Biovolume_unique=Biovol_summer$Biovolume[match(.$Taxa,Biovol_summer$Taxa)],
                Biovolume_n=Biovolume_unique*nb) %>%
  group_by(lake,month_full,Phylum,Taxa) %>%
  dplyr::summarise(median_BV_nb=median(Biovolume_n)) %>%
  dplyr::group_by(lake,month_full,Phylum) %>%
  dplyr::summarise(sum_BV_nb=sum(median_BV_nb)) %>%
  dplyr::group_by(lake,month_full) %>%
  dplyr::mutate(.,median_rel=(round(sum_BV_nb/sum(sum_BV_nb),4)*100),
                Phylum_legend=if_else(Phylum =="Miozoa","Miozoa (Dinoflagellates)",
                                      if_else(Phylum =="Bacillariophyta","Bacillariophyta (Diatoms)",Phylum)),
                Phylum_domain=if_else(Phylum != "Cyanobacteria","Eukaryotic","Cyanobacteria"),
                lake_month=paste0(lake,"_",month_full),
                month_order=factor(month_full,levels=c("Jun.","Jul.","Aug.","Sep.")),
                lake=
                  if_else(lake=="CHA","CSM",
                          if_else(lake=="VER","VSS",
                                  if_else(lake=="GDP","LGP",
                                          if_else(lake=="BOI","BLR",
                                                  if_else(lake=="CTL","CRE",
                                                          if_else(lake=="CRJ1","CER-S",
                                                                  if_else(lake=="CRJ2","CER-L",
                                                                          if_else(lake=="VAI","VSM","JAB")))))))),
                lake=factor(lake,levels=rev(c("VSS","CSM",
                                              "LGP","BLR","CRE",
                                              "CER-S","CER-L",
                                              "VSM","JAB"))))


write.csv(Fig_Phyto.df %>% subset(month_full=="Jul."),"phyto_july_2021.csv")
# 
# Fig_Phyto<- Fig_Phyto.df %>%
#   ggplot(.,aes(x=month_order, y=median_rel,
#                fill=fct_rev(fct_reorder(Phylum_legend,median_rel))))+
#   geom_col(width=0.8,color="black")+
#   theme_bw()+
#   theme(panel.grid = element_blank(),axis.title.y = element_blank(),
#         axis.title.x= element_blank(),
#         axis.text.x=element_text(size=10,hjust=0.5),
#         axis.ticks =element_blank(),legend.position="right",
#         legend.title = element_text(size=15,face = "bold"),
#         legend.text = element_text(size=12))+
#   scale_x_discrete(expand = c(0.085,0.085))+
#   scale_y_continuous(expand = c(0.02,0.02))+
#   #scale_fill_manual(values=c("#005a32","#238b45","#41ab5d","#6baed6","#74c476","#a1d99b","#c7e9c0","#e5f5e0","#f7fcf5"))+
#   scale_fill_manual(values=palette_phyto)+
#   guides(fill=guide_legend(ncol=3,overide.aes=list(shape=22,size=5)))+
#   labs(fill="Phytoplankton Phyla",y="Relative biovolume")+
#   facet_wrap2(~ lake, scales ="free_x",
#               strip = strip_themed(
#                 background_x = elem_list_rect(fill = rev(palette_lake_3T),
#                                               color = "black"),
#                 text_x = elem_list_text(colour = "white",
#                                         face = "bold",size=9)))

Fig_Phyto_domain<- Fig_Phyto.df %>%
  dplyr::mutate(.,
                Phyla_origin=if_else(Phylum=="Cyanobacteria","Procayotes (Cyanobacteria)","Eucaroytes"),
                month_full=factor(month_full,levels=c("Jun.","Jul.","Aug.","Sep."))) %>%
  dplyr::group_by(lake,month_full,Phyla_origin) %>%
  dplyr::summarise(.,median_Phyla_origin=sum(median_rel)) %>%
  ggplot(.,aes(x=month_full, y=median_Phyla_origin,
               fill=fct_rev(fct_reorder(Phyla_origin,median_Phyla_origin))))+
  geom_col(width=0.8,color="black")+
  theme_bw()+
  theme(panel.grid = element_blank(),axis.title.y = element_blank(),
        axis.title.x= element_blank(),
        axis.text.x=element_text(size=10,hjust=0.5),
        axis.ticks =element_blank(),legend.position="right",
        legend.title = element_text(size=15,face = "bold"),
        legend.text = element_text(size=12))+
  scale_x_discrete(expand = c(0.085,0.085))+
  scale_y_continuous(expand = c(0.02,0.02))+
  scale_fill_manual(values=c("#308238","#CCEBC5"))+
  guides(fill=guide_legend(ncol=3,overide.aes=list(shape=22,size=5)))+
  labs(fill="Phytoplankton domain",y="Relative biovolume")+
  facet_wrap2(~ lake, scales ="free_x",
              strip = strip_themed(
                background_x = elem_list_rect(fill = rev(palette_lake_3T),
                                              color = "black"),
                text_x = elem_list_text(colour = "white",
                                        face = "bold",size=9)))

# mean and sd of phyto domain over the four months for each lake
Phyto_summer %>%
  dplyr::mutate(.,Biovolume_unique=Biovol_summer$Biovolume[match(.$Taxa,Biovol_summer$Taxa)],
                Biovolume_n=Biovolume_unique*nb) %>%
  dplyr::mutate(domain=if_else(Phylum != "Cyanobacteria","Eukaryotic","Cyanobacteria")) %>%
  group_by(lake,month_full,column,domain) %>%
  dplyr::summarise(sum_BV_nb=sum(Biovolume_n)) %>%
  dplyr::group_by(lake,month_full,column) %>%
  dplyr::mutate(.,sum_rel=(round(sum_BV_nb/sum(sum_BV_nb),4)*100)) %>%
  group_by(lake,domain) %>%
  dplyr::summarise_at(vars("sum_rel"),list(mean=mean, sd=sd))

# mean and sd of phyto genus
Fig_genus$data %>%  subset(lake=="VSS") %>%
  subset(.,Cyano_legend %in% c("Ceratium")) %>%
  subset(month_full %in% c('Jul.','Aug.')) %>%
  group_by(lake) %>%
  dplyr::summarise_at(vars("median_cyano_legend"),list(mean=mean, sd=sd))

Fig_genus$data %>%  subset(lake=="CSM") %>%
  subset(.,Cyano_legend %in% c("Aphanizomenon")) %>%
  subset(month_full %in% c('Jun.','Jul.','Aug.')) %>%
  group_by(lake) %>%
  dplyr::summarise_at(vars("median_cyano_legend"),list(mean=mean, sd=sd))

Fig_genus$data %>%  subset(lake=="BLR") %>%
  subset(.,Cyano_legend %in% c("Cyanocatena-like")) %>%
  group_by(lake) %>%
  dplyr::summarise_at(vars("median_cyano_legend"),list(mean=mean, sd=sd))

Fig_genus$data %>%  subset(lake=="CSM") %>%
  subset(.,Cyano_legend %in% c("Dolichospermum")) %>%
  subset(month_full %in% c('Jun.','Jul.','Aug.')) %>%
  group_by(lake) %>%
  dplyr::summarise_at(vars("median_cyano_legend"),list(mean=mean, sd=sd))

Fig_Phyto_domain$data %>%
  subset(lake=="CSM") %>%
  subset(.,Phyla_origin %in% c("Procayotes (Cyanobacteria)")) %>%
  group_by(lake,month_full) %>%
  dplyr::summarise_at(vars("median_Phyla_origin"),list(mean=mean, sd=sd))


Fig_Phyto$data %>%  subset(lake=="VSS") %>%
  subset(.,Taxa %in% c("Aphanizomenon","Aphanizomenon_gracile")) %>%
  group_by(month_full) %>%
  dplyr::summarise(sum=sum(median_rel)) %>%
  group_by(lake) %>%
  dplyr::summarise_at(vars("sum"),list(median=median, sd=sd))
  
Fig_Cyano$data %>%  subset(lake=="VSS") %>%
  subset(.,Cyano_legend %in% c("Ceratium")) %>%
  subset(month_full %in% c('Jul.','Aug.')) %>%
  group_by(lake) %>%
  dplyr::summarise_at(vars("median_cyano_legend"),list(median=median, sd=sd))
View(Fig_Cyano.df)

Fig_Cyano.df <- Phyto_summer %>%
  dplyr::mutate(.,Biovolume_unique=Biovol_summer$Biovolume[match(.$Taxa,Biovol_summer$Taxa)],
                Biovolume_n=Biovolume_unique*nb) %>%
  group_by(lake,month_full,Phylum,Taxa) %>%
  dplyr::summarise(median_BV_nb=median(Biovolume_n)) %>%
  dplyr::group_by(lake,month_full) %>%
  dplyr::mutate(.,median_rel=(round(median_BV_nb/sum(median_BV_nb),4)*100)) %>%
  subset(.,Phylum=="Cyanobacteria" | Taxa=="Ceratium_hirundinella") %>% replace(is.na(.), "Others") %>%
  dplyr::mutate(.,
                lake_month=paste0(lake,"_",month_full),
                month_full=factor(month_full,levels=c("Jun.","Jul.","Aug.","Sep.")),
                Cyano_legend=
                  if_else(Taxa=="Aphanizomenon","Aphanizomenon",
                          if_else(Taxa=="Aphanizomenon_gracile","Aphanizomenon",
                                  if_else(Taxa=="Dolichospermum","Dolichospermum",
                                                          if_else(Taxa=="Cyanocatena","Cyanocatena-like",
                                                                  if_else(Taxa %in% c("Ceratium_hirundinella","Ceratium"),"Ceratium",
                                                                  "Others"))))),
                lake=
                  if_else(lake=="CHA","CSM",
                          if_else(lake=="VER","VSS",
                                  if_else(lake=="GDP","LGP",
                                          if_else(lake=="BOI","BLR",
                                                  if_else(lake=="CTL","CRE",
                                                          if_else(lake=="CRJ1","CER-S",
                                                                  if_else(lake=="CRJ2","CER-L",
                                                                          if_else(lake=="VAI","VSM","JAB")))))))),
                lake=factor(lake,levels=rev(c("VSS","CSM",
                                              "LGP","BLR","CRE",
                                              "CER-S","CER-L",
                                              "VSM","JAB")))) %>%
  dplyr::group_by(lake,month_full,Cyano_legend) %>%
  dplyr::summarise(.,median_cyano_legend=sum(median_rel))
  
Fig_genus<-Fig_Cyano.df %>%
  ggplot(.,aes(x=month_full, y=median_cyano_legend,
               fill=fct_rev(fct_reorder(Cyano_legend,median_cyano_legend))))+
  geom_col(width=0.8,color="black")+
  theme_bw()+
  theme(panel.grid = element_blank(),axis.title.y = element_blank(),
        axis.title.x= element_blank(),
        axis.text.x=element_text(size=10,hjust=0.5),
        axis.ticks =element_blank(),legend.position="right",
        legend.title = element_text(size=15,face = "bold"),
        legend.text.align = 0,
        legend.text = element_text(size=12))+
  scale_x_discrete(expand = c(0.085,0.085))+
  scale_y_continuous(expand = c(0.02,0.02),breaks = c(0,25,50,75,100),limits = c(0,100))+
  scale_fill_manual(values=c("#CCEBC5","#308238","#67B16B","lightgrey","#008080"),
  #scale_fill_manual(values=c("#173C5B","#225A88","#245D8E","#2E78B6","#97BCDB","#CBDDEC"),
                    labels =  c(expression(italic("Cyanocatena")),expression(italic("Ceratium")),
                               expression(italic("Aphanizomenon")),"Others",
                               expression(italic("Dolichospermum"))))+
  guides(fill=guide_legend(nrow=3,shape=22,size=5))+
  labs(fill="Phytoplankton Genus",y="Relative biovolume")+
  facet_wrap2(~ lake, scales ="free_x",
              strip = strip_themed(
                background_x = elem_list_rect(fill = rev(palette_lake_3T),
                                              color = "black"),
                text_x = elem_list_text(colour = "white",
                                        face = "bold",size=10)))
View(Fig_genus$data)          

Phyto_summer %>% subset(lake=="VER") %>%
  dplyr::mutate(.,Biovolume_unique=Biovol_summer$Biovolume[match(.$Taxa,Biovol_summer$Taxa)],
                Biovolume_n=Biovolume_unique*nb) %>%
  group_by(lake,month_full,Taxa) %>%
  dplyr::summarise(median_BV_nb=median(Biovolume_n)) %>%
  dplyr::group_by(lake,month_full) %>%
  dplyr::mutate(.,median_rel=(round(median_BV_nb/sum(median_BV_nb),4)*100)) %>%
  dplyr::group_by(month_full,Taxa) %>%
  dplyr::summarise_at(vars("median_rel"),list(taxa_sum=sum)) %>%
  dplyr::group_by(Taxa) %>%
  dplyr::summarise_at(vars("taxa_sum"),list(median=median, sd=sd))

####__ Corr PCA ####
corr_Fig<- melted_corr %>%  subset(.,Var2=="Dim.1") %>%
  ggplot(.,aes(x = Var2, y = fct_reorder(Var1,value), fill = value))+
  geom_point(color='black',shape=21,size=15) + theme_bw()+
  geom_text(aes(x = Var2, y = Var1, label = round(value, 2)), color = "white",size = 3.8)+
  theme(panel.grid=element_blank(),
        axis.title.x = element_text(face="bold",size = 15,hjust=0.5),
        axis.text.y = element_text(face="bold",size = 12),
        axis.text.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks = element_blank(),
        legend.title = element_text(face="bold",size = 10))+
  scale_fill_gradient2(midpoint = 0.5, mid ="grey70", low = "blue", high = "red",
                       limits = c(0, +1)) +
  labs(fill = expression(Cos^"2"),size="",x="PC1 - [29%]")
corr_Fig

####__ Core ####
Fig_shared <-
  rbind(shared_ASV.df[,-5],shared_KO.df) %>%
  dplyr::mutate(.,
                rel.abund=as.numeric(rel.abund),
                per_lake=as.numeric(per_lake),
                per_lake=round(per_lake,digits = 1),
                lake=
                  if_else(lake=="Champs","Champs-sur-Marne (CSM)",
                          if_else(lake=="Verneuil","Verneuil-sur-Seine (VSS)",
                                  if_else(lake=="Grande-Paroisse","La Grande-Paroisse (LGP)",
                                          if_else(lake=="Bois","Bois-le-Roi (BLR)",
                                                  if_else(lake=="Créteil","Créteil (CRE)",
                                                          if_else(lake=="Cergy-small","Cergy - small lake (CER-S)",
                                                                  if_else(lake=="Cergy-large","Cergy - large lake (CER-L)",
                                                                          if_else(lake=="Vaires","Vaires-sur-Marne (VSM)","Jablines (JAB)")))))))),
                lake=factor(lake,levels=rev(c("Verneuil-sur-Seine (VSS)","Champs-sur-Marne (CSM)",
                                              "La Grande-Paroisse (LGP)","Bois-le-Roi (BLR)","Créteil (CRE)",
                                              "Cergy - small lake (CER-S)","Cergy - large lake (CER-L)",
                                              "Vaires-sur-Marne (VSM)","Jablines (JAB)")))) %>%
  ggplot(.,aes(x=lake,y=features,size = rel.abund,label=per_lake,color=lake)) +
  geom_point(#color=c("#28A448","#1E7A36","#980000","#E40000","#E45000","#EB7947","#F2AB8C","#5C9BD6","#2F6B9D",
                     #"#28A448","#1E7A36","#980000","#E40000","#E45000","#EB7947","#F2AB8C","#5C9BD6","#2F6B9D"),
             show.legend = T)+
  geom_text(size=4,color="white",show.legend = F,fontface="bold")+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.ticks = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        legend.title = element_text(size=15,face = "bold"),
        legend.text = element_text(size=12),
        legend.position = "bottom",
        legend.direction = "vertical",
        axis.title.y = element_text(size=14,face = "bold",angle=90,hjust=0.5),
        axis.title.x = element_blank())+
  scale_color_manual(values=c(rev(palette_lake_3T),rev(palette_lake_3T)))+
  scale_size_continuous(range = c(6,17),breaks = c(30,50,90),limits = c(30,NA),labels = c("30%","50%","90%"))+
  guides(color=guide_legend(override.aes = list(shape=15,size=5),nrow=3),fill=F,shape=F,text=F,size=guide_legend(nrow=1))+
  labs(size="Relative abundance",y="ASVs   KOs",color="Lake")
Fig_shared

Fig_shared <- Fig_shared +theme(legend.position = "bottom",legend.direction = "vertical")+
  guides(size = guide_legend(nrow=1))

Fig_shared$data %>%  dplyr::mutate(.,status_features=paste0(status,"_",features),
                                   nb_features=as.numeric(nb_features)) %>%
  dplyr::group_by(status_features) %>%
  dplyr::summarise_at(vars("nb_features"),list(median=median, sd=sd))

Fig_shared$data %>% dplyr::mutate(.,status_features=paste0(status,"_",features),
                                  nb_features=as.numeric(nb_features)) %>%
  dplyr::group_by(status_features) %>%
  dplyr::summarise_at(vars("per_lake"),list(median=median, sd=sd))

Fig_shared$data %>%  dplyr::mutate(.,
                                   status=as.character(status),
                                   features=as.character(features),
                                   status_features=paste0(status,"_",features),
                                   status_features=as.character(status_features),
                                   nb_features=as.numeric(nb_features),
                                   per_lake=as.numeric(per_lake)) %>%
  dplyr::group_by(status_features) %>%
  dplyr::summarise_at(vars("nb_features"),list(median=median, sd=sd))
View(Fig_shared$data)

####__ Turnover ####

ASV_KO_turnover.df<- read.delim("ASV_KO_turnover_box.csv",sep = ";") %>%
  dplyr::mutate(.,
                lake=factor(lake,levels=c("CHA","VER","GDP","BOI","CTL","CRJ1","CRJ2","VAI","JAB")),
                Status=factor(Status,levels=c("Mesotrophic","Eutrophic","Hypereutrophic")),
                Status_features=factor(Status_features,levels=c("Mesotrophic_KOs","Mesotrophic_ASVs",
                                                                "Eutrophic_KOs","Eutrophic_ASVs",
                                                                "Hypereutrophic_KOs","Hypereutrophic_ASVs")))
ASV_KO_turnover.stats<-ASV_KO_turnover.df %>%
  dplyr::mutate(.,lake_features=paste0(lake,"_",features)) %>%
  wilcox_test(turnover ~ Status_features,p.adjust.method = "BH") %>%
  add_significance("p.adj")  %>% add_xy_position(x = "Status_features", dodge = 0.8) %>%
  mutate(g1=gsub("_","",as.character(.$group1)),g2=gsub("_","",as.character(.$group2))) %>%
  unite("comp",g1:g2,sep=" - ",remove = FALSE)

letters_ASV_KO_turnover<-tibble(Status_features=unique(sort(ASV_KO_turnover.df$Status_features)),
                              letters=cldList(p.adj ~ as.character(comp),ASV_KO_turnover.stats,threshold  = 0.05)$Letter)

letters_ASV_KO_turnover$letters<-c("a","b","a","b","a","c")

Fig_box_ASV_KO_turnover<- ASV_KO_turnover.df %>% 
  dplyr::group_by(Status_features) %>%
  dplyr::mutate(.,median_turnover=median(turnover)) %>%
  cbind(.,letter=letters_ASV_KO_turnover$letters[match(.$Status_features,letters_ASV_KO_turnover$Status_features)]) %>%
  ggplot(.,aes(x=Status_features,y=turnover,fill=features))+theme_bw()+
  geom_boxplot(width=0.6,color="black",show.legend = F)+
  geom_point(size=1,color="black",alpha=0.5,show.legend = F)+
  geom_text(mapping=aes(y=1.02*max(turnover),
                        label=letter),vjust=0,size=4.5,color="black",fontface="bold")+
  theme(panel.grid = element_blank(),plot.title = element_text(hjust=0.5,size=15,face = "bold"),
        axis.ticks = element_blank(),axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size=10,hjust=0.5),axis.text.y = element_text(size=12),
        legend.title = element_text(size=15,face = "bold"),legend.text = element_text(size=12))+
  scale_fill_manual(values = c("white","lightgrey"))+
  scale_y_continuous(limits =c(0,1.08*max(ASV_KO_turnover.df$turnover)),position = "right",
                     labels = scales::label_number(accuracy = 0.1))+
  scale_x_discrete(labels=c('KOs', 'ASVs',
                            'KOs', 'ASVs',
                            'KOs', 'ASVs'))+
  #labs(fill="Turnover")+
  facet_wrap2(~ Status,scale="free_x",nrow = 1,strip.position = "top",
              strip = strip_color_lake_type_afem<- strip_themed(
                background_x = elem_list_rect(fill = c("#2F6B9D","#E45000","#28A448"),
                                              color="black"),
                text_x = elem_list_text(colour = "white",
                                        face = "bold",size=10)))

Fig_box_ASV_KO_turnover$data %>%
  dplyr::group_by(Status_features) %>%
  dplyr::summarise_at(vars("turnover"),list(mean=mean, sd=sd))

Fig_box_ASV_KO_turnover$data %>% dplyr::group_by(features) %>%
  dplyr::summarise_at(vars("turnover"),list(median=median, sd=sd))

testlm<-lmer(turnover~ Status + month_number + (1|lake),data=subset(Fig_box_ASV_KO_turnover$data,features=="KOs"),REML = F)
summary(testlm)

####__ barplot ####
write.csv(A_D_rarW@sam_data,"metadata2021.csv")

Barplot_AFEM <- A_D_rarW %>% tax_glom(.,taxrank ="Phylum") %>%
  psmelt(.) %>% dplyr::group_by(lake_month,Phylum,.drop = FALSE,.add = TRUE) %>%
  summarise(median_count=median(Abundance),.groups = "keep") %>%
  dplyr::group_by(lake_month) %>%
  mutate(median_abundance=as.numeric(
    paste0((round(median_count/sum(median_count),4))*100))) %>%
  cbind(lake_order=
          (A_D_rarW@sam_data$lake_order[match(.$lake_month,
                                              A_D_rarW@sam_data$lake_month)]),
        month=
          (A_D_rarW@sam_data$month[match(.$lake_month,
                                         A_D_rarW@sam_data$lake_month)]),
        month_full=
          (A_D_rarW@sam_data$month_full[match(.$lake_month,
                                         A_D_rarW@sam_data$lake_month)]),
        lake_type=
          (A_D_rarW@sam_data$lake_type[match(.$lake_month,
                                             A_D_rarW@sam_data$lake_month)]),
        .) %>%
  dplyr::mutate(.,
                Phylum_legend=ifelse(median_abundance<1,"Phylum < 1%",Phylum),
                Phylum_legend=ifelse(Phylum_legend %in%
                                       c("Armatimonadota","Campilobacterota","SAR324_clade(Marine_group_B)","Bdellovibrionota","Desulfobacterota","Phylum < 1%"),
                                     "Others",Phylum_legend),
                lake_order=
                  if_else(lake_order=="CHA","CSM",
                          if_else(lake_order=="VER","VSS",
                                  if_else(lake_order=="GDP","LGP",
                                          if_else(lake_order=="BOI","BLR",
                                                  if_else(lake_order=="CTL","CRE",
                                                          if_else(lake_order=="CRJ1","CER-S",
                                                                  if_else(lake_order=="CRJ2","CER-L",
                                                                          if_else(lake_order=="VAI","VSM","JAB")))))))),
                lake_order=factor(lake_order,levels=rev(c("VSS","CSM",
                                              "LGP","BLR","CRE",
                                              "CER-S","CER-L",
                                              "VSM","JAB"))),
                month_short=if_else(month_full == "June","Jun.",
                                    if_else(month_full == "July","Jul.",
                                            if_else(month_full == "August","Aug.","Sep."))),
                month_short=factor(month_short,levels=c("Jun.","Jul.","Aug.","Sep."))) %>%
  ggplot(.,aes(x=month_short, y=(median_abundance/100),
               fill=fct_rev(fct_reorder(Phylum_legend,median_abundance))))+
  geom_col(width=0.8,color="black",size=0.2)+
  theme_bw()+theme(panel.grid = element_blank(),axis.title.y = element_blank(),
                   axis.title.x = element_blank(),
                   axis.text.x=element_text(size=10,hjust=0.5),
                   axis.ticks =element_blank(),
                   legend.title = element_text(size=15,face = "bold"),
                   legend.position = "bottom",
                   legend.direction = "vertical",
                   legend.text = element_text(size=12))+
  scale_x_discrete(expand = c(0.085,0.085))+
  scale_y_continuous(expand = c(0,0),position = "right",
                     labels = scales::label_number(accuracy = 0.1))+
  scale_fill_manual(values=c("#B3CDE3","#CCEBC5","#DECBE4","#FED9A6","#FFFFCC","#E5D8BD","#FDDAEC","#F2F2F2","firebrick"))+
  guides(fill=guide_legend(nrow=3,shape=22,size=5))+
  labs(fill="Prokaryotic Phyla")+
  facet_wrap2(~ lake_order,scale="fixed",
              strip = strip_color_lake<- strip_themed(
                background_x = elem_list_rect(fill = rev(palette_lake_3T),
                                              color = "black"),
                text_x = elem_list_text(colour = "white",face = "bold",size=9)))

Barplot_AFEM
unique(Barplot_AFEM$data$Phylum_legend)
Barplot_AFEM$data %>%
  dplyr::group_by(Phylum_legend) %>%
  dplyr::summarise_at(vars("median_abundance"),list(mean=mean, sd=sd))

subset_samples(A_D_rarW, lake_order=="BOI") %>% subset_taxa(.,Genus=="Cyanobium_PCC-6307")

barpot_cyanobium.df<-A_D_rarW %>%
  psmelt(.) %>% dplyr::group_by(lake_month,Genus,OTU,.drop = FALSE,.add = TRUE) %>%
  summarise(median_count=median(Abundance),.groups = "keep") %>%
  dplyr::group_by(lake_month,Genus) %>%
  dplyr::summarise(median_count=sum(median_count)) %>%
  dplyr::group_by(lake_month) %>%
  mutate(median_abundance=as.numeric(
    paste0((round(median_count/sum(median_count),4))*100))) %>%
  cbind(lake_order=
          (A_D_rarW@sam_data$lake_order[match(.$lake_month,
                                              A_D_rarW@sam_data$lake_month)]),
        month=
          (A_D_rarW@sam_data$month[match(.$lake_month,
                                         A_D_rarW@sam_data$lake_month)]),
        month_full=
          (A_D_rarW@sam_data$month_full[match(.$lake_month,
                                              A_D_rarW@sam_data$lake_month)]),
        lake_type=
          (A_D_rarW@sam_data$lake_type[match(.$lake_month,
                                             A_D_rarW@sam_data$lake_month)]),
        .)  %>% subset(Genus =="Cyanobium_PCC-6307") 



seq_cyanobium<-A_D_rarW %>%
  psmelt(.) %>%  dplyr::group_by(lake_month,lake,Genus,OTU,.drop = FALSE,.add = TRUE) %>%
  summarise(median_count=median(Abundance),.groups = "keep") %>%
  subset(Genus =="Cyanobium_PCC-6307") 

list_Cyanobium_BLR <- unique((seq_cyanobium %>% subset(lake =="BOI" & median_count >0))$OTU)

length(unique(seq_cyanobium$OTU))

View((seq_cyanobium %>% subset(lake !="BOI" & median_count >0)))

Fig_cyanobium<-
  barpot_cyanobium.df %>%
  dplyr::mutate(.,
                lake_order=
                  if_else(lake_order=="CHA","CSM",
                          if_else(lake_order=="VER","VSS",
                                  if_else(lake_order=="GDP","LGP",
                                          if_else(lake_order=="BOI","BLR",
                                                  if_else(lake_order=="CTL","CRE",
                                                          if_else(lake_order=="CRJ1","CER-S",
                                                                  if_else(lake_order=="CRJ2","CER-L",
                                                                          if_else(lake_order=="VAI","VSM","JAB")))))))),
                lake_order=factor(lake_order,levels=rev(c("VSS","CSM",
                                              "LGP","BLR","CRE",
                                              "CER-S","CER-L",
                                              "VSM","JAB"))),
                month_short=if_else(month_full == "June","Jun.",
                                    if_else(month_full == "July","Jul.",
                                            if_else(month_full == "August","Aug.","Sep."))),
                month_short=factor(month_short,levels=c("Jun.","Jul.","Aug.","Sep."))) %>%
  ggplot(.,aes(x=month_short, y=(median_abundance), fill=Genus))+
  geom_col(width=0.5,color="black",size=0.2,show.legend = F)+
  theme_bw()+theme(panel.grid = element_blank(),axis.title.y = element_blank(),
                   axis.title.x = element_blank(),
                   axis.text.x=element_text(size=10,hjust=0.5),
                   axis.ticks =element_blank(),
                   legend.position="bottom",
                   legend.direction ="vertical",
                   legend.title = element_text(size=15,face = "bold"),
                   legend.text = element_text(size=12))+
  scale_x_discrete(expand = c(0.085,0.085))+
  scale_y_continuous(expand = c(0.02,0.02),limits = c(0,100))+
  scale_fill_manual(values=c("#CCEBC5"))+
  guides(fill=guide_legend(nrow=1,shape=22,size=5),strip=T)+
  labs(fill="Microbial Genus")+
  facet_wrap2(~ lake_order,scale="fixed",nrow = 3,
              strip = strip_color_lake<- strip_themed(
                background_x = elem_list_rect(fill = rev(palette_lake_3T),
                                              color = "black"),
                text_x = elem_list_text(colour = "white",face = "bold",size=9)))
Fig_cyanobium

ggsave("/Users/pierre/Desktop/ms1/Figures/to_inkscape/FS6.pdf", height = 9.8, width = 13.4, units = "in",dpi = "retina")

designS5="
A
A
B"

Fig_simper_summer+Fig_cyanobium+
  plot_layout(design = designS5)+
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size=15,face="bold"),
        legend.position = "bottom",legend.direction = "vertical")   


barpot_cyanobium.df %>% dplyr::mutate(.,BOI=if_else(lake_order=="BOI","BOI","Others")) %>%
  dplyr::group_by(BOI) %>%
  dplyr::summarise_at(vars("median_abundance"),list(mean=mean, sd=sd))
View(barpot_cyanobium)

####|####

####__ KOs BC ####

Fig_hellinger_ko <-
  Ko_MDS_hellinger.df %>%
  dplyr::mutate(.,
                lake=
  if_else(lake=="CHA","Champs-sur-Marne (CSM)",
          if_else(lake=="VER","Verneuil-sur-Seine (VSS)",
                  if_else(lake=="GDP","La Grande-Paroisse (LGP)",
                          if_else(lake=="BOI","Bois-le-Roi (BLR)",
                                  if_else(lake=="CTL","Créteil (CRE)",
                                          if_else(lake=="CRJ1","Cergy - small lake (CER-S)",
                                                  if_else(lake=="CRJ2","Cergy - large lake (CER-L)",
                                                          if_else(lake=="VAI","Vaires-sur-Marne (VSM)","Jablines (JAB)")))))))),
  lake=factor(lake,levels=rev(c("Verneuil-sur-Seine (VSS)","Champs-sur-Marne (CSM)",
                                "La Grande-Paroisse (LGP)","Bois-le-Roi (BLR)","Créteil (CRE)",
                                "Cergy - small lake (CER-S)","Cergy - large lake (CER-L)",
                                "Vaires-sur-Marne (VSM)","Jablines (JAB)")))) %>%
  ggplot(.,aes(dim1,dim2,
               #shape=month_full,
               color = lake)) +
  geom_convexhull(aes(color = lake,fill=lake,
                      group =
                        factor(lake,levels=
                                 c("Verneuil-sur-Seine (VSS)","Champs-sur-Marne (CSM)",
                                   "La Grande-Paroisse (LGP)","Bois-le-Roi (BLR)","Créteil (CRE)",
                                   "Cergy - small lake  (CER-S)","Cergy - large lake (CER-L)",
                                   "Vaires-sur-Marne (VSM)","Jablines (JAB)"))),
                  alpha=0.4,show.legend = F,size=0.4)+
  geom_point(size=2,show.legend = F,
             aes(group =
                   factor(lake,levels=c("Verneuil-sur-Seine (VSS)","Champs-sur-Marne (CSM)",
                                        "La Grande-Paroisse (LGP)","Bois-le-Roi (BLR)","Créteil (CRE)",
                                        "Cergy - small lake  (CER-S)","Cergy - large lake (CER-L)",
                                        "Vaires-sur-Marne (VSM)","Jablines (JAB)")))) + theme_bw() +
  # geom_label(Ko_centroids_hellinger.df,mapping=aes(dim1,dim2,label = lake, fill = lake),
  #            color = "white",size = 4,fontface = "bold",show.legend = F,inherit.aes = F)+
  theme(panel.grid = element_blank(),
        plot.title = element_blank(),
        axis.text = element_text(size=10),
        axis.title = element_text(size=12,face = "bold"),
        axis.ticks = element_blank(),
        legend.title = element_text(size=15,face = "bold",hjust=0),
        legend.text = element_text(size=12),
        legend.position = "bottom",legend.direction = "vertical")+
  scale_color_manual(values=rev(palette_lake_3T))+
  scale_fill_manual(values=rev(palette_lake_3T))+
  scale_y_continuous(trans = 'reverse',breaks = c(-0.05,0,0.05))+
  scale_x_continuous(trans = 'reverse',breaks = c(-0.05,0,0.05,0.1))+
  guides(color = guide_legend(override.aes=list(size=5,shape=15),nrow=3,order = 1))+
  labs(color="Lake",#shape="Month",
       x=paste0("Axis 1 [",round(MDS_FUNC.hellinger.mds$eig[1]*100/sum(MDS_FUNC.hellinger.mds$eig),1),"%]"),
       y=paste0("Axis 2 [",round(MDS_FUNC.hellinger.mds$eig[2]*100/sum(MDS_FUNC.hellinger.mds$eig),1),"%]"))

Fig_hellinger_ko

####__ ASVs BC ####

Fig_BC_summer<-
  BC_AD_rarW$data %>%
  dplyr::mutate(.,
                lake=
                  if_else(lake=="CHA","Champs-sur-Marne (CSM)",
                          if_else(lake=="VER","Verneuil-sur-Seine (VSS)",
                                  if_else(lake=="GDP","La Grande-Paroisse (LGP)",
                                          if_else(lake=="BOI","Bois-le-Roi (BLR)",
                                                  if_else(lake=="CTL","Créteil (CRE)",
                                                          if_else(lake=="CRJ1","Cergy - small lake (CER-S)",
                                                                  if_else(lake=="CRJ2","Cergy - large lake (CER-L)",
                                                                          if_else(lake=="VAI","Vaires-sur-Marne (VSM)","Jablines (JAB)")))))))),
                lake=factor(lake,levels=rev(c("Verneuil-sur-Seine (VSS)","Champs-sur-Marne (CSM)",
                                              "La Grande-Paroisse (LGP)","Bois-le-Roi (BLR)","Créteil (CRE)",
                                              "Cergy - small lake (CER-S)","Cergy - large lake (CER-L)",
                                              "Vaires-sur-Marne (VSM)","Jablines (JAB)")))) %>%
  ggplot(.,aes(Axis.1,Axis.2,
               #shape=month_full
               color=lake))+
  geom_convexhull(aes(color = lake,fill=lake,
                      group =
                        factor(lake,levels=
                                 c("Verneuil-sur-Seine (VSS)","Champs-sur-Marne (CSM)",
                                   "La Grande-Paroisse (LGP)","Bois-le-Roi (BLR)","Créteil (CRE)",
                                   "Cergy - small lake (CER-S)","Cergy - large lake (CER-L)",
                                   "Vaires-sur-Marne (VSM)","Jablines (JAB)"))),
                  alpha=0.4,show.legend = F,size=0.4)+
  geom_point(size=2,show.legend = F,
             aes(group =
                   factor(lake,levels=c("Verneuil-sur-Seine (VSS)","Champs-sur-Marne (CSM)",
                                        "La Grande-Paroisse (LGP)","Bois-le-Roi (BLR)","Créteil (CRE)",
                                        "Cergy - small lake (CER-S)","Cergy - large lake (CER-L)",
                                        "Vaires-sur-Marne (VSM)","Jablines (JAB)")))) + theme_bw() +
  # geom_label(BC_AD_centroids.df,mapping=aes(Axis.1,Axis.2,label = lake, fill = lake),
  #            color = "white",size = 4,fontface = "bold",show.legend = F,inherit.aes = F)+
  theme(panel.grid = element_blank(),
        plot.title = element_blank(),
        axis.title = element_text(size=12,face = "bold"),
        axis.text = element_text(size=10),axis.ticks = element_blank(),
        legend.title = element_text(size=15,face = "bold"),
        legend.text = element_text(size=12),
        legend.position = "bottom",legend.direction = "vertical")+
  scale_color_manual(values=rev(palette_lake_3T))+
  scale_fill_manual(values=rev(palette_lake_3T))+
  guides(color = guide_legend(override.aes=list(size=5,shape=15)),order=1,nrow=3)+
  scale_y_continuous(breaks = c(-0.4,-0.2,0,0.2,0.4))+
  scale_x_continuous(trans = 'reverse',breaks = c(0.2,0,-0.2,-0.4))+
  labs(color="Lake",
       x=paste0("Axis 1 [",round(A_D_rarW.BC$values$Relative_eig[1],3)*100,"%]"),
       y=paste0("Axis 2 [",round(A_D_rarW.BC$values$Relative_eig[2],3)*100,"%]"))

Fig_BC_summer

####__ ASVs BC without Cyanobacteria ####

Fig_BC_summer_noCyano<-
  BC_AD_rarW_noCyano$data %>%
  dplyr::mutate(.,
                lake=
                  if_else(lake=="CHA","CSM",
                          if_else(lake=="VER","VSS",
                                  if_else(lake=="GDP","LGP",
                                          if_else(lake=="BOI","BLR",
                                                  if_else(lake=="CTL","CRE",
                                                          if_else(lake=="CRJ1","CER-S",
                                                                  if_else(lake=="CRJ2","CER-L",
                                                                          if_else(lake=="VAI","VSM","JAB")))))))),
                lake=factor(lake,levels=rev(c("VSS","CSM",
                                              "LGP","BLR","CRE",
                                              "CER-S","CER-L",
                                              "VSM","JAB")))) %>%
  ggplot(.,aes(Axis.1,Axis.2,
               #shape=month_full
               color=lake))+
  geom_convexhull(aes(color = lake,fill=lake,
                      group =
                        factor(lake,levels=
                                 c("VSS","CSM",
                                   "LGP","BLR","CRE",
                                   "CER-S","CER-L",
                                   "VSM","JAB"))),
                  alpha=0.4,show.legend = F,size=0.4)+
  geom_point(size=2,show.legend = T,
             aes(group =
                   factor(lake,levels=c("VSS","CSM",
                                        "LGP","BLR","CRE",
                                        "CER-S","CER-L",
                                        "VSM","JAB")))) + theme_bw() +
  # geom_label(BC_AD_centroids.df,mapping=aes(Axis.1,Axis.2,label = lake, fill = lake),
  #            color = "white",size = 4,fontface = "bold",show.legend = F,inherit.aes = F)+
  theme(panel.grid = element_blank(),aspect.ratio = 1,
        #plot.title= element_text(size=11.5,hjust=0),
        plot.title = element_blank(),
        axis.title = element_text(size=12,face = "bold"),
        axis.text = element_text(size=10),axis.ticks = element_blank(),
        legend.title = element_text(size=15,face = "bold"),
        legend.text = element_text(size=12),
        legend.position = "bottom",legend.direction = "vertical")+
  scale_color_manual(values=rev(palette_lake_3T))+
  scale_fill_manual(values=rev(palette_lake_3T))+
  guides(color = guide_legend(override.aes=list(size=5,shape=15),nrow = 3),fill=FALSE)+
  #scale_x_reverse()+
  labs(color="Lake",
       x=paste0("Axis 1 [",round(A_D_rarW_noCyano.BC$values$Relative_eig[1],3)*100,"%]"),
       y=paste0("Axis 2 [",round(A_D_rarW_noCyano.BC$values$Relative_eig[2],3)*100,"%]"))

Fig_BC_summer_noCyano

ggsave("/Users/pierre/Desktop/ms1/Figures/to_inkscape/FS7.pdf", height = 9.8, width = 13.4, units = "in",dpi = "retina")

####|####

####__ Distance-decay ####

####______ KOs ####

spearman_KO_geo<-cor.test(Hellinger_status.df$hellinger_dist,Hellinger_status.df$Havers_dist,method = "spearman")

var_KO_geo=paste0("p>0.05")

Fig_Hellinger_Geo<-Hellinger_status.df %>%
  dplyr::mutate(.,Havers_dist=(Havers_dist/1000)) %>%
  ggplot(.,aes(Havers_dist,hellinger_dist,group=1))+
  geom_point(size=2,shape=21,fill="grey",color="black")+
  stat_smooth(method = "lm",color="#4166F5",fill="lightblue",show.legend = F)+
  theme_classic()+
  theme(axis.ticks = element_blank(),
        axis.text = element_text(size=12),
        axis.title=element_text(size=12,face="bold"),aspect.ratio = 1,
        plot.caption = element_text(size=12,face="italic"))+
  #scale_color_manual(values=palette_slope)+
  #scale_fill_manual(values=palette_lake_type)+
  labs(x="Distance between lakes (km)",caption=var_KO_geo,
       y="Gene-content dissimilarity")
Fig_Hellinger_Geo

####______ ASVs ####

spearman_BC_geo<-cor.test(BC_Geo_Status.df$BC,BC_Geo_Status.df$Havers_dist,method = "spearman")

var_BC_geo=paste0("p<0.01 R2 0.22")

Fig_BC_Geo<-BC_Geo_Status.df %>%
  dplyr::mutate(.,Havers_dist=(Havers_dist/1000)) %>%
  ggplot(.,aes(x=Havers_dist,y=BC))+theme_classic()+
  geom_point(size=2,shape=21,fill="grey",color="black")+
  stat_smooth(method = "lm",show.legend = F,color="#4166F5",fill="lightblue",se = T)+
  theme(axis.ticks = element_blank(),
        axis.text = element_text(size=12),
        axis.title=element_text(size=12,face="bold"),aspect.ratio = 1,
        plot.caption = element_text(size=12,face="italic"))+
  #scale_color_manual(values=palette_lake)+
  labs(y="Taxa-content dissimilarity",color="Month",caption=var_BC_geo,
       x="Distance between lakes (km)")

####__ GHs BC ####

Fig_hellinger_GH <-
  MDS_GH_hellinger.df %>%
  dplyr::mutate(.,
                lake=
                  if_else(lake=="CHA","Champs",
                          if_else(lake=="VER","Verneuil",
                                  if_else(lake=="GDP","Grande-Paroisse",
                                          if_else(lake=="BOI","Bois",
                                                  if_else(lake=="CTL","Créteil",
                                                          if_else(lake=="CRJ1","Cergy-small",
                                                                  if_else(lake=="CRJ2","Cergy-large",
                                                                          if_else(lake=="VAI","Vaires","Jablines")))))))),
                lake=factor(lake,levels=rev(c("Verneuil","Champs",
                                              "Grande-Paroisse","Bois","Créteil","Cergy-small","Cergy-large",
                                              "Vaires","Jablines")))) %>%
  ggplot(.,aes(dim1,dim2,
               #shape=month_full,
               color = lake)) +
  geom_convexhull(aes(color = lake,fill=lake,
                      group =
                        factor(lake,levels=
                                 c("Verneuil","Champs",
                                   "Grande-Paroisse","Bois","Créteil","Cergy-small","Cergy-large",
                                   "Vaires","Jablines"))),
                  alpha=0.4,show.legend = F,size=0.4)+
  geom_point(size=2,
             aes(group =
                   factor(lake,levels=c("Verneuil","Champs",
                                        "Grande-Paroisse","Bois","Créteil","Cergy-small","Cergy-large",
                                        "Vaires","Jablines")))) + theme_bw() +
  # geom_label(Ko_centroids_hellinger.df,mapping=aes(dim1,dim2,label = lake, fill = lake),
  #            color = "white",size = 4,fontface = "bold",show.legend = F,inherit.aes = F)+
  theme(panel.grid = element_blank(),
        plot.title= element_text(size=11.5,hjust=0),
        #plot.title= element_text(size=14,hjust=0),
        axis.title = element_text(size=12,face = "bold"),
        axis.text = element_blank(),axis.ticks = element_blank(),
        legend.title = element_text(size=12,face = "bold",hjust=0),
        legend.text = element_text(size=10),legend.key.size = unit(0.5,"cm"),
        legend.position = "bottom",legend.direction = "vertical")+
  scale_color_manual(values=rev(palette_lake_3T))+
  scale_fill_manual(values=rev(palette_lake_3T))+
  #scale_shape_manual(values = Summer_2021)+
  # guides(color = guide_legend(override.aes=list(size=5)),
  #        shape = guide_legend(override.aes=list(size=5)))+
  scale_y_reverse()+
  labs(color="Lake",#shape="Month",
       title="GHs families community structure",
       x=paste0("Axis 1 [",round(MDS_GH.hellinger.mds$eig[1]*100/sum(MDS_GH.hellinger.mds$eig),1),"%]"),
       y=paste0("Axis 2 [",round(MDS_GH.hellinger.mds$eig[2]*100/sum(MDS_GH.hellinger.mds$eig),1),"%]"))

Fig_hellinger_GH<-Fig_hellinger_GH+
  guides(color = guide_legend(override.aes=list(shape=15,size=6),order = 1,
                              nrow=3))

GH.plsda <- MDS_GH.hellinger %>%
  as.matrix(.) %>% mixOmics::plsda(.,hellinger_factor$Status,scale = F)

summer_GH.plsda.df<-
  cbind(comp1=as.data.frame(GH.plsda$variates$X)$comp1,
        comp2=as.data.frame(GH.plsda$variates$X)$comp2) %>% as.data.frame(.) %>%
  dplyr::mutate(.,
                status=GH.plsda$Y,status=factor(status,levels =c("Hypereutrophic","Eutrophic",
                                                                        "Mesotrophic")))
# export the n first loading variable
n_VIP=10
Fig_summer_GH.plsda.loadings<-
  plotLoadings(GH.plsda,ndisplay = n_VIP,method = 'median',contrib = "max",
               title = paste0("PLS-DA \"Trophic status\" ",n_VIP,
                              " most disciminant GHs"),
               size.title = 1.2,legend.color = c("#2F6B9D","#E45000","#28A448"))
View(Fig_summer_GH.plsda.loadings)
list_5_GH_summer<-
  cbind(PLS_DA_5=rownames(as.data.frame(Fig_summer_GH.plsda.loadings)),
        importance=as.data.frame(Fig_summer_GH.plsda.loadings)$importance,
        status_contrib=as.data.frame(Fig_summer_GH.plsda.loadings)$GroupContrib) %>%
  as.data.frame()

View(list_5_GH_summer)

Fig_jaccard_GH <-
  MDS_GH.jaccard.df %>%
  dplyr::mutate(.,
                lake=
                  if_else(lake=="CHA","Champs",
                          if_else(lake=="VER","Verneuil",
                                  if_else(lake=="GDP","Grande-Paroisse",
                                          if_else(lake=="BOI","Bois",
                                                  if_else(lake=="CTL","Créteil",
                                                          if_else(lake=="CRJ1","Cergy-small",
                                                                  if_else(lake=="CRJ2","Cergy-large",
                                                                          if_else(lake=="VAI","Vaires","Jablines")))))))),
                lake=factor(lake,levels=rev(c("Verneuil","Champs",
                                              "Grande-Paroisse","Bois","Créteil","Cergy-small","Cergy-large",
                                              "Vaires","Jablines")))) %>%
  ggplot(.,aes(dim1,dim2,
               #shape=month_full,
               color = lake)) +
  geom_convexhull(aes(color = lake,fill=lake,
                      group =
                        factor(lake,levels=
                                 c("Verneuil","Champs",
                                   "Grande-Paroisse","Bois","Créteil","Cergy-small","Cergy-large",
                                   "Vaires","Jablines"))),
                  alpha=0.4,show.legend = F,size=0.4)+
  geom_point(size=2,
             aes(group =
                   factor(lake,levels=c("Verneuil","Champs",
                                        "Grande-Paroisse","Bois","Créteil","Cergy-small","Cergy-large",
                                        "Vaires","Jablines")))) + theme_bw() +
  # geom_label(Ko_centroids_hellinger.df,mapping=aes(dim1,dim2,label = lake, fill = lake),
  #            color = "white",size = 4,fontface = "bold",show.legend = F,inherit.aes = F)+
  theme(panel.grid = element_blank(),
        #plot.title= element_text(size=11.5,hjust=0),
        plot.title= element_text(size=20,hjust=0),
        axis.title = element_text(size=20,face = "bold"),
        axis.text = element_blank(),axis.ticks = element_blank(),
        legend.title = element_text(size=20,face = "bold",hjust=0),
        legend.text = element_text(size=15),legend.key.size = unit(0.5,"cm"),
        legend.position = "bottom",legend.direction = "vertical")+
  scale_color_manual(values=rev(palette_lake_3T))+
  scale_fill_manual(values=rev(palette_lake_3T))+
  #scale_shape_manual(values = Summer_2021)+
  # guides(color = guide_legend(override.aes=list(size=5)),
  #        shape = guide_legend(override.aes=list(size=5)))+
  scale_x_reverse()+
  scale_y_reverse()+
  labs(color="Lake",#shape="Month",
       title="Annotated GHs community content",
       x=paste0("Axis 1 [",round(MDS_GH.jaccard.mds$eig[1]*100/sum(MDS_GH.jaccard.mds$eig),1),"%]"),
       y=paste0("Axis 2 [",round(MDS_GH.jaccard.mds$eig[2]*100/sum(MDS_GH.jaccard.mds$eig),1),"%]"))

Fig_jaccard_GH<-Fig_jaccard_GH+
  guides(color = guide_legend(override.aes=list(shape=15,size=6),order = 1,
                              nrow=3))


####__ BGC BC ####

Fig_hellinger_BGC <-
  BGC_MDS_hellinger.df %>%
  dplyr::mutate(.,
                lake=
                  if_else(lake=="CHA","Champs-sur-Marne (CSM)",
                          if_else(lake=="VER","Verneuil-sur-Seine (VSS)",
                                  if_else(lake=="GDP","La Grande-Paroisse (LGP)",
                                          if_else(lake=="BOI","Bois-le-Roi (BLR)",
                                                  if_else(lake=="CTL","Créteil (CRE)",
                                                          if_else(lake=="CRJ1","Cergy - small lake (CER-S)",
                                                                  if_else(lake=="CRJ2","Cergy - large lake (CER-L)",
                                                                          if_else(lake=="VAI","Vaires-sur-Marne (VSM)","Jablines (JAB)")))))))),
                lake=factor(lake,levels=rev(c("Verneuil-sur-Seine (VSS)","Champs-sur-Marne (CSM)",
                                              "La Grande-Paroisse (LGP)","Bois-le-Roi (BLR)","Créteil (CRE)",
                                              "Cergy - small lake (CER-S)","Cergy - large lake (CER-L)",
                                              "Vaires-sur-Marne (VSM)","Jablines (JAB)")))) %>%
  ggplot(.,aes(dim1,dim2,
               color = lake)) +
  geom_convexhull(aes(color = lake,fill=lake,
                      group =
                        factor(lake,levels=
                                 c("Verneuil-sur-Seine (VSS)","Champs-sur-Marne (CSM)",
                                   "La Grande-Paroisse (LGP)","Bois-le-Roi (BLR)","Créteil (CRE)",
                                   "Cergy - small lake  (CER-S)","Cergy - large lake (CER-L)",
                                   "Vaires-sur-Marne (VSM)","Jablines (JAB)"))),
                  alpha=0.4,show.legend = F,size=0.4)+
  geom_point(size=2,show.legend = T,
             aes(group =
                   factor(lake,levels=c("Verneuil-sur-Seine (VSS)","Champs-sur-Marne (CSM)",
                                        "La Grande-Paroisse (LGP)","Bois-le-Roi (BLR)","Créteil (CRE)",
                                        "Cergy - small lake (CER-S)","Cergy - large lake (CER-L)",
                                        "Vaires-sur-Marne (VSM)","Jablines (JAB)")))) +
  # geom_segment(data = BGC_env.df,
  #              mapping=aes(x = 0, xend = Dim1*0.18, y = 0, yend = Dim2*0.18),
  #              color = "black",inherit.aes = F,arrow = arrow(length = unit(0.25,"cm"),type = "closed")) +
  # geom_text(data = BGC_env.df, mapping=aes(x = Dim1*0.19, y = Dim2*0.19, label = gene_metabolism),
  #            size = 3,fill="white",fontface = "bold",inherit.aes = F)+
  theme_bw() +
  theme(panel.grid = element_blank(),
        plot.title= element_blank(),
        axis.title = element_text(size=12,face = "bold"),
        axis.text = element_text(size=10),axis.ticks = element_blank(),
        legend.title = element_text(size=15,face = "bold",hjust=0),
        legend.text = element_text(size=12),
        legend.position = "bottom",legend.direction = "vertical")+
  scale_color_manual(values=rev(palette_lake_3T))+
  scale_fill_manual(values=rev(palette_lake_3T))+
  #scale_x_reverse()+
  guides(color=guide_legend(override.aes = list(shape=15,size=5),order = 1,nrow = 3),fill=F)+
  labs(color="Lake",
       x=paste0("Axis 1 [",round(MDS_BGC.hellinger.mds$eig[1]*100/sum(MDS_BGC.hellinger.mds$eig),1),"%]"),
       y=paste0("Axis 2 [",round(MDS_BGC.hellinger.mds$eig[2]*100/sum(MDS_BGC.hellinger.mds$eig),1),"%]"))

Fig_hellinger_BGC

####__ BGC short ####
BGC_short.boxplot<- BGC_short.df %>%
  dplyr::mutate(.,gene_metabolism=
                  if_else(gene_metabolism=="Carbon_fixation","C fixation",
                          if_else(gene_metabolism=="P_metabolism","P metabolism",
                                  if_else(gene_metabolism=="N_metabolism","N metabolism",
                                          if_else(gene_metabolism=="S_metabolism","S metabolism",gene_metabolism)))),
                gene_metabolism=factor(gene_metabolism,levels = c("C fixation","N metabolism","P metabolism","S metabolism","Iron-related"))) %>%
  ggplot(.,aes(Status,(CPM_sum),fill=gene_metabolism))+
  geom_boxplot(width=0.15,color="black",show.legend = F,outlier.size = 0.5)+theme_bw()+
  geom_point(show.legend = F,color="black",size=0.5)+
  theme(panel.grid = element_blank(),axis.ticks = element_blank(),
        axis.title.x = element_blank(),axis.title.y = element_blank(),
        axis.text.y = element_text(size=10),
        axis.text.x = element_blank())+
  scale_fill_manual(values=palette_BGC_short)+
  scale_y_continuous(limits = c(0,NA),name = NULL, sec.axis = sec_axis(~.))+
  guides(y = "none")+
  facet_grid2(gene_metabolism ~ Status,scales = "free",
              strip =
                strip_color_BGC_CP<- strip_themed(background_x = elem_list_rect(fill = c("#2F6B9D","#E45000","#28A448"),color = "black"),
                                                  text_x = elem_list_text(colour = "white",face = "bold",size=10),
                                                  background_y = elem_list_rect(fill = palette_BGC_short ,color = "black"),
                                                  text_y = elem_list_text(colour = "white",face = "bold",size=12)),switch = "y")


#### P####
P_short.df<- all_fun.df %>%
  subset(month %in% c("A","B","C","D")) %>%
  dplyr::group_by(sample) %>% dplyr::mutate(.,CPM_tot = sum(CPM)) %>%
  subset(KEGG_ko %in% BGC_short_list) %>%
  cbind(.,
        gene=BGC_short$Level_2[match(.$KEGG_ko,BGC_short$Kegg_ko)],
        gene_metabolism=BGC_short$Level_1[match(.$KEGG_ko,BGC_short$Kegg_ko)]) %>%
  dplyr::group_by(sample,gene_metabolism,KEGG_ko) %>% dplyr::mutate(.,CPM_rel = (CPM/CPM_tot)*100) %>%
  dplyr::group_by(sample,gene_metabolism,KEGG_ko) %>% dplyr::summarise(CPM_sum=sum(CPM_rel)) %>%
  cbind(.,
        month_full=sample_metadata$month_full[match(.$sample,sample_metadata$sample)],
        Status=sample_metadata$Status[match(.$sample,sample_metadata$sample)],
        lake=sample_metadata$lake[match(.$sample,sample_metadata$sample)],
        gene=BGC_short$Level_2[match(.$KEGG_ko,BGC_short$Kegg_ko)]) %>%
  dplyr::mutate(.,
                gene_metabolism=factor(gene_metabolism,
                                       levels = c("Carbon_fixation","N_metabolism",
                                                  "P_metabolism","S_metabolism",
                                                  "Methane_metabolism","Iron-related")),
                Status=if_else(lake %in% c("CHA","VER"),"Hypereutrophic",
                               if_else(lake %in% c("JAB","VAI"),"Mesotrophic",
                                       "Eutrophic")),
                Status=factor(Status,levels=c("Mesotrophic","Eutrophic","Hypereutrophic")),
                lake=
                  if_else(lake=="VER","Verneuil",
                          if_else(lake=="CHA","Champs",
                                  if_else(lake=="GDP","Grande-Paroisse",
                                          if_else(lake=="BOI","Bois",
                                                  if_else(lake=="CTL","Créteil",
                                                          if_else(lake=="CRJ1","Cergy-large",
                                                                  if_else(lake=="CRJ2","Cergy-small",
                                                                          if_else(lake=="VAI","Vaires","Jablines")))))))),
                lake=factor(lake,levels=c("Verneuil","Champs",
                                          "Grande-Paroisse","Bois","Créteil","Cergy-large","Cergy-small",
                                          "Vaires","Jablines"))) %>%
  subset(gene_metabolism == "P_metabolism")
  
P_short.boxplot<- P_short.df %>%
  ggplot(.,aes(Status,CPM_sum),fill="#3E842E")+
  geom_boxplot(width=0.15,color="black",show.legend = F,outlier.size = 0.5)+theme_bw()+
  geom_point(show.legend = F,color="black",size=0.5)+
  theme(panel.grid = element_blank(),axis.ticks = element_blank(),
        axis.title.x = element_blank(),axis.title.y = element_blank(),
        axis.text.y = element_text(size=10),
        axis.text.x = element_blank())+
  #scale_fill_manual(values=palette_BGC_short)+
  scale_y_continuous(limits = c(0,NA),name = NULL, sec.axis = sec_axis(~.))+
  guides(y = "none")+
  facet_grid2(gene~ Status,scales = "free",
              strip =
                strip_color_BGC_CP<- strip_themed(background_x = elem_list_rect(fill = c("#2F6B9D","#E45000","#28A448"),color = "black"),
                                                  text_x = elem_list_text(colour = "white",face = "bold",size=10),
                                                  background_y = elem_list_rect(fill ="#3E842E" ,color = "black"),
                                                  text_y = elem_list_text(colour = "white",face = "bold",size=12)),switch = "y")






palette_BGC_short<-c("#7BA8CF","#9B5C97","#3E842E","#C84D4C","black")

BGC_short.df %>% group_by(sample) %>%
  summarise(CPM=sum(CPM_sum)) %>%
  ungroup() %>% summarise(.,CPM=median(CPM))

BGC_short.df %>% group_by(sample) %>%
  summarise(CPM=sum(CPM_sum)) %>%
  ungroup() %>% summarise(.,SD=sd(CPM))

BGC_short.boxplot$data %>% subset(.,gene_metabolism =="C fixation" & lake !="Bois") %>%
  dplyr::group_by(gene_metabolism) %>%
  dplyr::summarise_at(vars("CPM_sum"),list(median=median, sd=sd))

BGC_short.boxplot$data %>% subset(.,gene_metabolism =="C fixation" & lake !="Bois") %>%
  dplyr::group_by(gene_metabolism) %>%
  dplyr::summarise_at(vars("CPM_sum"),list(median=median, sd=sd))

BGC_short.boxplot$data %>% subset(.,gene_metabolism =="P metabolism") %>%
  dplyr::group_by(Status) %>%
  dplyr::summarise_at(vars("CPM_sum"),list(median=median, sd=sd))

BGC_short.boxplot$data %>% subset(.,gene_metabolism =="C fixation" & lake =="Bois") %>%
  dplyr::group_by(gene_metabolism) %>%
  dplyr::summarise_at(vars("CPM_sum"),list(median=median, sd=sd))

BGC_short.df %>% pivot_wider(.,names_from=gene_metabolism,values_from=CPM_sum,values_fill = 0)


BGC_metab.mds$lake_month_full=paste0(BGC_metab.mds$lake,"_",BGC_metab.mds$month_full)
Fig_Chla_dim1$data$Chla

tmp<-Fig_Chla_dim1$data %>%
  dplyr::mutate(.,
                lake=
                  if_else(lake=="CHA","Champs",
                          if_else(lake=="VER","Verneuil",
                                  if_else(lake=="GDP","Grande-Paroisse",
                                          if_else(lake=="BOI","Bois",
                                                  if_else(lake=="CTL","Créteil",
                                                          if_else(lake=="CRJ1","Cergy-small",
                                                                  if_else(lake=="CRJ2","Cergy-large",
                                                                          if_else(lake=="VAI","Vaires","Jablines")))))))),
                lake=factor(lake,levels=rev(c("Verneuil","Champs",
                                              "Grande-Paroisse","Bois","Créteil","Cergy-small","Cergy-large",
                                              "Vaires","Jablines"))),
                lake_month_full=paste0(lake,"_",month_full))

BGC_Chla<-BGC_metab.mds %>% cbind(.,Chla=tmp$Chla[match(.$lake_month_full,tmp$lake_month_full)])

BGC_Chla$P_metabolism
lm_gradient_model <- lm(P_metabolism ~ Chla,data=BGC_Chla)
summary(lm_gradient_model)

####|####

####__ KOs vs ASVs dynamic ####

BC_Hellinger.df <- BC_Hellinger.df %>%
  dplyr::mutate(.,
                comp_ST=if_else(lake1 %in% c("CHA","VER"),"Hypereutrophic",
                                if_else(lake1 %in% c("JAB","VAI"),"Mesotrophic",
                                        "Eutrophic")),
                comp_ST=factor(comp_ST,levels=c("Hypereutrophic","Eutrophic","Mesotrophic")))

####____ KOs BC boxplot ####

# Stats KOs
box_Hellinger.stats<-BC_Hellinger.df %>%
  wilcox_test(hellinger_dist ~ comp_ST,p.adjust.method = "BH") %>%
  add_significance("p.adj")  %>% add_xy_position(x = "comp_ST", dodge = 0.8) %>%
  mutate(g1=gsub("_","",as.character(.$group1)),g2=gsub("_","",as.character(.$group2))) %>%
  unite("comp",g1:g2,sep=" - ",remove = FALSE)

letters_box_Hellinger<-tibble(comp_ST=unique(sort(BC_Hellinger.df$comp_ST)),
                              letters=cldList(p.adj ~ as.character(comp),box_Hellinger.stats,threshold  = 0.05)$Letter)

# Boxplot BC KOs
Fig_box_Hellinger<- BC_Hellinger.df %>% dplyr::group_by(comp_ST) %>%
  dplyr::mutate(.,mean_hellinger_dist=mean(hellinger_dist)) %>%
  cbind(.,letter=letters_box_Hellinger$letters[match(.$comp_ST,letters_box_Hellinger$comp_ST)]) %>%
  ggplot(.,aes(x=fct_reorder(comp_ST,mean_hellinger_dist),y=hellinger_dist,fill=comp_ST))+theme_bw()+
  geom_boxplot(show.legend = F,width=0.6,color="black")+
  geom_text(mapping=aes(y=1.01*max(hellinger_dist),
                        label=letter),vjust=0,size=4.5,color="black",fontface="bold")+
  theme(panel.grid = element_blank(),plot.title = element_blank(),
        axis.ticks = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size=10),axis.text.x = element_blank(),
        legend.title = element_text(size=14,face = "bold"),
        legend.text = element_text(size=12))+
  scale_y_continuous(limits = c(min(BC_Hellinger.df$hellinger_dist)/1.01,max(BC_Hellinger.df$hellinger_dist)*1.01),
                     breaks = c(0.06,0.08,0.1,0.12,0.14,0.16),
                     name = NULL, sec.axis = sec_axis(~.))+
  guides(y = "none")+labs(fill="Trophic Status")+
  scale_fill_manual(values = c("#28A448","#E45000","#2F6B9D"))+
  guides(fill=guide_legend(reverse=T,override.aes = list(shape=22,size=5)))

Fig_box_Hellinger$data %>%
  dplyr::group_by(comp_ST) %>%
  dplyr::summarise_at(vars("hellinger_dist"),list(median=median, sd=sd))


write.csv(Fig_box_BC$data,"Fig_box_BC.csv")


####____ ASVs BC boxplot ####

# Stats ASVs
box_BC.stats<- BC_Hellinger.df %>% wilcox_test(BC_dist ~ comp_ST,p.adjust.method = "BH") %>%
  add_significance("p.adj")  %>% add_xy_position(x = "comp_ST", dodge = 0.8) %>%
  mutate(g1=gsub("_","",as.character(.$group1)),g2=gsub("_","",as.character(.$group2))) %>%
  unite("comp",g1:g2,sep=" - ",remove = FALSE)

letters_box_BC<-tibble(comp_ST=unique(sort(BC_Hellinger.df$comp_ST)),
                       letters=cldList(p.adj ~ as.character(comp),box_BC.stats,threshold  = 0.05)$Letter)
# Boxplot BC ASVs
Fig_box_BC<- BC_Hellinger.df %>% dplyr::group_by(comp_ST) %>%
  dplyr::mutate(.,mean_BC_dist=mean(BC_dist)) %>%
  cbind(.,letter=letters_box_BC$letters[match(.$comp_ST,letters_box_BC$comp_ST)]) %>%
  ggplot(.,aes(x=comp_ST,y=BC_dist,fill=comp_ST))+theme_bw()+
  geom_boxplot(show.legend = T,width=0.6,color="black")+
  geom_text(mapping=aes(y=1.01*max(BC_dist),
                        label=letter),vjust=0,size=4.5,color="black",fontface="bold")+
  theme(panel.grid = element_blank(),plot.title = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        axis.text.x = element_text(size=10),
        axis.text.y = element_blank(),
        legend.title = element_text(size=14,face = "bold"),
        legend.text = element_text(size=12))+
  scale_y_continuous(limits = c(min(BC_Hellinger.df$BC_dist)/1.01,max(BC_Hellinger.df$BC_dist)*1.01),
                     breaks = c(0.30,0.4,0.5,0.6,0.7,0.8,0.9),position = "right")+
  scale_fill_manual(values=c("#28A448","#E45000","#2F6B9D"))+
  guides(fill = guide_legend(ncol=1))+
  labs(fill="Trophic Status")+coord_flip()

Fig_box_BC

Fig_box_BC$data %>%
  dplyr::group_by(comp_ST) %>%
  dplyr::summarise_at(vars("BC_dist"),list(median=median, sd=sd))

####____ KOs VS ASVs BC ####

spearman_BC<-cor(BC_Hellinger.df$BC_dist,BC_Hellinger.df$hellinger_dist,method = "spearman")

spearman_BC<-cor.test(BC_Hellinger.df$BC_dist,BC_Hellinger.df$hellinger_dist,method = "spearman")
cor.test

lm_BC_model <- lm(hellinger_dist~BC_dist*comp_ST,data=BC_Hellinger.df)

var_BC=paste0("Spearman ",round(spearman_BC,2),", R2 ",round(summary(lm_BC_model)[["adj.r.squared"]],2),", p <0.01")

BC_fun_centroids.df <- BC_Hellinger.df %>% dplyr::group_by(lake1) %>% dplyr::summarise_at(vars(3,14),median) %>%
  cbind(.,comp_ST=BC_Hellinger.df$comp_ST[match(.$lake1,BC_Hellinger.df$lake1)]) %>%
  dplyr::mutate(.,
                lake=
                  if_else(lake1=="CHA","CM",
                          if_else(lake1=="VER","VS",
                                  if_else(lake1=="GDP","GP",
                                          if_else(lake1=="BOI","BR",
                                                  if_else(lake1=="CTL","CR",
                                                          if_else(lake1=="CRJ1","CP-S",
                                                                  if_else(lake1=="CRJ2","CP-L",
                                                                          if_else(lake1=="VAI","VM","JA")))))))),
                lake=factor(lake,levels=rev(c("VS","CM",
                                              "GP","BR","CR","CP-S","CP-L",
                                              "VM","JA"))))

Fig_BC_Hellinger<- BC_Hellinger.df %>%
  dplyr::mutate(.,
                lake=
                  if_else(lake1=="CHA","CSM",
                          if_else(lake1=="VER","VSS",
                                  if_else(lake1=="GDP","LGP",
                                          if_else(lake1=="BOI","BLR",
                                                  if_else(lake1=="CTL","CRE",
                                                          if_else(lake1=="CRJ1","CER-S",
                                                                  if_else(lake1=="CRJ2","CER-L",
                                                                          if_else(lake1=="VAI","VSM","JAB")))))))),
                lake=factor(lake,levels=rev(c("VSS","CSM",
                                                        "LGP","BLR","CRE",
                                                        "CER-S","CER-L",
                                                        "VSM","JAB")))) %>%
  ggplot(.,aes(BC_dist,hellinger_dist,color=lake,fill=lake))+theme_bw()+
  geom_point(size=4,shape=16,show.legend = T )+
  annotate(geom="text", x=0.35, y=0.16,color="black",size=5,hjust=0,label="paste(italic(p), \"<0.01 \", italic(R) ^ 2, \" 0.65\")", parse = TRUE)+
  theme(
    aspect.ratio = 1,
    panel.grid = element_blank(),
    axis.title = element_text(size=14,face = "bold"),
    axis.text = element_text(size=12),
    axis.ticks = element_blank(),plot.caption = element_blank(),
    legend.title = element_text(size=15,face = "bold"),
    legend.text = element_text(size=12),
    legend.position = "bottom",legend.direction = "vertical")+
  scale_color_manual(values=rev(palette_lake_3T),
                     labels = rev(c("Verneuil-sur-Seine (VSS)","Champs-sur-Marne (CSM)",
                                "La Grande-Paroisse (LGP)","Bois-le-Roi (BLR)","Créteil (CRE)",
                                "Cergy - small lake (CER-S)","Cergy - large lake (CER-L)",
                                "Vaires-sur-Marne (VSM)","Jablines (JAB)")))+
  guides(color=guide_legend(reverse = F,nrow=3,override.aes=list(shape=15,size=5,label = "")),fill=FALSE)+
  scale_y_continuous(limits = c(min(BC_Hellinger.df$hellinger_dist)/1.01,0.16),
                     breaks = c(0.06,0.08,0.1,0.12,0.14,0.16))+
  scale_x_continuous(limits = c(min(BC_Hellinger.df$BC_dist)/1.01,max(BC_Hellinger.df$BC_dist)*1.01),
    breaks = c(0.30,0.4,0.5,0.6,0.7,0.8,0.9))+
  labs(x="Taxa-content dissimilarity",
       y="Gene-content dissimilarity",color="Lake")
 

Fig_BC_Hellinger

ggsave("/Users/pierre/Desktop/ms1/Figures/to_inkscape/F3.pdf", height = 9.8, width = 13.4, units = "in",dpi = "retina")


View(Fig_BC_Hellinger$data )
Fig_BC_Hellinger$data %>% dplyr::mutate(.,BOI=if_else(lake1=="BOI","BOI","Others")) %>%
  dplyr::group_by(BOI) %>%
  dplyr::summarise_at(vars("hellinger_dist"),list(min=min,max=max))

design_fun_tax="
AAAAB
AAAAB
AAAAB
AAAAB
AAAAB
CCCC#"

Panel_redun_BC<-Fig_BC_Hellinger+Fig_box_Hellinger+Fig_box_BC+
  plot_layout(design = design_fun_tax)+
  plot_annotation(tag_levels = 'A') & theme(plot.tag=element_text(size=20,face="bold",color="black"),legend.position = "bottom",legend.direction = "vertical")



####____ Overall disp ####
View(BC_Hellinger.df)
var_summer.df<- read_delim("Fig_box_var.csv",delim = ";")

var_summer.stats<- var_summer.df %>%
  wilcox_test(hellinger_dist ~ comp_ST_features,p.adjust.method = "BH") %>%
  add_significance("p.adj")  %>% add_xy_position(x = "comp_ST_features", dodge = 0.8) %>%
  mutate(g1=gsub("_","",as.character(.$group1)),g2=gsub("_","",as.character(.$group2))) %>%
  unite("comp",g1:g2,sep=" - ",remove = FALSE)

letters_var_summer<-tibble(comp_ST_features=unique(sort(var_summer.df$comp_ST_features)),
                            letters=cldList(p.adj ~ as.character(comp),var_summer.stats,threshold  = 0.05)$Letter)


Fig_var_summer <- var_summer.df %>%
  dplyr::mutate(.,comp_ST=factor(comp_ST,levels=rev(c("Hypereutrophic","Eutrophic","Mesotrophic"))),
                comp_ST_features=factor(comp_ST_features,levels=c("Mesotrophic_KOs","Mesotrophic_ASVs",
                                                                  "Eutrophic_KOs","Eutrophic_ASVs",
                                                                  "Hypereutrophic_KOs","Hypereutrophic_ASVs"))) %>%
  cbind(.,
        letter=letters_var_summer$letters[match(.$comp_ST_features,letters_var_summer$comp_ST_features)]) %>%
  ggplot(.,aes(comp_ST_features,hellinger_dist,fill=features))+
  geom_boxplot(width=0.6,color="black",show.legend = F)+theme_bw()+
  geom_point(size=1,color="black",alpha=0.5,show.legend = F)+
  geom_text(mapping=aes(y=1.02*max(hellinger_dist),
                        label=letter),vjust=0,size=4.5,color="black",fontface="bold")+
  theme(panel.grid = element_blank(),
        axis.ticks = element_blank(),axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size=10,hjust=0.5),axis.text.y = element_text(size=12),
        legend.title = element_text(size=14,face = "bold"),legend.text = element_text(size=12))+
  scale_fill_manual(values = c("lightgrey","white"))+
  scale_y_continuous(limits =c(0,1.06*max(var_summer.df$hellinger_dist)),position = "right",
                     labels = scales::label_number(accuracy = 0.1))+
  scale_x_discrete(labels=c('KOs', 'ASVs',
                            'KOs', 'ASVs',
                            'KOs', 'ASVs'))+
  labs(fill="Features")+
  facet_wrap2(~ comp_ST,scale="free_x",nrow = 1,strip.position = "top",
              strip = strip_color_lake_type_afem<- strip_themed(
                background_x = elem_list_rect(fill = c("#2F6B9D","#E45000","#28A448"),
                                              color="black"),
                text_x = elem_list_text(colour = "white",
                                        face = "bold",size=10)))
Fig_var_summer

Fig_var_summer$data %>%
  dplyr::group_by(comp_ST_features) %>%
  dplyr::summarise_at(vars("hellinger_dist"),list(mean=mean, sd=sd))

Fig_var_summer$data %>% dplyr::group_by(comp_ST_features) %>%
  dplyr::summarise_at(vars("hellinger_dist"),list(mean=mean, sd=sd))

testlm<-lmer(hellinger_dist ~ comp_ST + (1|lake1),data=subset(Fig_var_summer$data,features=="KOs"),REML = F)
summary(testlm)

Fig_disp_summer <- disp_summer.df %>%
  dplyr::group_by(Status_features) %>%
  dplyr::mutate(.,median_distance=median(distance)) %>%
  cbind(.,letter=letters_disp_summer$letters[match(.$Status_features,letters_disp_summer$Status_features)],) %>%
  ggplot(.,aes(Status_features,distance,fill=Features))+
  geom_boxplot(width=0.6,color="black",show.legend = F)+theme_bw()+
  geom_point(size=1,color="black",alpha=0.5,show.legend = F)+
  geom_text(mapping=aes(y=1.02*max(distance),
                        label=letter),vjust=0,size=4.5,color="black",fontface="bold")+
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(size=10,hjust=0.5),axis.text.y = element_text(size=12),
        axis.ticks = element_blank(),axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.title = element_blank(),legend.text = element_text(size=12))+
  scale_fill_manual(values = c("white","lightgrey"), breaks = c("KOs","ASVs"))+
  guides(fill=guide_legend(nrow=1,override.aes = list(size=5)))+
  scale_y_continuous(limits =c(0,1.06*max(disp_summer.df$distance)),position = "right",
                     labels = scales::label_number(accuracy = 0.1))+
  scale_x_discrete(labels=c('KOs', 'ASVs',
                            'KOs', 'ASVs',
                            'KOs', 'ASVs'))+
  facet_wrap2(~ Status,scale="free_x",nrow = 1,strip.position = "top",
              strip = strip_color_lake_type_afem<- strip_themed(
                background_x = elem_list_rect(fill = c("#2F6B9D","#E45000","#28A448"),
                                              color="black"),
                text_x = elem_list_text(colour = "white",
                                        face = "bold",size=10)))

Fig_disp_summer$data %>%
  dplyr::group_by(Status_features) %>%
  dplyr::summarise_at(vars("distance"),list(mean=mean, sd=sd))

Fig_disp_summer$data %>% dplyr::group_by(Status_features) %>%
  dplyr::summarise_at(vars("distance"),list(mean=mean, sd=sd))

testlm<-lmer(distance~ Status + (1|lake),data=subset(Fig_disp_summer$data,Features=="ASVs"),REML = F)
summary(testlm)

####|####

#### Overall panel ####

####______ F1 ####
design_afem_F1<-"
AC
AC
BC"

Panel_lake<- Fig_lake_median_Chla_hor+Fig_lake_median_Chla_hor+Fig_Phyto_domain+
  plot_layout(guides = 'collect',design = design_afem_F1)+plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size=15,face="bold"),
        legend.key.height= unit(0.5, 'cm'),
        legend.key.width= unit(0.5, 'cm'),
        legend.spacing.y = unit(0.5, 'cm'),
        legend.spacing.x = unit(0.2, 'cm'),
        legend.text = element_text(size=12),
        legend.title = element_text(size=15,face="bold"),
        legend.position = "bottom",legend.direction = "vertical")
Panel_lake

ggsave("/Users/pierre/Desktop/ms1/Figures/to_inkscape/F1.pdf", height = 9.8, width = 13.4, units = "in",dpi = "retina")

####______ F2 ####
design_afem_F2="
AACC
BBCC
DDEE
DDEE"

Fig_shared+Fig_box_ASV_KO_turnover+Barplot_AFEM+Fig_hellinger_ko+Fig_BC_summer+
  plot_layout(design=design_afem_F2,guides = 'collect')+
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size=15,face="bold"),
        legend.position = "bottom",legend.direction = "vertical")
  

ggsave("/Users/pierre/Desktop/ms1/Figures/to_inkscape/F2.pdf", height = 9.8, width = 13.4, units = "in",dpi = "retina")

####______ F3 ####
design_F3="
AACCC
BBCCC"

Fig_hellinger_BGC+Fig_hellinger_GH+Panel_redun_BC+
  plot_layout(design=design_F3,guides = 'collect')+
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size=15,face="bold"),
        legend.position = "bottom",legend.direction = "vertical")


Fig_BC_Hellinger+Fig_hellinger_BGC+
  plot_layout(guides = 'collect')+
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size=15,face="bold"),
        legend.position = "bottom",legend.direction = "vertical")


####______ F4 ####
design_F4="
AABB
AACC"

Fig_hellinger_BGC+Fig_simper_H_M+Fig_simper_BOI+
  plot_layout(guides = 'collect',design = design_F4)+
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size=15,face="bold"),
        legend.position = "bottom",legend.direction = "vertical")
  
ggsave("/Users/pierre/Desktop/ms1/Figures/to_inkscape/F4.pdf", height = 9.8, width = 13.4, units = "in",dpi = "retina")

####______ FS2 ####
design_S2="
AB
AB
AC"

Fig_Cyano+Fig_mds_param+Fig_Chla_dim1+
  plot_layout(guides = 'collect',design = design_S2)+
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size=15,face="bold"),
        legend.position = "bottom",legend.direction = "vertical")

ggsave("/Users/pierre/Desktop/ms1/Figures/to_inkscape/FS2.pdf", height = 9.8, width = 13.4, units = "in",dpi = "retina")

####______ FS3 ####

Fig_Hellinger_Geo+Fig_BC_Geo

ggsave("/Users/pierre/Desktop/ms1/Figures/to_inkscape/FS3.pdf", height = 9.8, width = 13.4, units = "in",dpi = "retina")

####______ FS4 ####
design_S4="
AB"

Fig_disp_summer+Fig_var_summer+
  plot_layout(design=design_S4,guides = 'collect')+
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size=15,face="bold"),
        legend.position = "bottom",legend.direction = "vertical")

####______ FS5 ####
design_S5="
A
B
C"

Fig_box_ASV_D_lake+Fig_box_ASV_S_lake+Fig_box_ASV_H_lake+
  plot_layout(guides = 'collect',design = design_S5) +
  plot_annotation(tag_levels = 'A') &
  theme(plot.tag = element_text(size=15,face="bold"),
        legend.position = "bottom",legend.direction = "vertical")



####______ FS6 ####

design_S6="
ADD
BDD
CDD"

Fig_simper_H_M_taxa+Fig_simper_E_M_taxa+Fig_simper_E_H_taxa+
  Fig_simper_ASV_gg+
  plot_layout(design=design_S6)+
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size=15,face="bold"))

  
  
  