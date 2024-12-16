#### Import packages ####

library(BiocManager)
library(phyloseq) 
library(mixOmics)
library(qiime2R) 
library(vegan) 
library(tidyverse)
library(tidytext)
library(ggh4x)
library(patchwork)
library(rstatix)
library(rcompanion)
library(reshape2)
library(codyn)
library(geosphere)
library(ggConvexHull)
library(foreign)
library(pairwiseAdonis)
library(scales) 

#for qiime2R :
#install.packages("devtools") #if you don't have devtools installed
#devtools::install_github("jbisanz/qiime2R")

#pour rcompanion :
#packageurl <- "https://cran.r-project.org/src/contrib/Archive/rcompanion/rcompanion_2.3.26.tar.gz"
#install.packages(packageurl, repos=NULL, type="source")
#need ‘DescTools’, ‘multcompView’, ‘EMT’, ‘lmtest’, ‘nortest’

####|####

#### Import data ####

all_rarW<-qza_to_phyloseq(features="rar_id_w_t_filtered_table.qza",
                          taxonomy="rar_f_w_taxonomy_C2L_A_F.qza",
                          metadata="rar_w_f_metadata_C2L_A_F_16S.txt",
                          tree="rar_w_f_rooted_tree.qza")

all_rarW@sam_data<-
  read_delim("ps_rar_w_f_samdata_ok.csv",";", escape_double = FALSE,
             trim_ws = TRUE,show_col_types = FALSE) %>%
  remove_rownames %>%
  column_to_rownames(var="Samples") %>%
  dplyr::mutate(.,
                samples=rownames(.),
                lake_col=paste0(lake,"_",column),
                lake_type=if_else(lake %in% c("CHA","VER"),"Hypereutrophic",
                                 if_else(lake %in% c("JAB","VAI"),"Mesotrophic",
                                         "Eutrophic")),
                lake_type=factor(lake_type,levels=c("Hypereutrophic","Eutrophic","Mesotrophic")),
                month_lake_type=paste0(month,"_",lake_type),
                lake_order=factor(lake,
                                  levels=c("CHA","VER","GDP","BOI","CTL","CRJ1","CRJ2","VAI","JAB")),
                month_letter=factor(month_letter,
                                    levels=c("Jun","Jul","Aug",
                                             "Sep","Oct","Nov")),
                month_full=factor(month_full,
                                  levels=c("June","July","August",
                                           "September","October","November")),
                Distance=
                  if_else(lake %in% 
                            c("CRJ1","CRJ2","VER"),"NW",
                          if_else(lake %in% c("CTL","GDP","BOI"),"S",
                                  "NE")),
                .after = 8
  ) %>% sample_data()

####__ summer subset ####

A_D_rarW<- all_rarW %>%
  subset_samples(samples !="VAI-B-W1-ADN1") %>%
  subset_samples(month %in% c("A","B","C","D"))

A_D_rarW_freq<- transform_sample_counts(A_D_rarW, function(x) 100 * x/sum(x)) 
View(A_D_rarW_freq@sam_data)


####|####

#### Palettes, strip & design ####

# To color the horizontal strips (=site name) according to the 3 levels of eutrophication 

palette_Charlotte<-c("#B3CDE3","#CCEBC5","#DECBE4","#FED9A6","#FFFFCC","#E5D8BD","#FDDAEC","#F2F2F2","#B3E2CD","#FDCDAC",
                     "#CBD5E8","#F4CAE4","#E6F5C9","#FFF2AE","#F1E2CC","#88F393","#2E6D74","#DCBCF3","#EFE0DC","#F18B6C",
                     "#B52239","#F5D2DB","#B7F7E6","#517BE6","#E0DBC1","#E19B8E","#BDCEC4","#99F8F3","#E45632","#D22A82",
                     "#70F4BB","#5D9BEF","#F391CC","#EEF4A3","#CAAB7D","#8AF3D4","#EAF7B9","#4263B2","#E0B333","#2CB550",
                     "#DA56EF","#EC97B3","#E0789F","#8CB8E2","#E89E34","#8D37A3","#71992A","#D8C5B2","#E75FA6","#EBE49A",
                     "#ADAC5A","#BDF6D3","#3BDCA8","#47D1C8","#F8EBD4","#F4F6D8","#B6A7BD","#D6F0F2","#3EA34E","#B8E2EC",
                     "#BE6AE7","#8B7664","#7161EC","#CACCA8","#F138F4","#8EAC73","#A899F0")


palette_Family<-c("#A76D42","#A5AEB9","#BFD190","#F4CAD4","#F170C3","#F09094","#568F72","#EDD5CF","#AEE5EF","#DDCFE9",
                  "#D2A23F","#DD699C","#79C995","#C4D7EE","#38AA4F","#72BFB6","#949EEF","#83ED8E","#D3F2D1","#5DBBEE",
                  "#EB8D43","#D5BDAC","#B1F2AB","#B8F47F","#9FCEF9","#E3F594","#EFC538","#EAF6B8","#92CE3F","#8064C1")

palette_lake_type<-c("#b2dcc0","#ffd9b2","#c9daf8")
palette_lake_type<-c("#2F6B9D","#E44600","#28A448")


pie(rep(1,length(palette_lake_type)), col=palette_lake_type)

palette_site<-c("#b2dcc0","#b2dcc0","#b2dcc0",
                "#ffd9b2","#ffd9b2","#ffd9b2",
                "#c9daf8","#c9daf8","#c9daf8")

palette_lake<-c("#2E6D74","#2CB550","#71992A",
                "#b19774","#D8C5B2",
                "#E19B8E","#C190A4",
                "#377EB8","#8CB8E2")
pie(rep(1,length(palette_lake)), col=palette_lake)

palette_lake<-c("#2E6D74","#71992A","#2CB550",
                "#E19B8E","#D8C5B2","#b19774",
                "#C190A4","#377EB8","#8CB8E2")

strip_color_lake<- strip_themed(
  background_x = elem_list_rect(fill = palette_lake,
                                color = "black"),
  text_x = elem_list_text(colour = "white",
                          face = "bold",size=10))

strip_color_lake_type<- strip_themed(
  background_x = elem_list_rect(fill = palette_lake_type,
                                color="black"),
  text_x = elem_list_text(colour = "#3b2c64",
                          face = "bold",size=10))

palette_lake_type

Fig_beta_rarW_design<- "
1122
1122
3344
3344"

####|####

#### Shared ASVs  ####

shared.df<-
  A_D_rarW %>% transform_sample_counts(., function(x) 100 * x/sum(x)) %>%
  psmelt(.) %>% .[,c(1:16,18,19,34:40)] %>%
  dplyr::group_by(OTU,lake_month) %>% dplyr::mutate(.,lake_rel_abund=median(Abundance),lake_rel_sd=sd(Abundance))
#View(shared_lake.df)

####__ Hypereutrophic ####
## VER
shared_lake.df <- shared.df %>% subset(.,lake=="VER")
lake_A=unique((subset(shared_lake.df,month=="A" & Abundance >0))$OTU)
lake_B=unique((subset(shared_lake.df,month=="B" & Abundance >0))$OTU)
lake_C=unique((subset(shared_lake.df,month=="C" & Abundance >0))$OTU)
lake_D=unique((subset(shared_lake.df,month=="D" & Abundance >0))$OTU)

Shared_VER<-as.data.frame(Reduce(intersect, list(lake_A,lake_B,lake_C,lake_D))) %>%
  rename("ASV"="Reduce(intersect, list(lake_A, lake_B, lake_C, lake_D))")

all_ASV<-length(unique((subset(shared_lake.df, Abundance >0))$OTU))

Shared_VER <- shared_lake.df %>% subset(.,OTU %in% as.vector(Shared_VER$ASV)) %>% .[,c(1,4,9,14,17,20,21,23,24,26,27)] %>%
  dplyr::group_by(OTU,month_full,Phylum) %>%
  dplyr::summarize(mean=mean(lake_rel_abund),sd=mean(lake_rel_sd)) %>%
  dplyr::group_by(month_full) %>%
  dplyr::mutate(.,sum_mean=sum(mean),sum_sd=sum(sd),nb_ASV=length(unique(.$OTU))) %>%
  ungroup() %>% dplyr::mutate(.,mean_lake=mean(sum_mean),mean_sd=mean(sum_sd),lake="VER",status="Hypereutrophic")

shared_ASV.df<-
  data.frame(lake=unique(Shared_VER$lake),status=unique(Shared_VER$status),
             nb_ASV=unique(Shared_VER$nb_ASV),
             rel.abund=unique(Shared_VER$mean_lake),
             rel.sd=unique(Shared_VER$mean_sd)) %>%
  dplyr::mutate(.,per_lake=(nb_ASV/all_ASV)*100)
unique(Shared_VER$Phylum)

## CHA
shared_lake.df <- shared.df %>% subset(.,lake=="CHA")
lake_A=unique((subset(shared_lake.df,month=="A" & Abundance >0))$OTU)
lake_B=unique((subset(shared_lake.df,month=="B" & Abundance >0))$OTU)
lake_C=unique((subset(shared_lake.df,month=="C" & Abundance >0))$OTU)
lake_D=unique((subset(shared_lake.df,month=="D" & Abundance >0))$OTU)

Shared_CHA<-as.data.frame(Reduce(intersect, list(lake_A,lake_B,lake_C,lake_D))) %>%
  rename("ASV"="Reduce(intersect, list(lake_A, lake_B, lake_C, lake_D))")

Shared_CHA <- shared_lake.df %>% subset(.,OTU %in% as.vector(Shared_CHA$ASV)) %>% .[,c(1,4,9,14,17,20,21,23,24,26,27)] %>%
  dplyr::group_by(OTU,month_full,Phylum) %>%
  dplyr::summarize(mean=mean(lake_rel_abund),sd=mean(lake_rel_sd)) %>%
  dplyr::group_by(month_full) %>%
  dplyr::mutate(.,sum_mean=sum(mean),sum_sd=sum(sd),nb_ASV=length(unique(.$OTU))) %>%
  ungroup() %>% dplyr::mutate(.,mean_lake=mean(sum_mean),mean_sd=mean(sum_sd),lake="CHA",status="Hypereutrophic")

all_ASV<- length(unique((subset(shared_lake.df, Abundance >0))$OTU))

unique(Shared_CHA$Phylum)
shared_ASV.df[nrow(shared_ASV.df) + 1,] = c(unique(Shared_CHA$lake),
                                            unique(Shared_CHA$status),
                                            unique(Shared_CHA$nb_ASV),
                                            unique(Shared_CHA$mean_lake),
                                            unique(Shared_CHA$mean_sd),
                                            (unique(Shared_CHA$nb_ASV)/all_ASV)*100)

####__ Eutrophic ####
## GDP
shared_lake.df <- shared.df %>% subset(.,lake=="GDP")
lake_A=unique((subset(shared_lake.df,month=="A" & Abundance >0))$OTU)
lake_B=unique((subset(shared_lake.df,month=="B" & Abundance >0))$OTU)
lake_C=unique((subset(shared_lake.df,month=="C" & Abundance >0))$OTU)
lake_D=unique((subset(shared_lake.df,month=="D" & Abundance >0))$OTU)

Shared_GDP<-as.data.frame(Reduce(intersect, list(lake_A,lake_B,lake_C,lake_D))) %>%
  rename("ASV"="Reduce(intersect, list(lake_A, lake_B, lake_C, lake_D))")

Shared_GDP <- shared_lake.df %>% subset(.,OTU %in% as.vector(Shared_GDP$ASV)) %>% .[,c(1,4,9,14,17,20,21,23,24,26,27)] %>%
  dplyr::group_by(OTU,month_full,Phylum) %>%
  dplyr::summarize(mean=mean(lake_rel_abund),sd=mean(lake_rel_sd)) %>%
  dplyr::group_by(month_full) %>%
  dplyr::mutate(.,sum_mean=sum(mean),sum_sd=sum(sd),nb_ASV=length(unique(.$OTU))) %>%
  ungroup() %>% dplyr::mutate(.,mean_lake=mean(sum_mean),mean_sd=mean(sum_sd),lake="GDP",status="Eutrophic")

unique(Shared_GDP$Phylum)
all_ASV<- length(unique((subset(shared_lake.df, Abundance >0))$OTU))
shared_ASV.df[nrow(shared_ASV.df) + 1,] = c(unique(Shared_GDP$lake),
                                            unique(Shared_GDP$status),
                                            unique(Shared_GDP$nb_ASV),
                                            unique(Shared_GDP$mean_lake),
                                            unique(Shared_GDP$mean_sd),
                                            (unique(Shared_GDP$nb_ASV)/all_ASV)*100)
## BOI
shared_lake.df <- shared.df %>% subset(.,lake=="BOI")
lake_A=unique((subset(shared_lake.df,month=="A" & Abundance >0))$OTU)
lake_B=unique((subset(shared_lake.df,month=="B" & Abundance >0))$OTU)
lake_C=unique((subset(shared_lake.df,month=="C" & Abundance >0))$OTU)
lake_D=unique((subset(shared_lake.df,month=="D" & Abundance >0))$OTU)

Shared_BOI<-as.data.frame(Reduce(intersect, list(lake_A,lake_B,lake_C,lake_D))) %>%
  rename("ASV"="Reduce(intersect, list(lake_A, lake_B, lake_C, lake_D))")

Shared_BOI <- shared_lake.df %>% subset(.,OTU %in% as.vector(Shared_BOI$ASV)) %>% .[,c(1,4,9,14,17,20,21,23,24,26,27)] %>%
  dplyr::group_by(OTU,month_full,Phylum) %>%
  dplyr::summarize(mean=mean(lake_rel_abund),sd=mean(lake_rel_sd)) %>%
  dplyr::group_by(month_full) %>%
  dplyr::mutate(.,sum_mean=sum(mean),sum_sd=sum(sd),nb_ASV=length(unique(.$OTU))) %>%
  ungroup() %>% dplyr::mutate(.,mean_lake=mean(sum_mean),mean_sd=mean(sum_sd),lake="BOI",status="Eutrophic")

unique(Shared_BOI$Phylum)
all_ASV<- length(unique((subset(shared_lake.df, Abundance >0))$OTU))

shared_ASV.df[nrow(shared_ASV.df) + 1,] = c(unique(Shared_BOI$lake),
                                            unique(Shared_BOI$status),
                                            unique(Shared_BOI$nb_ASV),
                                            unique(Shared_BOI$mean_lake),
                                            unique(Shared_BOI$mean_sd),
                                            (unique(Shared_BOI$nb_ASV)/all_ASV)*100)
## CTL
shared_lake.df <- shared.df %>% subset(.,lake=="CTL")
lake_A=unique((subset(shared_lake.df,month=="A" & Abundance >0))$OTU)
lake_B=unique((subset(shared_lake.df,month=="B" & Abundance >0))$OTU)
lake_C=unique((subset(shared_lake.df,month=="C" & Abundance >0))$OTU)
lake_D=unique((subset(shared_lake.df,month=="D" & Abundance >0))$OTU)

Shared_CTL<-as.data.frame(Reduce(intersect, list(lake_A,lake_B,lake_C,lake_D))) %>%
  rename("ASV"="Reduce(intersect, list(lake_A, lake_B, lake_C, lake_D))")

Shared_CTL <- shared_lake.df %>% subset(.,OTU %in% as.vector(Shared_CTL$ASV)) %>% .[,c(1,4,9,14,17,20,21,23,24,26,27)] %>%
  dplyr::group_by(OTU,month_full,Phylum) %>%
  dplyr::summarize(mean=mean(lake_rel_abund),sd=mean(lake_rel_sd)) %>%
  dplyr::group_by(month_full) %>%
  dplyr::mutate(.,sum_mean=sum(mean),sum_sd=sum(sd),nb_ASV=length(unique(.$OTU))) %>%
  ungroup() %>% dplyr::mutate(.,mean_lake=mean(sum_mean),mean_sd=mean(sum_sd),lake="CTL",status="Eutrophic")

all_ASV<- length(unique((subset(shared_lake.df, Abundance >0))$OTU))

shared_ASV.df[nrow(shared_ASV.df) + 1,] = c(unique(Shared_CTL$lake),
                                            unique(Shared_CTL$status),
                                            unique(Shared_CTL$nb_ASV),
                                            unique(Shared_CTL$mean_lake),
                                            unique(Shared_CTL$mean_sd),
                                            (unique(Shared_CTL$nb_ASV)/all_ASV)*100)
## CRJ1
shared_lake.df <- shared.df %>% subset(.,lake=="CRJ1")
lake_A=unique((subset(shared_lake.df,month=="A" & Abundance >0))$OTU)
lake_B=unique((subset(shared_lake.df,month=="B" & Abundance >0))$OTU)
lake_C=unique((subset(shared_lake.df,month=="C" & Abundance >0))$OTU)
lake_D=unique((subset(shared_lake.df,month=="D" & Abundance >0))$OTU)

Shared_CRJ1<-as.data.frame(Reduce(intersect, list(lake_A,lake_B,lake_C,lake_D))) %>%
  rename("ASV"="Reduce(intersect, list(lake_A, lake_B, lake_C, lake_D))")

Shared_CRJ1 <- shared_lake.df %>% subset(.,OTU %in% as.vector(Shared_CRJ1$ASV)) %>% .[,c(1,4,9,14,17,20,21,23,24,26,27)] %>%
  dplyr::group_by(OTU,month_full,Phylum) %>%
  dplyr::summarize(mean=mean(lake_rel_abund),sd=mean(lake_rel_sd)) %>%
  dplyr::group_by(month_full) %>%
  dplyr::mutate(.,sum_mean=sum(mean),sum_sd=sum(sd),nb_ASV=length(unique(.$OTU))) %>%
  ungroup() %>% dplyr::mutate(.,mean_lake=mean(sum_mean),mean_sd=mean(sum_sd),lake="CRJ1",status="Eutrophic")

all_ASV<- length(unique((subset(shared_lake.df, Abundance >0))$OTU))

shared_ASV.df[nrow(shared_ASV.df) + 1,] = c(unique(Shared_CRJ1$lake),
                                            unique(Shared_CRJ1$status),
                                            unique(Shared_CRJ1$nb_ASV),
                                            unique(Shared_CRJ1$mean_lake),
                                            unique(Shared_CRJ1$mean_sd),
                                            (unique(Shared_CRJ1$nb_ASV)/all_ASV)*100)
## CRJ2
shared_lake.df <- shared.df %>% subset(.,lake=="CRJ2")
lake_B=unique((subset(shared_lake.df,month=="B" & Abundance >0))$OTU)
lake_C=unique((subset(shared_lake.df,month=="C" & Abundance >0))$OTU)
lake_D=unique((subset(shared_lake.df,month=="D" & Abundance >0))$OTU)

Shared_CRJ2<-as.data.frame(Reduce(intersect, list(lake_A,lake_B,lake_C,lake_D))) %>%
  rename("ASV"="Reduce(intersect, list(lake_A, lake_B, lake_C, lake_D))")

Shared_CRJ2 <- shared_lake.df %>% subset(.,OTU %in% as.vector(Shared_CRJ2$ASV)) %>% .[,c(1,4,9,14,17,20,21,23,24,26,27)] %>%
  dplyr::group_by(OTU,month_full,Phylum) %>%
  dplyr::summarize(mean=mean(lake_rel_abund),sd=mean(lake_rel_sd)) %>%
  dplyr::group_by(month_full) %>%
  dplyr::mutate(.,sum_mean=sum(mean),sum_sd=sum(sd),nb_ASV=length(unique(.$OTU))) %>%
  ungroup() %>% dplyr::mutate(.,mean_lake=mean(sum_mean),mean_sd=mean(sum_sd),lake="CRJ2",status="Eutrophic")

all_ASV<- length(unique((subset(shared_lake.df, Abundance >0))$OTU))

shared_ASV.df[nrow(shared_ASV.df) + 1,] = c(unique(Shared_CRJ2$lake),
                                            unique(Shared_CRJ2$status),
                                            unique(Shared_CRJ2$nb_ASV),
                                            unique(Shared_CRJ2$mean_lake),
                                            unique(Shared_CRJ2$mean_sd),
                                            (unique(Shared_CRJ2$nb_ASV)/all_ASV)*100)
####__ Mesotrophic ####
## VAI
shared_lake.df <- shared.df %>% subset(.,lake=="VAI")
lake_A=unique((subset(shared_lake.df,month=="A" & Abundance >0))$OTU)
lake_B=unique((subset(shared_lake.df,month=="B" & Abundance >0))$OTU)
lake_C=unique((subset(shared_lake.df,month=="C" & Abundance >0))$OTU)
lake_D=unique((subset(shared_lake.df,month=="D" & Abundance >0))$OTU)

Shared_VAI<-as.data.frame(Reduce(intersect, list(lake_A,lake_B,lake_C,lake_D))) %>%
  rename("ASV"="Reduce(intersect, list(lake_A, lake_B, lake_C, lake_D))")

Shared_VAI <- shared_lake.df %>% subset(.,OTU %in% as.vector(Shared_VAI$ASV)) %>% .[,c(1,4,9,14,17,20,21,23,24,26,27)] %>%
  dplyr::group_by(OTU,month_full,Phylum) %>%
  dplyr::summarize(mean=mean(lake_rel_abund),sd=mean(lake_rel_sd)) %>%
  dplyr::group_by(month_full) %>%
  dplyr::mutate(.,sum_mean=sum(mean),sum_sd=sum(sd),nb_ASV=length(unique(.$OTU))) %>%
  ungroup() %>% dplyr::mutate(.,mean_lake=mean(sum_mean),mean_sd=mean(sum_sd),lake="VAI",status="Mesotrophic")

all_ASV<- length(unique((subset(shared_lake.df, Abundance >0))$OTU))

shared_ASV.df[nrow(shared_ASV.df) + 1,] = c(unique(Shared_VAI$lake),
                                            unique(Shared_VAI$status),
                                            unique(Shared_VAI$nb_ASV),
                                            unique(Shared_VAI$mean_lake),
                                            unique(Shared_VAI$mean_sd),
                                            (unique(Shared_VAI$nb_ASV)/all_ASV)*100)
## JAB
shared_lake.df <- shared.df %>% subset(.,lake=="JAB")
lake_A=unique((subset(shared_lake.df,month=="A" & Abundance >0))$OTU)
lake_B=unique((subset(shared_lake.df,month=="B" & Abundance >0))$OTU)
lake_C=unique((subset(shared_lake.df,month=="C" & Abundance >0))$OTU)
lake_D=unique((subset(shared_lake.df,month=="D" & Abundance >0))$OTU)

Shared_JAB<-as.data.frame(Reduce(intersect, list(lake_A,lake_B,lake_C,lake_D))) %>%
  rename("ASV"="Reduce(intersect, list(lake_A, lake_B, lake_C, lake_D))")

Shared_JAB <- shared_lake.df %>% subset(.,OTU %in% as.vector(Shared_JAB$ASV)) %>% .[,c(1,4,9,14,17,20,21,23,24,26,27)] %>%
  dplyr::group_by(OTU,month_full,Phylum) %>%
  dplyr::summarize(mean=mean(lake_rel_abund),sd=mean(lake_rel_sd)) %>%
  dplyr::group_by(month_full) %>%
  dplyr::mutate(.,sum_mean=sum(mean),sum_sd=sum(sd),nb_ASV=length(unique(.$OTU))) %>%
  ungroup() %>% dplyr::mutate(.,mean_lake=mean(sum_mean),mean_sd=mean(sum_sd),lake="JAB",status="Mesotrophic")

all_ASV<- length(unique((subset(shared_lake.df, Abundance >0))$OTU))

shared_ASV.df[nrow(shared_ASV.df) + 1,] = c(unique(Shared_JAB$lake),
                                            unique(Shared_JAB$status),
                                            unique(Shared_JAB$nb_ASV),
                                            unique(Shared_JAB$mean_lake),
                                            unique(Shared_JAB$mean_sd),
                                            (unique(Shared_JAB$nb_ASV)/all_ASV)*100)

VAI<-Shared_VAI$OTU
BOI<-Shared_BOI$OTU
CTL<-Shared_CTL$OTU
VER<-Shared_VER$OTU
CHA<-Shared_CHA$OTU
GDP<-Shared_GDP$OTU
JAB<-Shared_JAB$OTU
CRJ1<-Shared_CRJ1$OTU
CRJ2<-Shared_CRJ2$OTU

lake_shared_ASV<-as.data.frame(Reduce(intersect, list(VAI,BOI,CTL,VER,CHA,GDP,JAB,CRJ1,CRJ2))) %>%
  rename("ASV"="Reduce(intersect, list(VAI, BOI, CTL, VER, CHA, GDP, JAB, CRJ1, CRJ2))")
length(lake_shared_ASV$ASV)

#### Shared KOs  ####

shared.df<- all_fun.df %>% subset(.,month %in% c("A","B","C","D")) %>%
  dplyr::group_by(lake_month) %>%
  dplyr::mutate(.,sample_CPM=sum(CPM)) %>% ungroup() %>%
  dplyr::group_by(KEGG_ko,lake_month,sample_CPM) %>% dplyr::summarise(sum_CPM=sum(CPM))

shared.df<-
  shared.df %>% dplyr::group_by(KEGG_ko,lake_month) %>%
  dplyr::summarise(lake_rel_abund=(sum_CPM/sample_CPM)*100) %>%
  cbind(.,lake=all_fun.df$lake[match(.$lake_month,all_fun.df$lake_month)],
        month=all_fun.df$month[match(.$lake_month,all_fun.df$lake_month)]) #View(shared.df)

####__ Hypereutrophic ####
## VER
shared_lake.df <- shared.df %>% subset(.,lake=="VER")
lake_A=unique((subset(shared_lake.df,month=="A" & lake_rel_abund >0))$KEGG_ko)
lake_B=unique((subset(shared_lake.df,month=="B" & lake_rel_abund >0))$KEGG_ko)
lake_C=unique((subset(shared_lake.df,month=="C" & lake_rel_abund >0))$KEGG_ko)
lake_D=unique((subset(shared_lake.df,month=="D" & lake_rel_abund >0))$KEGG_ko)

Shared_VER<-as.data.frame(Reduce(intersect, list(lake_A,lake_B,lake_C,lake_D))) %>%
  rename("KO"="Reduce(intersect, list(lake_A, lake_B, lake_C, lake_D))")

Shared_VER <- shared_lake.df %>% subset(.,KEGG_ko %in% as.vector(Shared_VER$KO)) %>%
  dplyr::group_by(KEGG_ko,month) %>%
  dplyr::summarize(mean=mean(lake_rel_abund),sd=sd(lake_rel_abund)) %>%
  dplyr::group_by(month) %>%
  dplyr::mutate(.,sum_mean=sum(mean),sum_sd=sum(sd),nb_KO=length(unique(.$KEGG_ko))) %>%
  ungroup() %>% dplyr::mutate(.,mean_lake=median(sum_mean),mean_sd=mean(sum_sd),lake="VER",status="Hypereutrophic")

all_KO<- length(unique((subset(shared_lake.df, lake_rel_abund >0))$KEGG_ko))

shared_KO.df<-data.frame(lake=unique(Shared_VER$lake),status=unique(Shared_VER$status),
                         nb_KO=unique(Shared_VER$nb_KO),
                         rel.abund=unique(Shared_VER$mean_lake),
                         per_lake=(unique(Shared_VER$nb_KO)/all_KO)*100)

## CHA
shared_lake.df <- shared.df %>% subset(.,lake=="CHA")
lake_A=unique((subset(shared_lake.df,month=="A" & lake_rel_abund >0))$KEGG_ko)
lake_B=unique((subset(shared_lake.df,month=="B" & lake_rel_abund >0))$KEGG_ko)
lake_C=unique((subset(shared_lake.df,month=="C" & lake_rel_abund >0))$KEGG_ko)
lake_D=unique((subset(shared_lake.df,month=="D" & lake_rel_abund >0))$KEGG_ko)

Shared_CHA<-as.data.frame(Reduce(intersect, list(lake_A,lake_B,lake_C,lake_D))) %>%
  rename("KO"="Reduce(intersect, list(lake_A, lake_B, lake_C, lake_D))")

Shared_CHA <- shared_lake.df %>% subset(.,KEGG_ko %in% as.vector(Shared_CHA$KO)) %>%
  dplyr::group_by(KEGG_ko,month) %>%
  dplyr::summarize(mean=mean(lake_rel_abund),sd=sd(lake_rel_abund)) %>%
  dplyr::group_by(month) %>%
  dplyr::mutate(.,sum_mean=sum(mean),sum_sd=sum(sd),nb_KO=length(unique(.$KEGG_ko))) %>%
  ungroup() %>% dplyr::mutate(.,mean_lake=median(sum_mean),mean_sd=mean(sum_sd),lake="CHA",status="Hypereutrophic")

all_KO<- length(unique((subset(shared_lake.df, lake_rel_abund >0))$KEGG_ko))

shared_KO.df[nrow(shared_KO.df) + 1,] = c(unique(Shared_CHA$lake),
                                          unique(Shared_CHA$status),
                                          unique(Shared_CHA$nb_KO),
                                          unique(Shared_CHA$mean_lake),
                                          (unique(Shared_CHA$nb_KO)/all_KO)*100)
####__ Eutrophic ####
## GDP
shared_lake.df <- shared.df %>% subset(.,lake=="GDP")
lake_A=unique((subset(shared_lake.df,month=="A" & lake_rel_abund >0))$KEGG_ko)
lake_B=unique((subset(shared_lake.df,month=="B" & lake_rel_abund >0))$KEGG_ko)
lake_C=unique((subset(shared_lake.df,month=="C" & lake_rel_abund >0))$KEGG_ko)
lake_D=unique((subset(shared_lake.df,month=="D" & lake_rel_abund >0))$KEGG_ko)

Shared_GDP<-as.data.frame(Reduce(intersect, list(lake_A,lake_B,lake_C,lake_D))) %>%
  rename("KO"="Reduce(intersect, list(lake_A, lake_B, lake_C, lake_D))")

Shared_GDP <- shared_lake.df %>% subset(.,KEGG_ko %in% as.vector(Shared_GDP$KO)) %>%
  dplyr::group_by(KEGG_ko,month) %>%
  dplyr::summarize(mean=mean(lake_rel_abund),sd=sd(lake_rel_abund)) %>%
  dplyr::group_by(month) %>%
  dplyr::mutate(.,sum_mean=sum(mean),sum_sd=sum(sd),nb_KO=length(unique(.$KEGG_ko))) %>%
  ungroup() %>% dplyr::mutate(.,mean_lake=median(sum_mean),mean_sd=mean(sum_sd),lake="GDP",status="Eutrophic")

all_KO<- length(unique((subset(shared_lake.df, lake_rel_abund >0))$KEGG_ko))

shared_KO.df[nrow(shared_KO.df) + 1,] = c(unique(Shared_GDP$lake),
                                          unique(Shared_GDP$status),
                                          unique(Shared_GDP$nb_KO),
                                          unique(Shared_GDP$mean_lake),
                                          (unique(Shared_GDP$nb_KO)/all_KO)*100)

## BOI
shared_lake.df <- shared.df %>% subset(.,lake=="BOI")
lake_A=unique((subset(shared_lake.df,month=="A" & lake_rel_abund >0))$KEGG_ko)
lake_B=unique((subset(shared_lake.df,month=="B" & lake_rel_abund >0))$KEGG_ko)
lake_C=unique((subset(shared_lake.df,month=="C" & lake_rel_abund >0))$KEGG_ko)
lake_D=unique((subset(shared_lake.df,month=="D" & lake_rel_abund >0))$KEGG_ko)

Shared_BOI<-as.data.frame(Reduce(intersect, list(lake_A,lake_B,lake_C,lake_D))) %>%
  rename("KO"="Reduce(intersect, list(lake_A, lake_B, lake_C, lake_D))")

Shared_BOI <- shared_lake.df %>% subset(.,KEGG_ko %in% as.vector(Shared_BOI$KO)) %>%
  dplyr::group_by(KEGG_ko,month) %>%
  dplyr::summarize(mean=mean(lake_rel_abund),sd=sd(lake_rel_abund)) %>%
  dplyr::group_by(month) %>%
  dplyr::mutate(.,sum_mean=sum(mean),sum_sd=sum(sd),nb_KO=length(unique(.$KEGG_ko))) %>%
  ungroup() %>% dplyr::mutate(.,mean_lake=median(sum_mean),mean_sd=mean(sum_sd),lake="BOI",status="Eutrophic")

all_KO<- length(unique((subset(shared_lake.df, lake_rel_abund >0))$KEGG_ko))

shared_KO.df[nrow(shared_KO.df) + 1,] = c(unique(Shared_BOI$lake),
                                          unique(Shared_BOI$status),
                                          unique(Shared_BOI$nb_KO),
                                          unique(Shared_BOI$mean_lake),
                                          (unique(Shared_BOI$nb_KO)/all_KO)*100)
## CTL
shared_lake.df <- shared.df %>% subset(.,lake=="CTL")
lake_A=unique((subset(shared_lake.df,month=="A" & lake_rel_abund >0))$KEGG_ko)
lake_B=unique((subset(shared_lake.df,month=="B" & lake_rel_abund >0))$KEGG_ko)
lake_C=unique((subset(shared_lake.df,month=="C" & lake_rel_abund >0))$KEGG_ko)
lake_D=unique((subset(shared_lake.df,month=="D" & lake_rel_abund >0))$KEGG_ko)

Shared_CTL<-as.data.frame(Reduce(intersect, list(lake_A,lake_B,lake_C,lake_D))) %>%
  rename("KO"="Reduce(intersect, list(lake_A, lake_B, lake_C, lake_D))")

Shared_CTL <- shared_lake.df %>% subset(.,KEGG_ko %in% as.vector(Shared_CTL$KO)) %>%
  dplyr::group_by(KEGG_ko,month) %>%
  dplyr::summarize(mean=mean(lake_rel_abund),sd=sd(lake_rel_abund)) %>%
  dplyr::group_by(month) %>%
  dplyr::mutate(.,sum_mean=sum(mean),sum_sd=sum(sd),nb_KO=length(unique(.$KEGG_ko))) %>%
  ungroup() %>% dplyr::mutate(.,mean_lake=median(sum_mean),mean_sd=mean(sum_sd),lake="CTL",status="Eutrophic")

all_KO<- length(unique((subset(shared_lake.df, lake_rel_abund >0))$KEGG_ko))

shared_KO.df[nrow(shared_KO.df) + 1,] = c(unique(Shared_CTL$lake),
                                          unique(Shared_CTL$status),
                                          unique(Shared_CTL$nb_KO),
                                          unique(Shared_CTL$mean_lake),
                                          (unique(Shared_CTL$nb_KO)/all_KO)*100)
## CRJ1
shared_lake.df <- shared.df %>% subset(.,lake=="CRJ1")
lake_A=unique((subset(shared_lake.df,month=="A" & lake_rel_abund >0))$KEGG_ko)
lake_B=unique((subset(shared_lake.df,month=="B" & lake_rel_abund >0))$KEGG_ko)
lake_C=unique((subset(shared_lake.df,month=="C" & lake_rel_abund >0))$KEGG_ko)
lake_D=unique((subset(shared_lake.df,month=="D" & lake_rel_abund >0))$KEGG_ko)

Shared_CRJ1<-as.data.frame(Reduce(intersect, list(lake_A,lake_B,lake_C,lake_D))) %>%
  rename("KO"="Reduce(intersect, list(lake_A, lake_B, lake_C, lake_D))")

Shared_CRJ1 <- shared_lake.df %>% subset(.,KEGG_ko %in% as.vector(Shared_CRJ1$KO)) %>%
  dplyr::group_by(KEGG_ko,month) %>%
  dplyr::summarize(mean=mean(lake_rel_abund),sd=sd(lake_rel_abund)) %>%
  dplyr::group_by(month) %>%
  dplyr::mutate(.,sum_mean=sum(mean),sum_sd=sum(sd),nb_KO=length(unique(.$KEGG_ko))) %>%
  ungroup() %>% dplyr::mutate(.,mean_lake=median(sum_mean),mean_sd=mean(sum_sd),lake="CRJ1",status="Eutrophic")

all_KO<- length(unique((subset(shared_lake.df, lake_rel_abund >0))$KEGG_ko))

shared_KO.df[nrow(shared_KO.df) + 1,] = c(unique(Shared_CRJ1$lake),
                                          unique(Shared_CRJ1$status),
                                          unique(Shared_CRJ1$nb_KO),
                                          unique(Shared_CRJ1$mean_lake),
                                          (unique(Shared_CRJ1$nb_KO)/all_KO)*100)
## CRJ2
shared_lake.df <- shared.df %>% subset(.,lake=="CRJ2")
lake_B=unique((subset(shared_lake.df,month=="B" & lake_rel_abund >0))$KEGG_ko)
lake_C=unique((subset(shared_lake.df,month=="C" & lake_rel_abund >0))$KEGG_ko)
lake_D=unique((subset(shared_lake.df,month=="D" & lake_rel_abund >0))$KEGG_ko)

Shared_CRJ2<-as.data.frame(Reduce(intersect, list(lake_B,lake_C,lake_D))) %>%
  rename("KO"="Reduce(intersect, list(lake_B, lake_C, lake_D))")

Shared_CRJ2 <- shared_lake.df %>% subset(.,KEGG_ko %in% as.vector(Shared_CRJ2$KO)) %>%
  dplyr::group_by(KEGG_ko,month) %>%
  dplyr::summarize(mean=mean(lake_rel_abund),sd=sd(lake_rel_abund)) %>%
  dplyr::group_by(month) %>%
  dplyr::mutate(.,sum_mean=sum(mean),sum_sd=sum(sd),nb_KO=length(unique(.$KEGG_ko))) %>%
  ungroup() %>% dplyr::mutate(.,mean_lake=median(sum_mean),mean_sd=mean(sum_sd),lake="CRJ2",status="Eutrophic")

all_KO<- length(unique((subset(shared_lake.df, lake_rel_abund >0))$KEGG_ko))

shared_KO.df[nrow(shared_KO.df) + 1,] = c(unique(Shared_CRJ2$lake),
                                          unique(Shared_CRJ2$status),
                                          unique(Shared_CRJ2$nb_KO),
                                          unique(Shared_CRJ2$mean_lake),
                                          (unique(Shared_CRJ2$nb_KO)/all_KO)*100)
####__ Mesotrophic ####
## VAI
shared_lake.df <- shared.df %>% subset(.,lake=="VAI")
lake_A=unique((subset(shared_lake.df,month=="A" & lake_rel_abund >0))$KEGG_ko)
lake_B=unique((subset(shared_lake.df,month=="B" & lake_rel_abund >0))$KEGG_ko)
lake_C=unique((subset(shared_lake.df,month=="C" & lake_rel_abund >0))$KEGG_ko)
lake_D=unique((subset(shared_lake.df,month=="D" & lake_rel_abund >0))$KEGG_ko)

Shared_VAI<-as.data.frame(Reduce(intersect, list(lake_A,lake_B,lake_C,lake_D))) %>%
  rename("KO"="Reduce(intersect, list(lake_A, lake_B, lake_C, lake_D))")

Shared_VAI <- shared_lake.df %>% subset(.,KEGG_ko %in% as.vector(Shared_VAI$KO)) %>%
  dplyr::group_by(KEGG_ko,month) %>%
  dplyr::summarize(mean=mean(lake_rel_abund),sd=sd(lake_rel_abund)) %>%
  dplyr::group_by(month) %>%
  dplyr::mutate(.,sum_mean=sum(mean),sum_sd=sum(sd),nb_KO=length(unique(.$KEGG_ko))) %>%
  ungroup() %>% dplyr::mutate(.,mean_lake=median(sum_mean),mean_sd=mean(sum_sd),lake="VAI",status="Mesotrophic")

all_KO<- length(unique((subset(shared_lake.df, lake_rel_abund >0))$KEGG_ko))

shared_KO.df[nrow(shared_KO.df) + 1,] = c(unique(Shared_VAI$lake),
                                          unique(Shared_VAI$status),
                                          unique(Shared_VAI$nb_KO),
                                          unique(Shared_VAI$mean_lake),
                                          (unique(Shared_VAI$nb_KO)/all_KO)*100)
## JAB
shared_lake.df <- shared.df %>% subset(.,lake=="JAB")
lake_A=unique((subset(shared_lake.df,month=="A" & lake_rel_abund >0))$KEGG_ko)
lake_B=unique((subset(shared_lake.df,month=="B" & lake_rel_abund >0))$KEGG_ko)
lake_C=unique((subset(shared_lake.df,month=="C" & lake_rel_abund >0))$KEGG_ko)
lake_D=unique((subset(shared_lake.df,month=="D" & lake_rel_abund >0))$KEGG_ko)

Shared_JAB<-as.data.frame(Reduce(intersect, list(lake_A,lake_B,lake_C,lake_D))) %>%
  rename("KO"="Reduce(intersect, list(lake_A, lake_B, lake_C, lake_D))")

Shared_JAB <- shared_lake.df %>% subset(.,KEGG_ko %in% as.vector(Shared_JAB$KO)) %>%
  dplyr::group_by(KEGG_ko,month) %>%
  dplyr::summarize(mean=mean(lake_rel_abund),sd=sd(lake_rel_abund)) %>%
  dplyr::group_by(month) %>%
  dplyr::mutate(.,sum_mean=sum(mean),sum_sd=sum(sd),nb_KO=length(unique(.$KEGG_ko))) %>%
  ungroup() %>% dplyr::mutate(.,mean_lake=median(sum_mean),mean_sd=mean(sum_sd),lake="JAB",status="Mesotrophic")

all_KO<- length(unique((subset(shared_lake.df, lake_rel_abund >0))$KEGG_ko))

shared_KO.df[nrow(shared_KO.df) + 1,] = c(unique(Shared_JAB$lake),
                                          unique(Shared_JAB$status),
                                          unique(Shared_JAB$nb_KO),
                                          unique(Shared_JAB$mean_lake),
                                          (unique(Shared_JAB$nb_KO)/all_KO)*100)

## shared in all lakes and all months
VAI<-Shared_VAI$KEGG_ko
BOI<-Shared_BOI$KEGG_ko
CTL<-Shared_CTL$KEGG_ko
VER<-Shared_VER$KEGG_ko
CHA<-Shared_CHA$KEGG_ko
GDP<-Shared_GDP$KEGG_ko
JAB<-Shared_JAB$KEGG_ko
CRJ1<-Shared_CRJ1$KEGG_ko
CRJ2<-Shared_CRJ2$KEGG_ko

lake_shared_KO<-as.data.frame(Reduce(intersect, list(VAI,BOI,CTL,VER,CHA,GDP,JAB,CRJ1,CRJ2))) %>%
  rename("KO"="Reduce(intersect, list(VAI, BOI, CTL, VER, CHA, GDP, JAB, CRJ1, CRJ2))")
length(lake_shared_KO$KO)

#### Overall shared ####

shared_ASV.df <- shared_ASV.df %>%
  dplyr::mutate(.,lake=
                  if_else(lake=="VER","Verneuil",
                          if_else(lake=="CHA","Champs",
                                  if_else(lake=="GDP","Grande-Paroisse",
                                          if_else(lake=="BOI","Bois",
                                                  if_else(lake=="CTL","Créteil",
                                                          if_else(lake=="CRJ1","Cergy-small",
                                                                  if_else(lake=="CRJ2","Cergy-large",
                                                                          if_else(lake=="VAI","Vaires","Jablines")))))))),
                lake=factor(lake,levels=(rev(c("Verneuil","Champs",
                                           "Grande-Paroisse","Bois","Créteil","Cergy-small","Cergy-large",
                                           "Vaires","Jablines")))),
                rel.abund=as.numeric(rel.abund),
                rel.sd=as.numeric(rel.sd))
shared_ASV.df <- shared_ASV.df %>% dplyr::mutate(.,features="ASVs",lake_features=paste0(lake,"_",features))
shared_ASV.df<- shared_ASV.df %>% rename("nb_features"="nb_ASV")

shared_KO.df <- shared_KO.df %>%
  dplyr::mutate(.,lake=
                  if_else(lake=="VER","Verneuil",
                          if_else(lake=="CHA","Champs",
                                  if_else(lake=="GDP","Grande-Paroisse",
                                          if_else(lake=="BOI","Bois",
                                                  if_else(lake=="CTL","Créteil",
                                                          if_else(lake=="CRJ1","Cergy-small",
                                                                  if_else(lake=="CRJ2","Cergy-large",
                                                                          if_else(lake=="VAI","Vaires","Jablines")))))))),
                lake=factor(lake,levels=(rev(c("Verneuil","Champs",
                                               "Grande-Paroisse","Bois","Créteil","Cergy-small","Cergy-large",
                                               "Vaires","Jablines")))),
                rel.abund=as.numeric(rel.abund),
                features="KOs",lake_features=paste0(lake,"_",features))
shared_KO.df<- shared_KO.df %>% rename("nb_features"="nb_KO")

Fig_shared <-
  rbind(shared_ASV.df[,-5],shared_KO.df) %>%
  dplyr::mutate(.,rel.abund=as.numeric(rel.abund),
                lake=factor(lake,levels=(rev(c("Verneuil","Champs",
                                                 "Grande-Paroisse","Bois","Créteil","Cergy-small","Cergy-large",
                                                 "Vaires","Jablines"))))) %>%
  ggplot(.,aes(x=lake,y=features,size = rel.abund,label=nb_features)) +
  geom_point(color=c("#28A448","#1E7A36","#980000","#E40000","#E45000","#EB7947","#F2AB8C","#5C9BD6","#2F6B9D",
                     "#28A448","#1E7A36","#980000","#E40000","#E45000","#EB7947","#F2AB8C","#5C9BD6","#2F6B9D"))+
  geom_text(size=6,color="white",show.legend = F,fontface="bold")+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.ticks = element_blank(),
        axis.text.y = element_text(size=15,face = "bold",angle=90,hjust=0.5),
        axis.text.x = element_blank(),
        legend.title = element_text(size=15,face = "bold"),
        legend.text = element_text(size=12),
        plot.title = element_text(hjust=0,size=15),
        axis.title.y = element_blank(),
        axis.title.x = element_blank())+
  scale_size_continuous(range = c(10,35),breaks = c(30,50,90),limits = c(30,NA))+
  labs(size="Relative\nabundance")

Fig_shared
#### Alpha diversity  ####

####__ Indices df ####
alphaDiv_summer<- as.data.frame(A_D_rarW@sam_data) %>%
  cbind(.,
        estimate_richness(A_D_rarW, measures=c("Shannon")),
        estimate_richness(A_D_rarW, measures=c("Observed")),
        estimate_richness(A_D_rarW, measures=c("Chao1"))) %>%
  dplyr::rename("shannon"="Shannon","richness"="Observed","chao1"="Chao1") %>%
  dplyr::mutate(.,
                evenness=as.numeric(shannon/log(richness)),
                Gap_CR=(chao1-richness),
                Gap_CR_percent=(((chao1-richness)/chao1)*100),
                completude=(chao1/richness),
                lake=
                  if_else(lake=="VER","Verneuil",
                          if_else(lake=="CHA","Champs",
                                  if_else(lake=="GDP","Grande-Paroisse",
                                          if_else(lake=="BOI","Bois",
                                                  if_else(lake=="CTL","Creteil",
                                                          if_else(lake=="CRJ1","Cergy-small",
                                                                  if_else(lake=="CRJ2","Cergy-large",
                                                                          if_else(lake=="VAI","Vaires","Jablines")))))))),
                lake=factor(lake,levels=(c("Verneuil","Champs",
                                           "Grande-Paroisse","Bois","Creteil","Cergy-small","Cergy-large",
                                           "Vaires","Jablines"))))

#View(alphaDiv_summer)

write.csv(alphaDiv_summer,
          "alphaDiv_summer_tmp.csv")
wilcox.test(richness ~ lake, data = alphaDiv_summer, paired = TRUE)


alphaDiv_summer_tmp<-read.delim("alphaDiv_summer_tmp.csv",sep = ";")
View(alphaDiv_summer_tmp)
lm(richness~lake+ (1|lake:month),data=alphaDiv_summer)
library(emmeans)

testlm<-lmer(richness~ month + lake + (1|lake),data=alphaDiv_summer,REML = F)
summary(testlm)

?emmeans
emm_lake<-pairs(emmeans(testlm,"lake"))
summ_lmw<-summary(emm_lake,adjust="BH")
View(summ_lmw)

emmeans(testlm, specs = pairwise ~ lake:lake, type = "response")

#Stats richness
alphaDiv_summer.stat<-alphaDiv_summer %>%
  wilcox_test(richness ~ lake,p.adjust.method = "BH",paired = T) %>%
  add_significance("p.adj")  %>% add_xy_position(x = "lake_order", dodge = 0.8) %>%
  mutate(g1=gsub("_","",as.character(.$group1)),g2=gsub("_","",as.character(.$group2))) %>%
  unite("comp",g1:g2,sep=" - ",remove = FALSE)

letters_richness<-tibble(lake_order=unique(sort(alphaDiv_summer$lake_order)),
                         letters=cldList(p.adj ~ as.character(comp),alphaDiv_summer.stat,threshold  = 0.05)$Letter)
alphaDiv_summer$lake_order
#Stats shannon
alphaDiv_summer.stat<-alphaDiv_summer %>% wilcox_test(shannon ~ lake_order,p.adjust.method = "BH") %>%
  add_significance("p.adj")  %>% add_xy_position(x = "lake_order", dodge = 0.8) %>%
  mutate(g1=gsub("_","",as.character(.$group1)),g2=gsub("_","",as.character(.$group2))) %>%
  unite("comp",g1:g2,sep=" - ",remove = FALSE)

letters_shannon<-tibble(lake_order=unique(sort(alphaDiv_summer$lake_order)),
                        letters=cldList(p.adj ~ as.character(comp),alphaDiv_summer.stat,threshold  = 0.05)$Letter)

#Stats evenness
alphaDiv_summer.stat<-alphaDiv_summer %>% wilcox_test(evenness ~ lake_order,p.adjust.method = "BH") %>%
  add_significance("p.adj")  %>% add_xy_position(x = "lake_order", dodge = 0.8) %>%
  mutate(g1=gsub("_","",as.character(.$group1)),g2=gsub("_","",as.character(.$group2))) %>%
  unite("comp",g1:g2,sep=" - ",remove = FALSE)

letters_evenness<-tibble(lake_order=unique(sort(alphaDiv_summer$lake_order)),
                         letters=cldList(p.adj ~ as.character(comp),alphaDiv_summer.stat,threshold  = 0.05)$Letter)

####__ Indices plot ####

Fig_box_ASV_H_lake<- alphaDiv_summer %>%
  dplyr::mutate(.,lake=
                  if_else(lake=="Champs","CSM",
                          if_else(lake=="Verneuil","VSS",
                                  if_else(lake=="Grande-Paroisse","LGP",
                                          if_else(lake=="Bois","BLR",
                                                  if_else(lake=="Creteil","CRE",
                                                          if_else(lake=="Cergy-small","CER-S",
                                                                  if_else(lake=="Cergy-large","CER-L",
                                                                          if_else(lake=="Vaires","VSM","JAB")))))))),
                lake=factor(lake,levels=rev(c("VSS","CSM",
                                                        "LGP","BLR","CRE",
                                                        "CER-S","CER-L",
                                                        "VSM","JAB")))) %>%
  #cbind(.,letter=letters_shannon$letters[match(.$lake_order,letters_shannon$lake_order)]) %>%
  ggplot(.,aes(x=lake,y=shannon,fill=lake))+theme_bw()+
  geom_boxplot(show.legend = T,width=0.6,color="black")+
  # geom_text(mapping=aes(y=1.01*max(shannon),
  #                       label=letter),vjust=0,size=4.5,color="black",fontface="bold")+
  theme(panel.grid = element_blank(),
        axis.ticks = element_blank(),axis.title.x = element_blank(),
        axis.title.y = element_text(hjust=0.5,size=12,face = "bold"),
        axis.text.y = element_text(size=10),axis.text.x = element_blank(),
        legend.title = element_text(size=15,face = "bold"),legend.text = element_text(size=12),
        legend.position = "bottom",legend.direction = "vertical")+
  scale_fill_manual(values = rev(palette_lake_3T))+
  scale_y_continuous(limits = c(min(alphaDiv_summer$shannon)/1.01,max(alphaDiv_summer$shannon)*1.02),
                     breaks = c(3,3.5,4,4.5))+
  guides(fill = guide_legend(nrow=3,override.aes = list(size=5,shape=22)))+
  labs(y="Shannon",fill="Lake")

Fig_box_ASV_H_lake$data %>% dplyr::group_by(lake) %>%
  dplyr::summarise_at(vars("shannon"),list(median=median, sd=sd))

Fig_box_ASV_H_lake$data %>% subset(.,lake != "Jablines" & lake != "Bois") %>%
  dplyr::summarise_at(vars("shannon"),list(median=median, sd=sd))

Fig_box_ASV_S_lake<- alphaDiv_summer %>%
  dplyr::mutate(.,lake=
                  if_else(lake=="Champs","CSM",
                          if_else(lake=="Verneuil","VSS",
                                  if_else(lake=="Grande-Paroisse","LGP",
                                          if_else(lake=="Bois","BLR",
                                                  if_else(lake=="Creteil","CRE",
                                                          if_else(lake=="Cergy-small","CER-S",
                                                                  if_else(lake=="Cergy-large","CER-L",
                                                                          if_else(lake=="Vaires","VSM","JAB")))))))),
                lake=factor(lake,levels=rev(c("VSS","CSM",
                                              "LGP","BLR","CRE",
                                              "CER-S","CER-L",
                                              "VSM","JAB")))) %>%
  #cbind(.,letter=letters_evenness$letters[match(.$lake_order,letters_evenness$lake_order)]) %>%
  ggplot(.,aes(x=lake,y=evenness,fill=lake))+theme_bw()+
  geom_boxplot(show.legend = F,width=0.6,color="black")+
  # geom_text(mapping=aes(y=1.01*max(evenness),
  #                       label=letter),vjust=0,size=4.5,color="black",fontface="bold")+
  theme(panel.grid = element_blank(),
        axis.ticks = element_blank(),axis.title.x = element_blank(),
        axis.title.y = element_text(hjust=0.5,size=12,face = "bold"),
        axis.text.y = element_text(size=10),axis.text.x = element_blank(),
        legend.title = element_text(size=15,face = "bold"),legend.text = element_text(size=12))+
  scale_fill_manual(values = rev(palette_lake_3T))+
  scale_y_continuous(limits = c(min(alphaDiv_summer$evenness)/1.01,max(alphaDiv_summer$evenness)*1.02),
                     breaks = c(0.5,0.6,0.7,0.8,0.9))+
  guides(fill = guide_legend(nrow=3,override.aes = list(size=5,shape=22)))+
  labs(y="Evenness",fill="Lake")

Fig_box_ASV_S_lake$data %>% dplyr::group_by(lake) %>%
  dplyr::summarise_at(vars("evenness"),list(median=median, sd=sd))

Fig_box_ASV_S_lake$data %>% subset(.,lake != "Bois") %>%
  dplyr::summarise_at(vars("evenness"),list(median=median, sd=sd))

Fig_box_ASV_D_lake<- alphaDiv_summer %>%
  dplyr::mutate(.,lake=
                  if_else(lake=="Champs","CSM",
                          if_else(lake=="Verneuil","VSS",
                                  if_else(lake=="Grande-Paroisse","LGP",
                                          if_else(lake=="Bois","BLR",
                                                  if_else(lake=="Creteil","CRE",
                                                          if_else(lake=="Cergy-small","CER-S",
                                                                  if_else(lake=="Cergy-large","CER-L",
                                                                          if_else(lake=="Vaires","VSM","JAB")))))))),
                lake=factor(lake,levels=rev(c("VSS","CSM",
                                              "LGP","BLR","CRE",
                                              "CER-S","CER-L",
                                              "VSM","JAB")))) %>%
  #cbind(.,letter=letters_richness$letters[match(.$lake_order,letters_richness$lake_order)]) %>%
  ggplot(.,aes(x=lake,y=richness,fill=lake))+theme_bw()+
  geom_boxplot(show.legend = F,width=0.6,color="black")+
  # geom_text(mapping=aes(y=1.01*max(richness),
  #                       label=letter),vjust=0,size=4.5,color="black",fontface="bold")+
  theme(panel.grid = element_blank(),
        axis.ticks = element_blank(),axis.title.x = element_blank(),
        axis.title.y = element_text(hjust=0.5,size=12,face = "bold"),
        axis.text.y = element_text(size=10),axis.text.x = element_blank(),
        legend.title = element_text(size=15,face = "bold"),legend.text = element_text(size=12))+
  scale_fill_manual(values = rev(palette_lake_3T))+
  scale_y_continuous(limits = c(min(alphaDiv_summer$richness)/1.01,max(alphaDiv_summer$richness)*1.02),
                     breaks = c(100,200,300,400))+
  guides(fill = guide_legend(nrow=3,override.aes = list(size=5,shape=22)))+
  labs(y="ASVs Richness",fill="Lake")

Fig_box_ASV_D_lake$data %>%
  #dplyr::group_by(lake) %>%
  dplyr::summarise(median=,
                   sd=sd(richness)
                   )

median(Fig_box_ASV_S_lake$data$evenness)
sd(Fig_box_ASV_S_lake$data$evenness)
max(Fig_box_ASV_S_lake$data$evenness)
min(Fig_box_ASV_S_lake$data$evenness)
median(Fig_box_ASV_H_lake$data$shannon)
sd(Fig_box_ASV_H_lake$data$shannon)
median(Fig_box_ASV_D_lake$data$richness)
sd(Fig_box_ASV_D_lake$data$richness)
 # dplyr::summarise_at(vars("richness"),list(median=median, sd=sd))

Fig_box_ASV_D_lake$data %>% subset(.,lake != "Jablines") %>%
  dplyr::summarise_at(vars("richness"),list(median=median, sd=sd))

####__ overall plot ####

design_alpha<-"
A
B
C"

Fig_box_ASV_D_lake+Fig_box_ASV_S_lake+Fig_box_ASV_H_lake+
  plot_layout(guides = 'collect',design = design_alpha) +
  plot_annotation(tag_levels = 'A') &
  theme(plot.tag = element_text(size=15,face="bold"),
        legend.position = "bottom",legend.direction = "vertical")

####|####

#### ASVs Turnover ####

Summer_codyn.df <- A_D_rarW %>%
  psmelt(.) %>% dplyr::group_by(lake_month,OTU,.drop = FALSE,.add = TRUE) %>%
  summarise(median_count=median(Abundance),.groups = "keep") %>% dplyr::group_by(lake_month) %>%
  mutate(median_abundance=as.numeric(paste0((round(median_count/sum(median_count),4))*100))) %>%
  cbind(lake=
          (A_D_rarW@sam_data$lake_order[match(.$lake_month,A_D_rarW@sam_data$lake_month)]),
        month_full=
          (A_D_rarW@sam_data$month_full[match(.$lake_month,A_D_rarW@sam_data$lake_month)]),
        lake_type=
          (A_D_rarW@sam_data$lake_type[match(.$lake_month,A_D_rarW@sam_data$lake_month)]),.) %>%
  dplyr::mutate(.,
                ASV_codyn=ifelse(median_abundance<1,"ASV_<_1%",OTU),
                sample_codyn=paste0(lake,"_",month_full,"_",ASV_codyn))

ASV_codyn<-as.data.frame(cbind(
  sample_codyn=Summer_codyn.df$sample_codyn,
  median_abund=Summer_codyn.df$median_abundance)) %>%
  dplyr::mutate(.,median_abund=as.numeric(median_abund)) %>%
  dplyr::group_by(sample_codyn) %>%
  dplyr::summarise(median_rel_abund=sum(median_abund)) %>%
  cbind(lake=Summer_codyn.df$lake[match(.$sample_codyn,Summer_codyn.df$sample_codyn)],
        month_full=Summer_codyn.df$month_full[match(.$sample_codyn,Summer_codyn.df$sample_codyn)],
        ASV_codyn=Summer_codyn.df$ASV_codyn[match(.$sample_codyn,Summer_codyn.df$sample_codyn)],
        median_rel_abund=.$median_rel_abund) %>% .[,-c(1,2)] %>%
  dplyr::mutate(.,
                month_number=A_D_rarW@sam_data$month_number[match(.$month_full,A_D_rarW@sam_data$month_full)])

ASV_turnover<-turnover(ASV_codyn,time.var = "month_number",species.var = "ASV_codyn",
                       abundance.var = "median_rel_abund",replicate.var = "lake",metric = "total") %>%
  dplyr::rename("ASV_turnover"="total") %>%
  dplyr::mutate(.,
                lake_month=paste0(lake,"_",month_number),
                Status=A_D_rarW@sam_data$lake_type[match(.$lake,A_D_rarW@sam_data$lake)],
                lake=factor(.$lake,levels=c("CHA","VER","GDP","BOI","CTL","CRJ1","CRJ2","VAI","JAB"))) #
View(ASV_turnover)

compare_turnover.stats<-ASV_turnover %>% wilcox_test(ASV_turnover ~ Status,p.adjust.method = "BH") %>%
  add_significance("p.adj")  %>% add_xy_position(x = "Status", dodge = 0.8) %>%
  mutate(g1=gsub("_","",as.character(.$group1)),g2=gsub("_","",as.character(.$group2))) %>%
  unite("comp",g1:g2,sep=" - ",remove = FALSE)

letters_box_ASV_turnover<-tibble(Status=unique(sort(ASV_turnover$Status)),
                                 letters=cldList(p.adj ~ as.character(comp),compare_turnover.stats,threshold  = 0.05)$Letter)

Fig_box_ASV_turnover<- ASV_turnover %>%
  dplyr::group_by(Status) %>% dplyr::mutate(.,mean_turnover=mean(ASV_turnover)) %>%
  cbind(.,letter=letters_box_ASV_turnover$letters[match(.$Status,letters_box_ASV_turnover$Status)]) %>%
  ggplot(.,aes(x=fct_rev(fct_reorder(Status,mean_turnover)),y=ASV_turnover,fill=Status))+theme_bw()+
  geom_boxplot(show.legend = T,width=0.6,color="black")+
  geom_text(mapping=aes(y=1.01*max(ASV_turnover),
                        label=letter),vjust=0,size=4,color="darkred",fontface="bold")+
  theme(panel.grid = element_blank(),plot.title = element_text(hjust=0.5,size=15,face = "bold"),
        axis.ticks = element_blank(),axis.title.y = element_blank(),
        axis.title.x = element_text(hjust=0.5,size=12,face = "bold"),
        axis.text.y = element_blank(),axis.text.x = element_text(size=10),
        plot.caption = element_text(size=12,face="italic",hjust = 1),
        legend.title = element_text(size=15,face = "bold"),legend.text = element_text(size=12))+
  scale_fill_manual(values = palette_status)+
  labs(y="ASV turnover")+coord_flip()
Fig_box_ASV_turnover

View(Summer_codyn.df)

write.csv(ASV_KO_turnover,"ASV_KO_turnover.csv")

Fig_box_ASV_KO_turnover<- read.delim("ASV_KO_turnover_box.csv",sep = ";") %>%
  dplyr::mutate(.,lake=factor(lake,levels=c("CHA","VER","GDP","BOI","CTL","CRJ1","CRJ2","VAI","JAB")),
                  Status=factor(Status,levels=c("Hypereutrophic","Eutrophic","Mesotrophic"))) %>%
  ggplot(.,aes(x=Status_features,y=turnover,fill=features))+theme_bw()+
  geom_boxplot(width=0.6,color="black")+
  theme(panel.grid = element_blank(),plot.title = element_text(hjust=0.5,size=15,face = "bold"),
        axis.ticks = element_blank(),axis.title.x = element_blank(),
        axis.title.y = element_text(hjust=0.5,size=15,face = "bold"),
        axis.text.x = element_blank(),axis.text.y = element_text(size=12),
        legend.title = element_text(size=15,face = "bold"),legend.text = element_text(size=12))+
  scale_fill_manual(values = c("blue","orange"))+
  facet_wrap2(~ Status,scale="free_x",
              strip = strip_color_lake_afem,nrow = 1,strip.position = "bottom")+
  labs(fill="",y="Turnover")
Fig_box_ASV_KO_turnover

####|####

#### Taxa Barplot ####
Fig_Summer <- A_D_rarW %>% tax_glom(.,taxrank ="Phylum") %>%
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
        lake_type=
          (A_D_rarW@sam_data$lake_type[match(.$lake_month,
                                             A_D_rarW@sam_data$lake_month)]),
        .) %>%
  dplyr::mutate(.,
                Phylum_legend=ifelse(median_abundance<1,"Phylum < 1%",Phylum)) %>%
  ggplot(.,aes(x=month, y=median_abundance,
               fill=fct_rev(fct_reorder(Phylum_legend,median_abundance))))+
  geom_col(width=0.8,color="black",size=0.2)+
  theme_bw()+theme(panel.grid = element_blank(),axis.title.y = element_text(size=12,face="bold"),
                   axis.title.x = element_blank(),axis.text.x=element_text(size=10),
                   axis.ticks =element_blank(),legend.position="right",
                   legend.title = element_text(face="bold",size=15),
                   legend.text = element_text(size=10),legend.key.size = unit(0.8,"cm"))+
  scale_x_discrete(expand = c(0.085,0.085))+
  scale_y_continuous(expand = c(0.02,0.02))+
  scale_fill_manual(values=palette_Charlotte)+
  guides(fill=guide_legend(ncol=1))+
  labs(fill="Microbial phyla",y="Relative abundance")+
  facet_wrap2(~ lake_order,scale="fixed",
              strip = strip_color_lake)

####__ overall plot ####
Fig_Summer<-Fig_Summer+theme(legend.position = "bottom",legend.direction = "vertical")+
  guides(fill=guide_legend(nrow=4))

Fig_BC_summer+Fig_Summer+
  #plot_layout(guides = 'collect')+
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face="bold",size=15))

####|####

#### Beta diversity ####

####__ Bray-Curtis ####
A_D_rarW.BC=ordinate(A_D_rarW, "PCoA","bray")

BC_AD_rarW<- plot_ordination(A_D_rarW,A_D_rarW.BC,shape="month",color="lake_order")

BC_AD_rarW$data$lake=factor(BC_AD_rarW$data$lake,levels=c("CHA","VER","GDP","BOI","CTL","CRJ1","CRJ2","VAI","JAB"))

####______ PCoA ####
Fig_AD_wout.BC_polygon<-ggplot(BC_AD_rarW$data,aes(Axis.1,Axis.2,color=lake,shape=month_full))+
  geom_convexhull(aes(color = lake,fill=lake,group = lake),alpha=0.1,show.legend = F,size=0.6)+
  geom_point(size=4.5) + theme_bw() +
  geom_label(BC_AD_centroids.df,mapping=aes(Axis.1,Axis.2,label = lake, fill = lake),
             color = "white",size = 7,fontface = "bold",show.legend = F,inherit.aes = F)+
  theme(aspect.ratio = 1,panel.grid = element_blank(),
        axis.title = element_text(size=20,face = "bold"), axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.title = element_text(size=25,face = "bold"),
        legend.text = element_text(size=20))+
  scale_color_manual(values=palette_lake)+
  scale_fill_manual(values=palette_lake)+
  scale_shape_manual(values = Summer_2021)+
  guides(color = guide_legend(override.aes=list(size=7)),
         shape = guide_legend(override.aes=list(size=7)))+
  labs(color="Lake",shape="Month - 2021",
       x=paste0("PCoA Axis 1 [",round(A_D_rarW.BC$values$Relative_eig[1],3)*100,"%]"),
       y=paste0("PCoA Axis 2 [",round(A_D_rarW.BC$values$Relative_eig[2],3)*100,"%]"))

Fig_AD_wout.BC_polygon

Fig_BC_summer<-Fig_AD_wout.BC_polygon+theme(legend.position = "bottom",legend.direction = "vertical")+
  guides(color = guide_legend(override.aes=list(size=6),nrow=3))

#### Without Cyanobacteria ####

A_D_rarW_noCyano.BC<-ordinate(A_D_rarW %>% subset_taxa(Phylum !="Cyanobacteria"), "PCoA","bray")

BC_AD_rarW_noCyano<- plot_ordination(A_D_rarW %>% subset_taxa(Phylum !="Cyanobacteria"),A_D_rarW_noCyano.BC,shape="month",color="lake_order")

BC_AD_rarW_noCyano$data$lake=factor(BC_AD_rarW_noCyano$data$lake,levels=c("CHA","VER","GDP","BOI","CTL","CRJ1","CRJ2","VAI","JAB"))

####______ disp ####

summer_2021.bd <- betadisper(A_D_rarW.BC.dist,group=ASV_factor$lake)
??betadisper
ASV_summer.df<-boxplot(summer_2021.bd, xlab="lake")

ASV_summer.df<- ASV_summer.df$stats %>% as.data.frame(.) 

colnames(ASV_summer.df)<-c("Champs","Verneuil","Grande-Paroisse","Bois","Créteil","Cergy-large","Cergy-small","Vaires","Jablines")
rownames(ASV_summer.df)<-c("d1","d2","d3","d4","d5")

ASV_summer.df<- ASV_summer.df %>% melt(.,value.name = "distance") %>%
  dplyr::mutate(.,
                lake=factor(variable,levels=rev(c("Verneuil","Champs",
                                          "Grande-Paroisse","Bois","Créteil","Cergy-large","Cergy-small",
                                          "Vaires","Jablines"))),
                Status=if_else(lake %in% c("Champs","Verneuil"),"Hypereutrophic",if_else(lake %in% c("Jablines","Vaires"),"Mesotrophic","Eutrophic")),
                Status=factor(Status,levels=c("Mesotrophic","Eutrophic","Hypereutrophic")),
                Features=paste("ASVs"),Status_features=paste0(Status,"_",Features),
                Status_features=factor(Status_features,levels=c("Mesotrophic_KOs","Mesotrophic_ASVs",
                                                                "Eutrophic_KOs","Eutrophic_ASVs",
                                                                "Hypereutrophic_KOs","Hypereutrophic_ASVs")))

summer_2021.bd <- betadisper(MDS_FUNC.hellinger,group=hellinger_factor$lake)

KO_summer.df<-boxplot(summer_2021.bd, xlab="lake")

KO_summer.df<- KO_summer.df$stats %>% as.data.frame(.) 

colnames(KO_summer.df)<-c("Champs","Verneuil","Grande-Paroisse","Bois","Créteil","Cergy-large","Cergy-small","Vaires","Jablines")
rownames(KO_summer.df)<-c("d1","d2","d3","d4","d5")

KO_summer.df<- KO_summer.df %>% melt(.,value.name = "distance") %>%
  dplyr::mutate(.,
                lake=factor(variable,levels=rev(c("Verneuil","Champs",
                                                  "Grande-Paroisse","Bois","Créteil","Cergy-large","Cergy-small",
                                                  "Vaires","Jablines"))),
                Status=if_else(lake %in% c("Champs","Verneuil"),"Hypereutrophic",if_else(lake %in% c("Jablines","Vaires"),"Mesotrophic","Eutrophic")),
                Status=factor(Status,levels=c("Mesotrophic","Eutrophic","Hypereutrophic")),
                Features=paste("KOs"),Status_features=paste0(Status,"_",Features),
                Status_features=factor(Status_features,levels=c("Mesotrophic_KOs","Mesotrophic_ASVs",
                                                                "Eutrophic_KOs","Eutrophic_ASVs",
                                                                "Hypereutrophic_KOs","Hypereutrophic_ASVs")))

disp_summer.df <- rbind(ASV_summer.df,KO_summer.df)

disp_summer.stats<- disp_summer.df %>%
  wilcox_test(distance ~ Status_features,p.adjust.method = "BH") %>%
  add_significance("p.adj")  %>% add_xy_position(x = "Status_features", dodge = 0.8) %>%
  mutate(g1=gsub("_","",as.character(.$group1)),g2=gsub("_","",as.character(.$group2))) %>%
  unite("comp",g1:g2,sep=" - ",remove = FALSE)

letters_disp_summer<-tibble(Status_features=unique(sort(disp_summer.df$Status_features)),
                                letters=cldList(p.adj ~ as.character(comp),disp_summer.stats,threshold  = 0.05)$Letter)

Fig_disp_summer <- disp_summer.df %>%
  dplyr::group_by(Status_features) %>%
  dplyr::mutate(.,median_distance=median(distance)) %>%
  cbind(.,letter=letters_disp_summer$letters[match(.$Status_features,letters_disp_summer$Status_features)]) %>%
  ggplot(.,aes(Status_features,distance,fill=Features))+
  geom_boxplot(width=0.6,color="black",show.legend = T)+theme_bw()+
  geom_point(size=1,color="black",alpha=0.5,show.legend = F)+
  geom_text(mapping=aes(y=1.02*max(distance),
                        label=letter),vjust=0,size=4.5,color="black",fontface="bold")+
  theme(panel.grid = element_blank(),plot.title = element_text(hjust=0.5,size=15,face = "bold"),
        axis.ticks = element_blank(),axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),axis.text.y = element_text(size=12),
        legend.title = element_text(size=15,face = "bold"),legend.text = element_text(size=12))+
  scale_fill_manual(values = c("white","lightgrey"))+
  scale_y_continuous(limits =c(0,1.025*max(disp_summer.df$distance)),position = "right",
                     labels = scales::label_number(accuracy = 0.1))+
  labs(fill="Dispersion")+
  facet_wrap2(~ Status,scale="free_x",nrow = 1,strip.position = "top",
              strip = strip_color_lake_type_afem<- strip_themed(
                background_x = elem_list_rect(fill = c("#2F6B9D","#E44600","#28A448"),
                                              color="black"),
                text_x = elem_list_text(colour = "white",
                                        face = "bold",size=12)))
Fig_disp_summer

####|####

#### SIMPER ####

####__ 0.1% Filter ####

total_depth_summer <- sum(taxa_sums(A_D_rarW))
threshold_summer <- 0.001 * total_depth_summer #0.1% abundance
abundant.taxa_summer <-(taxa_sums(A_D_rarW) > threshold_summer)
View(abundant.taxa_summer)

# tax table simper
tax_simper<- A_D_rarW@tax_table %>% t(.) %>% as.data.frame(.) %>% t(.)
tax_simper <- tax_simper %>% as.data.frame()

ps_summer_01_total<-prune_taxa(abundant.taxa_summer,A_D_rarW_freq)

# otu table simper
simper_01_total_com.df <- ps_summer_01_total@otu_table %>% as.data.frame() %>% t()

#### BLR ####

simper_01_total.df<-simper(simper_01_total_com.df,A_D_rarW_freq@sam_data$BOI,permutations = 999,parallel = 5)

# Contraste BOI
simper_BOI<- summary(simper_01_total.df, ordered =T)[[1]] %>% as.data.frame() %>%
  subset(.,cumsum <= 0.7 & p <=0.001)
list_simper_BOI<-rownames(simper_BOI)
#View(simper_BOI %>% dplyr::mutate(.,Phylum=tax_simper$Phylum[match(rownames(.),rownames(tax_simper))]))


#### Trophic Status ####
simper_01_total.df<-simper(simper_01_total_com.df,A_D_rarW_freq@sam_data$lake_type,permutations = 999,parallel = 5)

####__ Contrast H _ M ####
simper_H_M <- summary(simper_01_total.df)[[3]] %>% as.data.frame() %>%
  subset(.,cumsum <= 0.7 & p <=0.001)

list_simper_H_M<-rownames(simper_H_M)

simper_per_taxa<- simper_H_M
percent_taxa_H_M<-simper_per_taxa$cumsum[1]*100

for (i in 2:nrow(simper_per_taxa)) {
  percent<-(simper_per_taxa$cumsum[i]*100-simper_per_taxa$cumsum[i-1]*100)
  percent_taxa_H_M<-c(percent_taxa_H_M,percent)
}
simper_per_taxa$percent<-percent_taxa_H_M

Fig_simper_H_M_taxa<-simper_per_taxa %>% rownames_to_column(var="ASV") %>%
  cbind(.,Phylum=tax_simper$Phylum[match(.$ASV,rownames(tax_simper))]) %>%
  group_by(Phylum) %>% dplyr::summarise(sum_var=round(sum(percent),1)) %>%
  arrange(sum_var) %>%
  add_row(.,Phylum = '< 70 % cum. var.', sum_var = 100-sum(.$sum_var)) %>%
  dplyr::mutate(.,
                ypos = cumsum(sum_var)- 0.5*sum_var,
                xpos= if_else(Phylum =="Proteobacteria",0.8,0.5)) %>%
  ggplot(., aes(x=0.5, y=sum_var, fill=fct_rev(fct_reorder(Phylum,sum_var))))+
  geom_bar(stat="identity", width=1, color="black",show.legend = F) +
  geom_text(aes(y = ypos,x = xpos,
                label = sum_var),fontface = "bold",color = "black", size=4) +
  theme_void()+
  theme(legend.title = element_text(size=15,face="bold"),
        legend.text = element_text(size=12),
        legend.position = "bottom",legend.direction = "vertical")+
  scale_fill_manual(values = c("gray38","#B3CDE3","#CCEBC5","#DECBE4","#FED9A6","#FFFFCC","#E5D8BD","#FDDAEC","#F2F2F2","white")) +
  guides(fill=guide_legend(override.aes=list(size=5,shape=15),nrow=3))+
  coord_polar("y", start=0) +
  labs(fill="Microbial Phyla")  # remove background, grid, numeric labels
Fig_simper_H_M_taxa

####__ Contrast E _ M ####

simper_E_M <- summary(simper_01_total.df, ordered =T)[[2]] %>% as.data.frame() %>%
  subset(.,cumsum <= 0.7 & p <=0.001)
list_simper_E_M<-rownames(simper_E_M)

simper_per_taxa<- simper_E_M
percent_taxa_E_M<-simper_per_taxa$cumsum[1]*100

for (i in 2:nrow(simper_per_taxa)) {
  percent<-(simper_per_taxa$cumsum[i]*100-simper_per_taxa$cumsum[i-1]*100)
  percent_taxa_E_M<-c(percent_taxa_E_M,percent)
}
simper_per_taxa$percent<-percent_taxa_E_M

Fig_simper_E_M_taxa<-simper_per_taxa %>% rownames_to_column(var="ASV") %>%
  cbind(.,Phylum=tax_simper$Phylum[match(.$ASV,rownames(tax_simper))]) %>%
  group_by(Phylum) %>% dplyr::summarise(sum_var=round(sum(percent),1)) %>%
  add_row(.,Phylum = '< 70 % cum. var.', sum_var = 100-sum(.$sum_var)) %>%
  arrange(sum_var) %>%
  dplyr::mutate(.,
                ypos = cumsum(sum_var)- 0.5*sum_var,
                xpos= if_else(Phylum =="Bacteroidota",0.8,0.5)) %>%
  ggplot(., aes(x=0.5, y=sum_var, fill=fct_rev(fct_reorder(Phylum,sum_var))))+
  geom_bar(stat="identity", width=1, color="black",show.legend = F) +
  geom_text(aes(y = ypos, x = xpos,
                label = sum_var), fontface = "bold",color = "black", size=4) +
  theme_void()+
  theme(legend.title = element_text(size=15,face="bold"),
        legend.text = element_text(size=12),
        legend.position = "bottom",legend.direction = "vertical")+
  scale_fill_manual(values = c("#CCEBC5","gray38","#B3CDE3","#DECBE4","#FED9A6","#FFFFCC","#E5D8BD","#FDDAEC","#F2F2F2","white")) +
  guides(fill=guide_legend(override.aes=list(size=5,shape=15),nrow=3))+
  coord_polar("y", start=0) +
  labs(fill="Microbial Phyla")  # remove background, grid, numeric labels
Fig_simper_E_M_taxa

####__ Contrast E _ H ####

simper_E_H <- summary(simper_01_total.df, ordered =T)[[1]] %>% as.data.frame() %>%
  subset(.,cumsum <= 0.7 & p <=0.001)
list_simper_E_H<-rownames(simper_E_H)

simper_per_taxa<- simper_E_H
percent_taxa_E_H<-simper_per_taxa$cumsum[1]*100

for (i in 2:nrow(simper_per_taxa)) {
  percent<-(simper_per_taxa$cumsum[i]*100-simper_per_taxa$cumsum[i-1]*100)
  percent_taxa_E_H<-c(percent_taxa_E_H,percent)
}
simper_per_taxa$percent<-percent_taxa_E_H

Fig_simper_E_H_taxa<-simper_per_taxa %>% rownames_to_column(var="ASV") %>%
  cbind(.,Phylum=tax_simper$Phylum[match(.$ASV,rownames(tax_simper))]) %>%
  group_by(Phylum) %>% dplyr::summarise(sum_var=round(sum(percent),1)) %>%
  add_row(.,Phylum = '< 70 % cum. var.', sum_var = 100-sum(.$sum_var)) %>%
  arrange(sum_var) %>%
  dplyr::mutate(.,
                ypos = cumsum(sum_var)- 0.5*sum_var,
                xpos = if_else(Phylum =="Verrucomicrobiota",0.9,
                               if_else(Phylum =="Chloroflexi",0.7,0.5))) %>%
  ggplot(., aes(x=0.5, y=sum_var, fill=fct_rev(fct_reorder(Phylum,sum_var))))+
  geom_bar(stat="identity", width=1, color="black",show.legend = F) +
  geom_text(aes(y = ypos,x = xpos,
                label = sum_var),fontface = "bold",color = "black", size=4) +
  theme_void()+
  theme(legend.title = element_text(size=15,face="bold"),
        legend.text = element_text(size=12),
        legend.position = "bottom",legend.direction = "vertical")+
  scale_fill_manual(values = c("#FFFFCC","gray38","#DECBE4","#CCEBC5","lightgray","#E5D8BD","#B3CDE3"))+
  guides(fill=guide_legend(override.aes=list(size=5,shape=15),nrow=2))+
  coord_polar("y", start=0)+
  labs(fill="Microbial Phyla")  # remove background, grid, numeric labels
Fig_simper_E_H_taxa

# simper taxa list cumsum 0.7 p 0.05
list_simper<-unique(c(list_simper_H_M, list_simper_E_M,list_simper_E_H)) %>% as.data.frame() %>%
  rename("taxa_simper"=".") %>%
  dplyr::mutate(.,Phylum=tax_simper$Phylum[match(.$taxa_simper,rownames(tax_simper))],
                Class=tax_simper$Class[match(.$taxa_simper,rownames(tax_simper))],
                Family=tax_simper$Family[match(.$taxa_simper,rownames(tax_simper))],
                Genus=tax_simper$Genus[match(.$taxa_simper,rownames(tax_simper))],
                ASV_name=ASV_summer_simper$ASV_names[match(.$taxa_simper,ASV_summer_simper$original_names)])

#View(list_simper)
# new ps
View(A_D_rarW_freq@sam_data)
A_D_rarW_freq@sam_data <- A_D_rarW_freq@sam_data %>% as.matrix(.) %>% as.data.frame(.) %>%
  dplyr::mutate(.,
                lake_type=if_else(lake %in% c("CHA","VER"),"Hypereutrophic",
                                  if_else(lake %in% c("JAB","VAI"),"Mesotrophic",
                                          "Eutrophic")),
                lake_type=factor(lake_type,levels=c("Hypereutrophic","Eutrophic","Mesotrophic")),
                month_lake_type=paste0(month,"_",lake_type),
                lake_order=
                  if_else(lake_order=="CHA","Champs",
                          if_else(lake_order=="VER","Verneuil",
                                  if_else(lake_order=="GDP","Grande-Paroisse",
                                          if_else(lake_order=="BOI","Bois",
                                                  if_else(lake_order=="CTL","Créteil",
                                                          if_else(lake_order=="CRJ1","Cergy-small",
                                                                  if_else(lake_order=="CRJ2","Cergy-large",
                                                                          if_else(lake_order=="VAI","Vaires","Jablines")))))))),
                lake_order=factor(lake_order,levels=rev(c("Verneuil","Champs",
                                              "Grande-Paroisse","Bois","Créteil","Cergy-small","Cergy-large",
                                              "Vaires","Jablines"))),
                BOI=if_else(lake_order=="Bois","Bois","Others")) %>% sample_data(.)

ASV_summer_simper <- cbind(original_names=rownames(A_D_rarW_freq@otu_table),
                    ASV_names=paste0("ASV_",1:ntaxa(A_D_rarW_freq))) %>% as.data.frame()

simper_summer<-prune_taxa(list_simper$taxa_simper, A_D_rarW_freq)

length(list_simper$taxa_simper)

simper_summer.df<- psmelt(simper_summer) %>% .[,c(1:3,16,11,14,18)] %>% group_by(OTU,lake_order) %>%
  dplyr::summarise_at(vars("Abundance"),list(median=median)) %>% as.data.frame() %>%
  dplyr::mutate(.,Phylum=tax_simper$Phylum[match(.$OTU,rownames(tax_simper))],
                Class=tax_simper$Class[match(.$OTU,rownames(tax_simper))],
                Family=tax_simper$Family[match(.$OTU,rownames(tax_simper))],
                Genus=tax_simper$Genus[match(.$OTU,rownames(tax_simper))],
                ASV_name=ASV_summer_simper$ASV_names[match(.$OTU,ASV_summer_simper$original_names)],
                lake_type=simper_summer@sam_data$lake_type[match(.$lake_order,simper_summer@sam_data$lake_order)])

#View(simper_summer.df)
colnames(simper_summer.df)
Pheatmap_ASV.df<- simper_summer.df %>% .[c(8,3,2)] %>%
  pivot_wider(names_from = lake_order,values_from =median ,values_fill = 0) %>%
  column_to_rownames(var = "ASV_name") %>% as.matrix(.)

#View(Fig_simper_ASV_gg.df)

Pheatmap_ASV_norm <- t(apply(Pheatmap_ASV.df, 1, cal_z_score))

Fig_simper_ASV_gg.df<-Pheatmap_ASV_norm %>% melt(value.name = "Zscore") %>%
  rename("ASV"="Var1","lake"="Var2") %>%
  dplyr::mutate(.,
                Phylum=list_simper$Phylum[match(.$ASV,list_simper$ASV_name)],
                Class=list_simper$Class[match(.$ASV,list_simper$ASV_name)],
                Family=list_simper$Family[match(.$ASV,list_simper$ASV_name)],
                Genus=list_simper$Genus[match(.$ASV,list_simper$ASV_name)],
                ASV=factor(ASV,levels=rev(c("ASV_3468","ASV_689","ASV_635","ASV_670","ASV_681","ASV_739","ASV_2051","ASV_4170","ASV_4848","ASV_840",
                                            "ASV_5313","ASV_5291","ASV_726","ASV_835","ASV_2304","ASV_645","ASV_5434","ASV_2083","ASV_5138","ASV_660",
                                            "ASV_2012","ASV_1726","ASV_679","ASV_4929","ASV_1949","ASV_1650","ASV_1703","ASV_5501","ASV_1904","ASV_2040",
                                            "ASV_5310","ASV_1955"))),
                ASV_full=paste0(ASV,"_",Phylum,"_",Family,"_",Genus))

####__ heatmap ####
Fig_simper_ASV_gg<- Fig_simper_ASV_gg.df %>%
  ggplot(., aes(lake,ASV_full, fill= Zscore))+ 
  geom_tile(color="black",width=1,size=0.1)+
  theme_classic()+
  theme(panel.grid = element_blank(),aspect.ratio = 1.6,
        panel.spacing = unit(0, "cm"),
        axis.line=element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text.y = element_text(hjust = 0,size=8),
        axis.text.x = element_blank(),
        legend.text = element_text(size=12),
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.title = element_text(size=12,face="bold",vjust=1))+
  scale_y_discrete(position = "right")+
  scale_fill_gradient2(high="#8302DC", mid="#ECAAD8",low = "#F9F871")+
  labs(fill="Z-score")+
  coord_cartesian(xlim=c(0.5,9.31),ylim=c(0.5,33))+
  #row background (Phyla)
  annotate("rect",xmin=9.6,xmax=9.9,ymin=21.48,ymax=32.5,color="black",fill="#CCEBC5",size=0.5)+
  annotate("rect",xmin=9.6,xmax=9.9,ymin=19.48,ymax=21.48,color="black",fill="#F2F2F2",size=0.5)+
  annotate("rect",xmin=9.6,xmax=9.9,ymin=12.48,ymax=19.48,color="black",fill="#DECBE4",size=0.5)+
  annotate("rect",xmin=9.6,xmax=9.9,ymin=11.48,ymax=12.48,color="black",fill="#FED9A6",size=0.5)+
  annotate("rect",xmin=9.6,xmax=9.9,ymin=10.48,ymax=11.48,color="black",fill="#E5D8BD",size=0.5)+
  annotate("rect",xmin=9.6,xmax=9.9,ymin=3.48,ymax=10.48,color="black",fill="#B3CDE3",size=0.5)+
  annotate("rect",xmin=9.6,xmax=9.9,ymin=0.52,ymax=3.48,color="black",fill="#FFFFCC",size=0.5)+
  #row line 
  annotate("segment",x=0.5,xend=9.5,y=32.5,yend=32.5,size=0.5,color="black")+
  annotate("rect",xmin=7.5,xmax=9.5,ymin=29.5,ymax=32.5,size=0.8,color="black",fill=NA)+
  annotate("rect",xmin=0.5,xmax=2.5,ymin=21.48,ymax=29.5,size=0.8,color="black",fill=NA)+
  annotate("rect",xmin=0.5,xmax=2.5,ymin=14.48,ymax=19.48,size=0.8,color="black",fill=NA)+
  annotate("rect",xmin=7.5,xmax=9.5,ymin=12.48,ymax=14.48,size=0.8,color="black",fill=NA)+
  annotate("rect",xmin=0.5,xmax=2.5,ymin=11.48,ymax=12.48,size=0.8,color="black",fill=NA)+
  annotate("rect",xmin=7.5,xmax=9.5,ymin=9.48,ymax=11.48,size=0.8,color="black",fill=NA)+
  annotate("rect",xmin=7.5,xmax=9.5,ymin=9.48,ymax=10.48,size=0.8,color="black",fill=NA)+
  annotate("rect",xmin=0.5,xmax=2.5,ymin=3.48,ymax=9.48,size=0.8,color="black",fill=NA)+
  annotate("rect",xmin=7.5,xmax=9.5,ymin=0.5,ymax=3.48,size=0.8,color="black",fill=NA)+
  annotate("segment",x=0.5,xend=9.5,y=0.5,yend=0.5,size=0.5,color="black")+
  #col background (trophic status)
  annotate("rect",xmin=0.5,xmax=1.5,ymin=32.7,ymax=33.5,
           color="black",fill="#2F6B9D")+
  annotate("rect",xmin=1.5,xmax=2.5,ymin=32.7,ymax=33.5,
           color="black",fill="#5C9BD6")+
    annotate("rect",xmin=2.5,xmax=3.5,ymin=32.7,ymax=33.5,
             color="black",fill="#F2AB8C")+
  annotate("rect",xmin=3.5,xmax=4.5,ymin=32.7,ymax=33.5,
           color="black",fill="#EB7947")+
  annotate("rect",xmin=4.5,xmax=5.5,ymin=32.7,ymax=33.5,
           color="black",fill="#E45000")+
  annotate("rect",xmin=5.5,xmax=6.5,ymin=32.7,ymax=33.5,
            color="black",fill="#E40000")+
  annotate("rect",xmin=6.5,xmax=7.5,ymin=32.7,ymax=33.5,
            color="black",fill="#980000")+
  annotate("rect",xmin=7.5,xmax=8.5,ymin=32.7,ymax=33.5,
            color="black",fill="#1E7A36")+
  annotate("rect",xmin=8.5,xmax=9.5,ymin=32.7,ymax=33.5,
            color="black",fill="#28A448")+
  #col text (lake)
  annotate("text",x=1,y=33.1,label="JAB",color="white",fontface="bold",size=3)+
  annotate("text",x=2,y=33.1,label="VSM",color="white",fontface="bold",size=3)+
  annotate("text",x=3,y=33.1,label="CER-L",color="white",fontface="bold",size=3)+
  annotate("text",x=4,y=33.1,label="CER-S",color="white",fontface="bold",size=3)+
  annotate("text",x=5,y=33.1,label="CRE",color="white",fontface="bold",size=3)+
  annotate("text",x=6,y=33.1,label="BLR",color="white",fontface="bold",size=3)+
  annotate("text",x=7,y=33.1,label="LGP",color="white",fontface="bold",size=3)+
  annotate("text",x=8,y=33.1,label="CSM",color="white",fontface="bold",size=3)+
  annotate("text",x=9,y=33.1,label="VSS",color="white",fontface="bold",size=3)+
  #col line
  annotate("segment",x=0.5,xend=0.5,y=32.5,yend=0.5,size=0.5,color="black")+
  annotate("segment",x=9.5,xend=9.5,y=32.5,yend=0.5,size=0.5,color="black")
Fig_simper_ASV_gg

# perfome a PLS-DA analysis with status as factor
ps_summer<-A_D_rarW_freq
ASV_summer <- cbind(original_names=rownames(ps_summer@otu_table),
                    ASV_names=paste0("X",1:ntaxa(ps_summer))) %>% as.data.frame()
# update phyloseq rownames
rownames(ps_summer@otu_table)=
  ASV_summer$ASV_names[match(rownames(ps_summer@otu_table),ASV_summer$original_names)]

rownames(ps_summer@tax_table)=
  ASV_summer$ASV_names[match(rownames(ps_summer@tax_table),ASV_summer$original_names)]

summer.plsda <- ps_summer@otu_table %>% as.tibble(.) %>% t(.) %>%
  as.matrix(.) %>% mixOmics::plsda(.,ps_summer@sam_data$Status,scale = F)

# Create a data.frame to export MixOmics data to a clean ggplot
summer.plsda.df<-
  cbind(comp1=as.data.frame(summer.plsda$variates$X)$comp1,
        comp2=as.data.frame(summer.plsda$variates$X)$comp2) %>% as.data.frame(.) %>%
  dplyr::mutate(.,
                Status=summer.plsda$Y,Status=factor(Status,levels =c("Hypereutrophic","Eutrophic",
                                                                     "Mesotrophic")))

#### Other  ####
# PLS-DA #
Fig_summer.plsda <- ggplot(summer.plsda.df,aes(comp1,comp2,color = Status))+
  geom_point(size= 3)+ theme_bw()+
  theme(plot.title = element_text(hjust = 0.5,size = 12,face = "bold"),
        legend.title = element_text(size = 12,face= "bold"), panel.grid = element_blank(),
        axis.title = element_text(size = 12),axis.ticks = element_blank())+
  stat_ellipse(geom = "polygon",
               aes(color=lake,fill=lake),show.legend = F,alpha=0.3) +
  scale_color_manual(values= c("#2F6B9D","#E45000","#28A448"))+
  scale_fill_manual(values= c("#2F6B9D","#E45000","#28A448"))+
  labs(x=paste0("\nX - comp1 [",round(summer.plsda$prop_expl_var$X[[1]],3)*100,"%]"),
       y=paste0("X - comp2 [",round(summer.plsda$prop_expl_var$X[[2]],3)*100,"%]\n"))
Fig_summer.plsda

# export the n first loading variable
n_VIP=20
Fig_summer.plsda.loadings<-
  plotLoadings(summer.plsda,ndisplay = n_VIP,method = 'median',contrib = "max",
               title = paste0("PLS-DA \"Trophic status\" ",n_VIP,
                              " most disciminant ASVs"),
               size.title = 1.2,legend.color = c("#2F6B9D","#E45000","#28A448"))

tax_summer<-ps_summer@tax_table %>% as.data.frame()
View(ps_summer@tax_table)
View(Fig_summer.plsda.loadings)

# VIP table #
list_50_ASVs_summer<-
  cbind(PLS_DA_50=rownames(as.data.frame(Fig_summer.plsda.loadings)),
        importance=as.data.frame(Fig_summer.plsda.loadings)$importance,
        status_contrib=as.data.frame(Fig_summer.plsda.loadings)$GroupContrib) %>%
  as.data.frame() %>%
  dplyr::mutate(.,
                Phylum=tax_summer$Phylum[match(.$PLS_DA_50,rownames(tax_summer))],
                Class=tax_summer$Class[match(.$PLS_DA_50,rownames(tax_summer))],
                Genus=tax_summer$Genus[match(.$PLS_DA_50,rownames(tax_summer))])
View(list_50_ASVs_summer)

write.csv(list_50_ASVs_summer,
          "list_50_ASVs_summer.csv")

# Verrucomicrobiota #
Verrucomicrobiota_summer.df<-A_D_rarW_freq %>%
  tax_glom(taxrank = "Phylum") %>% subset_taxa(.,Phylum=="Verrucomicrobiota") %>%
  psmelt() %>% .[,c(1:18)]

Verrucomicrobiota_summer.stat<-Verrucomicrobiota_summer.df %>% wilcox_test(Abundance ~lake_type ,p.adjust.method = "BH") %>%
  add_significance("p.adj")  %>% add_xy_position(x = "lake_type", dodge = 0.8) %>%
  mutate(g1=gsub("_","",as.character(.$group1)),g2=gsub("_","",as.character(.$group2))) %>%
  unite("comp",g1:g2,sep=" - ",remove = FALSE)

letters_Verrucomicrobiota<-tibble(lake_type=unique(sort(Verrucomicrobiota_summer.df$lake_type)),
                                  letters=cldList(p.adj~as.character(comp),Verrucomicrobiota_summer.stat,threshold=0.05)$Letter)

Verrucomicrobiota_summer_barplot<-
  Verrucomicrobiota_summer.df %>%
  cbind(.,letter=letters_Verrucomicrobiota$letters[match(.$lake_type,letters_Verrucomicrobiota$lake_type)]) %>%
  ggplot(.,aes(x=lake_type,y=Abundance))+
  geom_boxplot(width=0.6,color="black",show.legend = F,fill="#E5D8BD")+theme_bw()+
  geom_text(mapping=aes(y=1.05*max(Abundance),label=letter),vjust=0,size=5.5,color="black",fontface="bold")+
  theme(panel.grid = element_blank(),axis.ticks = element_blank(),
        axis.title.x = element_blank(),axis.title.y = element_text(size=15,face="bold"),
        axis.text.y = element_text(size=12),axis.text.x = element_text(size=12),
        plot.title=element_text(size=15,hjust=0,face="bold"))+
  scale_x_discrete(expand = c(0.085,0.085))+
  labs(y="relative abundance",title="Verrucomicrobiota")
Verrucomicrobiota_summer_barplot

# Alphaproteobacteria #
Alphaproteobacteria_summer.df<-A_D_rarW_freq %>%
  tax_glom(taxrank = "Class") %>% subset_taxa(.,Class=="Alphaproteobacteria") %>%
  psmelt() %>% .[,c(1:18)]

Alphaproteobacteria_summer.stat<-Alphaproteobacteria_summer.df %>% wilcox_test(Abundance~lake_type ,p.adjust.method = "BH") %>%
  add_significance("p.adj")  %>% add_xy_position(x = "lake_type", dodge = 0.8) %>%
  mutate(g1=gsub("_","",as.character(.$group1)),g2=gsub("_","",as.character(.$group2))) %>%
  unite("comp",g1:g2,sep=" - ",remove = FALSE)

letters_Alphaproteobacteria<-tibble(lake_type=unique(sort(Alphaproteobacteria_summer.df$lake_type)),
                                    letters=cldList(p.adj~as.character(comp),Alphaproteobacteria_summer.stat,threshold=0.05)$Letter)

Alphaproteobacteria_summer_barplot<-
  Alphaproteobacteria_summer.df %>%
  cbind(.,letter=letters_Alphaproteobacteria$letters[match(.$lake_type,letters_Alphaproteobacteria$lake_type)]) %>%
  ggplot(.,aes(x=lake_type,y=Abundance))+
  geom_boxplot(width=0.6,color="black",show.legend = F,fill="#FED9A6")+theme_bw()+
  geom_text(mapping=aes(y=1.05*max(Abundance),label=letter),vjust=0,size=5.5,color="black",fontface="bold")+
  theme(panel.grid = element_blank(),axis.ticks = element_blank(),
        axis.title.x = element_blank(),axis.title.y = element_text(size=15,face="bold"),
        axis.text.y = element_text(size=12),axis.text.x = element_text(size=12),
        plot.title=element_text(size=15,hjust=0,face="bold"))+
  scale_x_discrete(expand = c(0.085,0.085))+
  labs(y="relative abundance",title="Alphaproteobacteria")
Alphaproteobacteria_summer_barplot

# Planctomycetota #
Planctomycetota_summer.df<-A_D_rarW_freq %>%
  tax_glom(taxrank = "Phylum") %>% subset_taxa(.,Phylum=="Planctomycetota") %>%
  psmelt() %>% .[,c(1:18)]

Planctomycetota_summer.stat<-Planctomycetota_summer.df %>% wilcox_test(Abundance ~lake_type ,p.adjust.method = "BH") %>%
  add_significance("p.adj")  %>% add_xy_position(x = "lake_type", dodge = 0.8) %>%
  mutate(g1=gsub("_","",as.character(.$group1)),g2=gsub("_","",as.character(.$group2))) %>%
  unite("comp",g1:g2,sep=" - ",remove = FALSE)

letters_Planctomycetota<-tibble(lake_type=unique(sort(Planctomycetota_summer.df$lake_type)),
                                letters=cldList(p.adj~as.character(comp),Planctomycetota_summer.stat,threshold=0.05)$Letter)

Planctomycetota_summer_barplot<-
  Planctomycetota_summer.df %>%
  cbind(.,letter=letters_Planctomycetota$letters[match(.$lake_type,letters_Planctomycetota$lake_type)]) %>%
  ggplot(.,aes(x=lake_type,y=Abundance))+
  geom_boxplot(width=0.6,color="black",show.legend = F,fill="#FFFFCC")+theme_bw()+
  geom_text(mapping=aes(y=1.05*max(Abundance),label=letter),vjust=0,size=5.5,color="black",fontface="bold")+
  theme(panel.grid = element_blank(),axis.ticks = element_blank(),
        axis.title.x = element_blank(),axis.title.y = element_text(size=15,face="bold"),
        axis.text.y = element_text(size=12),axis.text.x = element_text(size=12),
        plot.title=element_text(size=15,hjust=0,face="bold"))+
  scale_x_discrete(expand = c(0.085,0.085))+
  labs(y="relative abundance",title="Planctomycetota")
Planctomycetota_summer_barplot


Verrucomicrobiota_summer_barplot/Planctomycetota_summer_barplot

####|####

#### Spatial distance ####
BC_Geo_Status.df<- as.matrix(phyloseq::distance(A_D_rarW_freq,method = "bray")) %>%
  melt(.) %>% rename(G1=Var1,G2=Var2,BC=value) %>%
  filter(as.character(G1) != as.character(G2)) %>%
  mutate_if(is.factor,as.character) %>%
  cbind(.,
        G1_month_num=as.numeric(
          A_D_rarW_freq@sam_data$month_number[match(.$G1,
                                                    rownames(A_D_rarW_freq@sam_data))]),
        G2_month_num=as.numeric(
          A_D_rarW_freq@sam_data$month_number[match(.$G2,
                                                    rownames(A_D_rarW_freq@sam_data))]),
        G1_month=A_D_rarW_freq@sam_data$month_full[match(.$G1,
                                                         rownames(A_D_rarW_freq@sam_data))],
        G2_month=A_D_rarW_freq@sam_data$month_full[match(.$G2,
                                                         rownames(A_D_rarW_freq@sam_data))],
        G1_lake_month=A_D_rarW_freq@sam_data$lake_month[match(.$G1,
                                                              rownames(A_D_rarW_freq@sam_data))],
        G2_lake_month=A_D_rarW_freq@sam_data$lake_month[match(.$G2,
                                                              rownames(A_D_rarW_freq@sam_data))],
        G1_lake=A_D_rarW_freq@sam_data$lake[match(.$G1,
                                                  rownames(A_D_rarW_freq@sam_data))],
        G2_lake=A_D_rarW_freq@sam_data$lake[match(.$G2,
                                                  rownames(A_D_rarW_freq@sam_data))],
        G1_column=A_D_rarW_freq@sam_data$column[match(.$G1,
                                                      rownames(A_D_rarW_freq@sam_data))],
        G2_column=A_D_rarW_freq@sam_data$column[match(.$G2,
                                                      rownames(A_D_rarW_freq@sam_data))],
        
        Havers_dist=Havers_df$Haver_dist[match(.$G1,Havers_df$sample2)]) %>%
  subset(.,G1 != G2) %>%
  cbind(.,
        Chla1=param_table_AD$Chla_median[match(.$G1_lake_month,rownames(param_table_AD))],
        Chla2=param_table_AD$Chla_median[match(.$G2_lake_month,rownames(param_table_AD))]) %>%
  #subset(.,G1_lake == G2_lake) %>%
  dplyr::mutate(.,
                Status=sample_metadata$Status[match(.$G1_lake,sample_metadata$lake)],
                Chla_dist= rowMeans(.[,c(15,16)]),
                G1_lake=factor(G1_lake,levels=c("CHA","VER","GDP","BOI","CTL","CRJ1","CRJ2","VAI","JAB"))) %>%
  dplyr::filter(duplicated(BC) == FALSE)

spearman_BC_geo<-cor(BC_Geo_Status.df$BC,BC_Geo_Status.df$Havers_dist,method = "spearman")

spearman_BC_geo<-cor(BC_Geo_Status.df$BC,BC_Geo_Status.df$Havers_dist,method = "spearman")

lm_BC_status_geo <- lm(BC~Havers_dist,data=BC_Geo_Status.df)

var_BC_geo=paste0("Spearman ",round(spearman_BC_geo,2),", R2 ",round(summary(lm_BC_status_geo)[["adj.r.squared"]],2),", p <0.05")

Fig_BC_Geo<-BC_Geo_Status.df %>%
  ggplot(.,aes(x=Havers_dist,y=BC))+theme_classic()+
  geom_point(size=2,shape=21,fill="grey",color="black")+
  stat_smooth(method = "lm",show.legend = F,color="#4166F5",fill="lightblue",se = T)+
  theme(axis.ticks = element_blank(),
        axis.text = element_text(size=12),
        axis.title=element_text(size=12,face="bold"),aspect.ratio = 1,
        plot.caption = element_text(size=12,face="italic"))+
  #scale_color_manual(values=palette_lake)+
  labs(y="Taxonomical dissimilarity",color="Month",caption=var_BC_geo,
       x="Distance between lakes (m)")


spearman_BC_status<-cor(BC_Geo_Status.df$BC,BC_Geo_Status.df$Chla_dist,method = "spearman")

lm_BC_status_model <- lm(BC~Chla_dist,data=BC_Geo_Status.df)

var_BC_status=paste0("Spearman ",round(spearman_BC_status,2),", R2 ",round(summary(lm_BC_status_model)[["adj.r.squared"]],2),", p <0.01")

Fig_BC_Status<-BC_Geo_Status.df %>%
  ggplot(.,aes(x=Chla_dist,y=BC))+theme_classic()+
  geom_point(size=2,shape=21,fill="grey",color="black")+
  stat_smooth(method = "lm",show.legend = F,color="#4166F5",fill="lightblue",se = T)+
  theme(axis.ticks = element_blank(),
        axis.text = element_text(size=12),
        axis.title=element_text(size=12,face="bold"),aspect.ratio = 1,
        plot.caption = element_text(size=12,face="italic"))+
  scale_x_continuous(trans="log10",breaks = c(0,5,10,40,100,200))+
  #scale_color_manual(values=palette_lake)+
  scale_fill_manual(values=palette_month)+
  labs(y="Taxonomical dissimilarity",color="Month",caption = var_BC_status,
       x="Trophic gradient")

####____ overall plot ####
design_ST_Geo="
AB"

Fig_Hellinger_Geo+Fig_BC_Geo+
  plot_layout(design=design_ST_Geo)+
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size=15,face="bold"))

####|####

####____ SIMPER ####


lake_simper<- subset_samples(all_rarW_wout,lake=="CHA") %>% 
  simper_CHA.df <- vegan::simper(t(lake_simper@otu_table),group = lake_simper@sam_data$month)


View(as.data.frame(simper_CHA.df$A_F))

####|####

#### PERMANOVA ####

A_D_rarW.BC.dist=phyloseq::distance(A_D_rarW,method = "bray")
View(ASV_factor)
View(ASV_factor$Status)
####____ Status & Time ####
#merge the factor to the metabolites table
ASV_factor<-  data.frame(t(A_D_rarW_freq@otu_table)) %>%
  dplyr::mutate(.,
                lake=A_D_rarW_freq@sam_data$lake_order[match(rownames(.),A_D_rarW_freq@sam_data$samples)],
                Status=A_D_rarW_freq@sam_data$lake_type[match(rownames(.),A_D_rarW_freq@sam_data$samples)],
                month=A_D_rarW_freq@sam_data$month[match(rownames(.),A_D_rarW_freq@sam_data$samples)],
                LON=A_D_rarW_freq@sam_data$LON[match(rownames(.),A_D_rarW_freq@sam_data$samples)],
                LAT=A_D_rarW_freq@sam_data$LAT[match(rownames(.),A_D_rarW_freq@sam_data$samples)])
####______ permanova ####
permanova_ASV<-adonis2(A_D_rarW.BC.dist ~ Status + month + Status*month ,data=ASV_factor,permutations=999,method="bray")
permanova_ASV$`Pr(>F)`[[1]] #General P-value (if ANY sig. diff.)
permanova_ASV
####______ Betadisp ####
#Betadisp GENERAL (if NO sig., InterVar >> IntraVar)
betadisper_ASV<-as.data.frame(anova(betadisper(d=A_D_rarW.BC.dist,type="centroid",group=ASV_factor$Status))) #View(betadisper_ASV)
betadisper_ASV$`Pr(>F)`[[1]]
####______ Pairwise ####
# sorted Pairwise (adonis between group) 
pairwise_perm_ASV_sort<-pairwise.adonis(A_D_rarW.BC.dist,ASV_factor$Status) %>% .[order(.$pairs),] #View(pairwise_perm_ko_sort)
# sorted Pairwise_Betadisp IF NEDDED (betadisp between group)
pairwise_betadisper_ASV_sort<-permutest(betadisper(A_D_rarW.BC.dist,ASV_factor$Status), pairwise = TRUE) %>% .[[2]] %>% data.frame(.$permuted)
####______ df final ####
final_stats_ASV<-tibble(pairwise_perm_ASV_sort,betadisp=pairwise_betadisper_ASV_sort$permuted)
#View(final_stats_ASV)
