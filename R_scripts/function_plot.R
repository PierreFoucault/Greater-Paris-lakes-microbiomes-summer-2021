#### Import data ####

sample_metadata<-read_delim("/Users/pierre/Desktop/MetaG/50_10_bins_data/bin_metadata.csv",";", escape_double = FALSE,trim_ws = TRUE,show_col_types = FALSE) %>%
  dplyr::mutate(.,
                lake=factor(lake,levels=c("CHA","VER","GDP","BOI","CTL","CRJ1","CRJ2","VAI","JAB")),
                month_full=factor(month_full,levels=c("June","July","August","September","October","November")),
                Status=if_else(lake %in% c("VER","CHA"),"Hypereutrophic",
                       if_else(lake %in% c("VAI","JAB"),"Mesotrophic","Eutrophic")),
                Status=factor(Status,levels=c("Hypereutrophic","Eutrophic","Mesotrophic")),
                Distance=
                  if_else(lake %in% 
                            c("CRJ1","CRJ2","VER"),"NW",
                          if_else(lake %in% c("CTL","GDP","BOI"),"S",
                                  "NE")))
View(sample_metadata)
all_fun.df<- read_delim("all_fun.csv",",", escape_double = FALSE,trim_ws = TRUE,show_col_types = FALSE,num_threads = 4)%>%
  dplyr::mutate(.,lake=factor(lake,levels=c("CHA","VER","GDP","BOI","CTL","CRJ1","CRJ2","VAI","JAB"))) %>% .[,-1]

BGC_fun.df<- read_delim("BGC_fun.csv",",", escape_double = FALSE,trim_ws = TRUE,show_col_types = FALSE,num_threads=4) %>%
  dplyr::mutate(.,lake=factor(lake,levels=c("CHA","VER","GDP","BOI","CTL","CRJ1","CRJ2","VAI","JAB")),
                month_full=sample_metadata$month_full[match(.$sample,sample_metadata$sample)]) %>% .[,-1]

MAG_fun.df<- read_delim("MAG_func.csv",",", escape_double = FALSE,trim_ws = TRUE,show_col_types = FALSE,num_threads = 4)%>% .[,-1] %>%
  dplyr::mutate(.,sample_bin=paste0(.$sample,"_",.$bin),.after=2) %>%
  dplyr::mutate(.,
                lake=factor(lake,levels=c("CHA","VER","GDP","BOI","CTL","CRJ1","CRJ2","VAI","JAB")),
                month_full=sample_metadata$month_full[match(.$sample,sample_metadata$sample)],
                Domain=GTDB_MAG.df$Domain[match(.$sample_bin,GTDB_MAG.df$sample_bin)],
                Phylum=GTDB_MAG.df$Phylum[match(.$sample_bin,GTDB_MAG.df$sample_bin)],
                Class=GTDB_MAG.df$Class[match(.$sample_bin,GTDB_MAG.df$sample_bin)],
                Order=GTDB_MAG.df$Order[match(.$sample_bin,GTDB_MAG.df$sample_bin)],
                Family=GTDB_MAG.df$Family[match(.$sample_bin,GTDB_MAG.df$sample_bin)],
                Genus=GTDB_MAG.df$Genus[match(.$sample_bin,GTDB_MAG.df$sample_bin)]) %>%
  subset(.,KEGG_ko != "-") %>% group_by(KEGG_ko,sample,lake,month_full,bin,Domain,Phylum,Class,Order,Family,Genus,) %>%
  summarize(CPM = sum(CPM))

MAG_BGC.df<- read_delim("MAG_BGC.csv",",", escape_double = FALSE,trim_ws = TRUE,show_col_types = FALSE,num_threads=4) %>% .[,-1] %>%
  dplyr::mutate(.,sample_bin=paste0(.$sample,"_",.$bin)) %>%
  dplyr::mutate(.,
                lake=factor(lake,levels=c("CHA","VER","GDP","BOI","CTL","CRJ1","CRJ2","VAI","JAB")),
                month_full=sample_metadata$month_full[match(.$sample,sample_metadata$sample)],
                Domain=GTDB_MAG.df$Domain[match(.$sample_bin,GTDB_MAG.df$sample_bin)],
                Phylum=GTDB_MAG.df$Phylum[match(.$sample_bin,GTDB_MAG.df$sample_bin)],
                Class=GTDB_MAG.df$Class[match(.$sample_bin,GTDB_MAG.df$sample_bin)],
                Order=GTDB_MAG.df$Order[match(.$sample_bin,GTDB_MAG.df$sample_bin)],
                Family=GTDB_MAG.df$Family[match(.$sample_bin,GTDB_MAG.df$sample_bin)],
                Genus=GTDB_MAG.df$Genus[match(.$sample_bin,GTDB_MAG.df$sample_bin)]) %>%
  subset(.,KEGG_ko != "-") %>% group_by(KEGG_ko,sample,lake,month_full,bin,Domain,Phylum,Class,Order,Family,Genus,level_1,level_2,level_3,Vigneron) %>% summarize(CPM = sum(CPM))
#View(unique(MAG_BGC.df$Vigneron))

####. ####

#### All KO ####

####____ Dissimilarity ####

####______ vs ASV ####
#KO pairwise dissimilarity values
Hellinger_matrix_ML<- MDS_FUNC.hellinger %>% as.matrix() %>%
  as.data.frame() %>%  rownames_to_column(.,var="sample1_hellinger") %>%
  pivot_longer(-sample1_hellinger,names_to = "sample2_hellinger",values_to = "hellinger_dist") %>%
  dplyr::mutate(.,tmp1=sample1_hellinger,tmp2=sample2_hellinger,
                sample1_hellinger=sample_metadata$lake_month[match(tmp1,
                                                                   sample_metadata$sample)],
                sample2_hellinger=sample_metadata$lake_month[match(tmp2,
                                                                   sample_metadata$sample)],
                comp_hellinger=paste0(sample1_hellinger,"vs",sample2_hellinger))

#BC and pairwise dissimilarity table
BC_Hellinger.df<- A_D_rarW %>% subset_samples(.,column=="W2") %>%
  phyloseq::distance(.,method = "bray") %>% as.matrix() %>%
  as.data.frame() %>% rownames_to_column(.,var="sample1_BC") %>%
  pivot_longer(-sample1_BC,names_to = "sample2_BC",values_to = "BC_dist") %>%
  dplyr::mutate(.,
                sample1_BC = all_rarW_wout@sam_data$lake_month[match(.$sample1_BC,all_rarW_wout@sam_data$samples)],
                sample2_BC = all_rarW_wout@sam_data$lake_month[match(.$sample2_BC,all_rarW_wout@sam_data$samples)],
                comp_BC = paste0(sample1_BC,"vs",sample2_BC),
                lake1 = all_rarW_wout@sam_data$lake_order[match(.$sample1_BC,all_rarW_wout@sam_data$samples)],
                lake2 = all_rarW_wout@sam_data$lake_order[match(.$sample2_BC,all_rarW_wout@sam_data$samples)],
                comp_lake = paste0(lake1,"vs",lake2),
                trophicST1=if_else(lake1 %in% c("GDP","VER","CHA"),"Hypereutrophic",
                               if_else(lake1 %in% c("VAI","JAB"),"Oligotrophic",
                                       if_else(lake1 %in% c("CRJ1","CRJ2"),"Mesotrophic",
                                               "Eutrophic"))),
                trophicST2=if_else(lake2 %in% c("GDP","VER","CHA"),"Hypereutrophic",
                                   if_else(lake2 %in% c("VAI","JAB"),"Oligotrophic",
                                           if_else(lake2 %in% c("CRJ1","CRJ2"),"Mesotrophic",
                                                   "Eutrophic"))),
                comp_ST = paste0(trophicST1,"vs",trophicST2),
                month_number1 = all_rarW_wout@sam_data$month_number[match(.$sample1_BC,all_rarW_wout@sam_data$samples)],
                month_number2 = all_rarW_wout@sam_data$month_number[match(.$sample2_BC,all_rarW_wout@sam_data$samples)],
                comp_month = paste0(month_number1,"vs",month_number2),
                month_gap=abs(month_number1-month_number2)) %>%
  dplyr::mutate(.,
                hellinger_dist=Hellinger_matrix_ML$hellinger_dist[match(.$comp_BC,Hellinger_matrix_ML$comp_hellinger)]) %>%
  subset(.,sample1_BC != sample2_BC) %>% subset(.,lake1 == lake2) %>%
  dplyr::filter(duplicated(BC_dist) == FALSE) %>%
  dplyr::mutate(.,
                comp_ST=if_else(comp_ST=="HypereutrophicvsHypereutrophic","Hypereutrophic",
                                   if_else(comp_ST== "EutrophicvsEutrophic","Eutrophic",
                                           if_else(comp_ST=="OligotrophicvsOligotrophic","Oligotrophic",
                                                   "Mesotrophic"))),
                comp_ST=factor(comp_ST,levels=c("Hypereutrophic","Eutrophic",
                                                "Mesotrophic",
                                                "Oligotrophic")),
                lake1=factor(lake1,levels=c("CHA","VER","GDP","BOI","CTL","CRJ1","CRJ2","VAI","JAB")),
                comp_lake=factor(comp_lake,levels=c("CHAvsCHA","VERvsVER","GDPvsGDP",
                                                    "BOIvsBOI","CTLvsCTL",
                                                    "CRJ1vsCRJ1","CRJ2vsCRJ2",
                                                    "VAIvsVAI","JABvsJAB")))
               
#View(BC_Hellinger.df)

spearman_BC<-cor(BC_Hellinger.df$BC_dist,BC_Hellinger.df$hellinger_dist,method = "spearman")

lm_BC_model <- lm(hellinger_dist~BC_dist,data=BC_Hellinger.df)

var_BC=paste0("p<0.01, R2 ",round(summary(lm_BC_model)[["adj.r.squared"]],2))

BC_fun_centroids.df <- BC_Hellinger.df %>% dplyr::group_by(lake1) %>% dplyr::summarise_at(vars(3,14),median) %>%
  cbind(.,comp_ST=BC_Hellinger.df$comp_ST[match(.$lake1,BC_Hellinger.df$lake1)])



Fig_BC_Hellinger<- BC_Hellinger.df %>%
  ggplot(.,aes(BC_dist,hellinger_dist,group=1))+theme_classic()+
  geom_point(size=2,shape=21,fill="grey",color="black")+
  stat_smooth(method = "lm",color="#4166F5",fill="lightblue",show.legend = F)+
  theme_classic()+
  theme(axis.ticks = element_blank(),
        axis.text = element_text(size=10),
        axis.title=element_text(size=12,face="bold"),
        plot.caption = element_text(size=12,face="italic"))+
  scale_y_continuous(#limits = c(min(BC_Hellinger.df$hellinger_dist)/1.01,max(BC_Hellinger.df$hellinger_dist)*1.01),
                     breaks = c(0,0.06,0.08,0.1,0.12,0.14,0.16,0.18,0.2))+
  scale_x_continuous(#limits = c(min(BC_Hellinger.df$BC_dist)/1.01,max(BC_Hellinger.df$BC_dist)*1.01),
                     breaks = c(0.30,0.4,0.5,0.6,0.7,0.8,0.9,1))+
  labs(x="Taxonomical dissimilarity",y="Gene-content dissimilarity",caption = var_BC)
Fig_BC_Hellinger


design_fun_tax="
AAAAB
AAAAB
AAAAB
AAAAB
AAAAB
CCCC#"

Panel_redun_BC<-Fig_BC_Hellinger+Fig_box_Hellinger+Fig_box_BC+
  plot_layout(design = design_fun_tax)+
  plot_annotation(tag_levels = 'A') & theme(plot.tag=element_text(size=20,face="bold",color="black"))
Panel_redun_BC


#WU and pairwise dissimilarity table

WU_Hellinger.df<- A_D_rarW %>% subset_samples(.,column=="W2") %>%
  phyloseq::distance(.,method = "wunifrac") %>% as.matrix() %>%
  as.data.frame() %>% rownames_to_column(.,var="sample1_BC") %>%
  pivot_longer(-sample1_BC,names_to = "sample2_BC",values_to = "BC_dist") %>%
  dplyr::mutate(.,
                sample1_BC = all_rarW_wout@sam_data$lake_month[match(.$sample1_BC,all_rarW_wout@sam_data$samples)],
                sample2_BC = all_rarW_wout@sam_data$lake_month[match(.$sample2_BC,all_rarW_wout@sam_data$samples)],
                comp_BC = paste0(sample1_BC,"vs",sample2_BC),
                lake1 = all_rarW_wout@sam_data$lake_order[match(.$sample1_BC,all_rarW_wout@sam_data$samples)],
                lake2 = all_rarW_wout@sam_data$lake_order[match(.$sample2_BC,all_rarW_wout@sam_data$samples)],
                comp_lake = paste0(lake1,"vs",lake2),
                trophicST1=if_else(lake1 %in% c("GDP","VER","CHA"),"Hypereutrophic",
                                   if_else(lake1 %in% c("VAI","JAB"),"Oligotrophic",
                                           if_else(lake1 %in% c("CRJ1","CRJ2"),"Mesotrophic",
                                                   "Eutrophic"))),
                trophicST2=if_else(lake2 %in% c("GDP","VER","CHA"),"Hypereutrophic",
                                   if_else(lake2 %in% c("VAI","JAB"),"Oligotrophic",
                                           if_else(lake2 %in% c("CRJ1","CRJ2"),"Mesotrophic",
                                                   "Eutrophic"))),
                comp_ST = paste0(trophicST1,"vs",trophicST2),
                month_number1 = all_rarW_wout@sam_data$month_number[match(.$sample1_BC,all_rarW_wout@sam_data$samples)],
                month_number2 = all_rarW_wout@sam_data$month_number[match(.$sample2_BC,all_rarW_wout@sam_data$samples)],
                comp_month = paste0(month_number1,"vs",month_number2),
                month_gap=abs(month_number1-month_number2)) %>%
  dplyr::mutate(.,
                hellinger_dist=Hellinger_matrix_ML$hellinger_dist[match(.$comp_BC,Hellinger_matrix_ML$comp_hellinger)]) %>%
  subset(.,sample1_BC != sample2_BC) %>%  subset(.,lake1 == lake2) %>%
  dplyr::filter(duplicated(BC_dist) == FALSE) %>%
  dplyr::mutate(.,
                comp_ST=if_else(comp_ST=="HypereutrophicvsHypereutrophic","Hypereutrophic",
                                if_else(comp_ST== "EutrophicvsEutrophic","Eutrophic",
                                        if_else(comp_ST=="OligotrophicvsOligotrophic","Oligotrophic",
                                                "Mesotrophic"))),
                comp_ST=factor(comp_ST,levels=c("Hypereutrophic","Eutrophic","Mesotrophic","Oligotrophic")),
                lake1=factor(lake1,levels=c("CHA","VER","GDP","BOI","CTL","CRJ1","CRJ2","VAI","JAB")),
                comp_lake=factor(comp_lake,levels=c("CHAvsCHA","VERvsVER","GDPvsGDP",
                                                    "BOIvsBOI","CTLvsCTL",
                                                    "CRJ1vsCRJ1","CRJ2vsCRJ2",
                                                    "VAIvsVAI","JABvsJAB")))

spearman_WU<-cor(WU_Hellinger.df$BC_dist,WU_Hellinger.df$hellinger_dist,method = "spearman")

lm_WU_model <- lm(hellinger_dist~BC_dist*comp_ST,data=WU_Hellinger.df)

var_WU=paste0("Spearman ",round(spearman_WU,2),", R2 ",round(summary(lm_WU_model)[["adj.r.squared"]],2),", p<0.01")

WU_fun_centroids.df <- WU_Hellinger.df %>% dplyr::group_by(lake1) %>% dplyr::summarise_at(vars(3,14),median) %>%
  cbind(.,comp_ST=WU_Hellinger.df$comp_ST[match(.$lake1,WU_Hellinger.df$lake1)])

Fig_WU_Hellinger<- WU_Hellinger.df %>%
  ggplot(.,aes(BC_dist,hellinger_dist,color=lake1,fill=lake1))+theme_classic()+
  geom_convexhull(alpha=0.1,show.legend = F,size=0.6)+
  geom_point(size=3,shape=16,show.legend = F)+
  geom_label(WU_fun_centroids.df,mapping=aes(BC_dist,hellinger_dist,label = lake1, fill = lake1),color = "white",size = 4.5,
             fontface = "bold",show.legend = F,inherit.aes = F)+
  scale_color_manual(values=palette_lake)+
  scale_fill_manual(values=palette_lake)+ 
  theme(axis.ticks = element_blank(),axis.text = element_text(size=12),
        axis.title=element_text(size=15,face="bold"),aspect.ratio = 1,
        plot.caption = element_text(size=12,face="italic",hjust = 1))+
  scale_y_continuous(limits = c(min(BC_Hellinger.df$hellinger_dist)/1.01,max(BC_Hellinger.df$hellinger_dist)*1.01),
                     breaks = c(0.06,0.08,0.1,0.12,0.14,0.16))+
  scale_x_continuous(limits = c(min(WU_Hellinger.df$BC_dist)/1.01,max(WU_Hellinger.df$BC_dist)*1.01),
                     breaks = c(0.02,0.06,0.1,0.14,0.18))+
  labs(x="Taxonomical dissimilarity within lake",y="Functional dissimilarity within lake",
       caption=paste0("Method: Bray-Curtis & W.Unifrac\n",var_WU),
       color="Lake")
Fig_WU_Hellinger

Panel_redun_WU<-Fig_WU_Hellinger+Fig_lake_median_Chla+
  plot_layout(guides = 'collect',design = design_fun_tax)+plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size=15,face="bold"))
Panel_redun_WU

####____ boxplot ####
View(BC_Hellinger.df)
#Stats KO
BC_Hellinger.df <- BC_Hellinger.df %>%
  dplyr::mutate(.,
              comp_ST=if_else(lake1 %in% c("CHA","VER"),"Hypereutrophic",
                              if_else(lake1 %in% c("JAB","VAI"),"Mesotrophic",
                                      "Eutrophic")),
              comp_ST=factor(comp_ST,levels=c("Hypereutrophic","Eutrophic","Mesotrophic")))

box_Hellinger.stats<-BC_Hellinger.df %>%
  wilcox_test(hellinger_dist ~ comp_ST,p.adjust.method = "BH") %>%
  add_significance("p.adj")  %>% add_xy_position(x = "comp_ST", dodge = 0.8) %>%
  mutate(g1=gsub("_","",as.character(.$group1)),g2=gsub("_","",as.character(.$group2))) %>%
  unite("comp",g1:g2,sep=" - ",remove = FALSE)

letters_box_Hellinger<-tibble(comp_ST=unique(sort(BC_Hellinger.df$comp_ST)),
                              letters=cldList(p.adj ~ as.character(comp),box_Hellinger.stats,threshold  = 0.05)$Letter)

# boxplot
Fig_box_Hellinger<- BC_Hellinger.df %>% dplyr::group_by(comp_ST) %>%
  dplyr::mutate(.,mean_hellinger_dist=mean(hellinger_dist)) %>%
  cbind(.,letter=letters_box_Hellinger$letters[match(.$comp_ST,letters_box_Hellinger$comp_ST)]) %>%
  ggplot(.,aes(x=fct_reorder(comp_ST,mean_hellinger_dist),y=hellinger_dist,fill=comp_ST))+theme_bw()+
  geom_boxplot(show.legend = F,width=0.6,color="black")+
  geom_text(mapping=aes(y=1.01*max(hellinger_dist),
                        label=letter),vjust=0,size=7,color="darkred",fontface="bold")+
  theme(panel.grid = element_blank(),plot.title = element_blank(),
        axis.ticks = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size=15),axis.text.x = element_blank(),
        legend.title = element_text(size=25,face = "bold"),legend.text = element_text(size=20))+
  scale_y_continuous(limits = c(min(BC_Hellinger.df$hellinger_dist)/1.01,max(BC_Hellinger.df$hellinger_dist)*1.01),
                     breaks = c(0.06,0.08,0.1,0.12,0.14,0.16),
                     name = NULL, sec.axis = sec_axis(~.))+
  guides(y = "none")+
  scale_fill_manual(values = palette_lake_type_afem)

#Stats BC
box_BC.stats<- BC_Hellinger.df %>% wilcox_test(BC_dist ~ comp_ST,p.adjust.method = "BH") %>%
  add_significance("p.adj")  %>% add_xy_position(x = "comp_ST", dodge = 0.8) %>%
  mutate(g1=gsub("_","",as.character(.$group1)),g2=gsub("_","",as.character(.$group2))) %>%
  unite("comp",g1:g2,sep=" - ",remove = FALSE)

letters_box_BC<-tibble(comp_ST=unique(sort(BC_Hellinger.df$comp_ST)),
                       letters=cldList(p.adj ~ as.character(comp),box_BC.stats,threshold  = 0.05)$Letter)
# boxplot
Fig_box_BC<- BC_Hellinger.df %>% dplyr::group_by(comp_ST) %>%
  dplyr::mutate(.,mean_BC_dist=mean(BC_dist)) %>%
  cbind(.,letter=letters_box_BC$letters[match(.$comp_ST,letters_box_BC$comp_ST)]) %>%
  ggplot(.,aes(x=fct_reorder(comp_ST,mean_BC_dist),y=BC_dist,fill=comp_ST))+theme_bw()+
  geom_boxplot(show.legend = F,width=0.6,color="black")+
  geom_text(mapping=aes(y=1.01*max(BC_dist),
                        label=letter),vjust=0,size=7,color="darkred",fontface="bold")+
  theme(panel.grid = element_blank(),plot.title = element_blank(),
        axis.ticks = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size=15),axis.text.y = element_blank(),
        legend.title = element_text(size=25,face = "bold"),legend.text = element_text(size=20))+
  scale_y_continuous(limits = c(min(BC_Hellinger.df$BC_dist)/1.01,max(BC_Hellinger.df$BC_dist)*1.01),
                    breaks = c(0.30,0.4,0.5,0.6,0.7,0.8,0.9))+
  scale_fill_manual(values=palette_lake_type_afem)+
  guides(fill = guide_legend(ncol=1))+
  labs(fill="Trophic Status")+coord_flip()
Fig_box_BC

####____ Richness ####
colnames(all_fun.df)
# Compute ko richness (number of unique KO by sample)
fun_richness.df <- all_fun.df %>%  group_by(sample) %>%
  dplyr::summarise(ko_richness=length(unique(KEGG_ko))) %>%
  dplyr::mutate(.,
                month = sample_metadata$month[match(sample,sample_metadata$sample)],
                month_full = sample_metadata$month_full[match(sample,sample_metadata$sample)],
                lake_month = sample_metadata$lake_month[match(sample,sample_metadata$sample)],
                Status=sample_metadata$Status[match(sample,sample_metadata$sample)],
                lake = sample_metadata$lake[match(sample,sample_metadata$sample)],
                lake = factor(lake,levels = c("CHA","VER","GDP","BOI","CTL","CRJ1","CRJ2","VAI","JAB")),
                Status = factor(Status,levels = c("Hypereutrophic","Eutrophic","Mesotrophic","Oligotrophic")))

compare_richness<-as.data.frame(Fig_KO_ASV_richness$data)

#View(fun_richness.df)

spearman_S<-cor(compare_richness$ko_richness,compare_richness$ASV_richness,method = "spearman")

lm_S_model <- lm(ko_richness~ASV_richness*lake,data=compare_richness)

var_S=paste0("Spearman=",round(spearman_S,2),", R2=",round(summary(lm_S_model)[["adj.r.squared"]],2),", p<0.01")

S_fun_centroids.df <- compare_richness %>% dplyr::group_by(lake) %>% dplyr::summarise_at(vars(2,3),median)

Fig_KO_ASV_richness<- alphaDiv_rarW %>% subset(.,samples !="VAI-B-W1-ADN1") %>%
  subset(.,month %in% c("A","B","C","D")) %>% group_by(lake_month) %>%
  dplyr::summarise_at(.,.var="richness",.funs = median) %>% rename("ASV_richness"="richness") %>%
  dplyr::mutate(.,
                ko_richness=fun_richness.df$ko_richness[match(.$lake_month,fun_richness.df$lake_month)],
                lake=fun_richness.df$lake[match(.$lake_month,fun_richness.df$lake_month)],
                month_full=fun_richness.df$month_full[match(.$lake_month,fun_richness.df$lake_month)],
                Status=fun_richness.df$Status[match(.$lake_month,fun_richness.df$lake_month)]) %>%
  ggplot(.,aes(ASV_richness,ko_richness,fill=lake,color=lake))+theme_classic()+
  geom_convexhull(alpha=0.1,show.legend = F,size=0.6)+
  geom_point(size=3,shape=16,show.legend = F)+
  geom_label(S_fun_centroids.df,mapping=aes(ASV_richness,ko_richness,label = lake, fill = lake),color = "white",size = 4.5,
             fontface = "bold",show.legend = F,inherit.aes = F)+
  geom_point(size=2,shape=16,alpha=0.5,show.legend = F)+
  scale_color_manual(values=palette_lake)+
  scale_fill_manual(values=palette_lake)+
  theme(axis.ticks = element_blank(),axis.text = element_text(size=12),
        axis.title=element_text(size=15,face="bold"),aspect.ratio = 1,
        plot.caption = element_text(size=12,face="italic",hjust = 1))+
  expand_limits(y = c(min(fun_richness.df$ko_richness),max(fun_richness.df$ko_richness)))+
  labs(x="Prokaryotic ASV richness",y="Prokaryotic KO richness",
       caption=var_S)

Fig_KO_ASV_richness

Fig_KO_ASV_richness$data %>%
  dplyr::summarise_at(vars("ko_richness"),list(median=median, sd=sd))

#Stats KO
compare_richness_KO.stats<-compare_richness %>% wilcox_test(ko_richness ~ lake,p.adjust.method = "BH") %>%
  add_significance("p.adj")  %>% add_xy_position(x = "lake", dodge = 0.8) %>%
  mutate(g1=gsub("_","",as.character(.$group1)),g2=gsub("_","",as.character(.$group2))) %>%
  unite("comp",g1:g2,sep=" - ",remove = FALSE)

letters_box_KO<-tibble(lake=unique(sort(compare_richness$lake)),
                        letters=cldList(p.adj ~ as.character(comp),compare_richness_KO.stats,threshold  = 0.05)$Letter)

Fig_box_KO_S<- compare_richness %>%
  dplyr::group_by(lake) %>% dplyr::mutate(.,mean_S=mean(ko_richness)) %>%
  cbind(.,letter=letters_box_KO$letters[match(.$lake,letters_box_KO$lake)]) %>%
  ggplot(.,aes(x=fct_reorder(lake,mean_S),y=ko_richness,fill=lake))+theme_bw()+
  geom_boxplot(show.legend = T,width=0.6,color="black")+
  geom_text(mapping=aes(y=1.01*max(ko_richness),
                        label=letter),vjust=0,size=4,color="darkred",fontface="bold")+
  theme(panel.grid = element_blank(),plot.title = element_text(hjust=0.5,size=15,face = "bold"),
        axis.ticks = element_blank(),axis.title.x = element_blank(),
        axis.title.y = element_text(hjust=0.5,size=12,face = "bold"),
        axis.text.x = element_blank(),axis.text.y = element_text(size=10),
        plot.caption = element_text(size=12,face="italic",hjust = 1),
        legend.title = element_text(size=15,face = "bold"),legend.text = element_text(size=12))+
  scale_fill_manual(values = palette_lake)+
  labs(y="KOs richness")
Fig_box_KO_S

Fig_box_KO_S$data %>%
  dplyr::summarise_at(vars("ko_richness"),list(median=median, sd=sd))

#Stats ASV
compare_richness.stats<-compare_richness %>% wilcox_test(ASV_richness ~ Status,p.adjust.method = "BH") %>%
  add_significance("p.adj")  %>% add_xy_position(x = "Status", dodge = 0.8) %>%
  mutate(g1=gsub("_","",as.character(.$group1)),g2=gsub("_","",as.character(.$group2))) %>%
  unite("comp",g1:g2,sep=" - ",remove = FALSE)

letters_box_ASV<-tibble(Status=unique(sort(compare_richness$Status)),
                              letters=cldList(p.adj ~ as.character(comp),compare_richness.stats,threshold  = 0.05)$Letter)

Fig_box_ASV_S<- compare_richness %>%
  dplyr::group_by(Status) %>% dplyr::mutate(.,mean_S=mean(ASV_richness)) %>%
  cbind(.,letter=letters_box_ASV$letters[match(.$Status,letters_box_ASV$Status)]) %>%
  ggplot(.,aes(x=fct_rev(fct_reorder(Status,mean_S)),y=ASV_richness,fill=Status))+theme_bw()+
  geom_boxplot(show.legend = T,width=0.6,color="black")+
  geom_text(mapping=aes(y=1.01*max(ASV_richness),
                        label=letter),vjust=0,size=4,color="darkred",fontface="bold")+
  theme(panel.grid = element_blank(),plot.title = element_text(hjust=0.5,size=15,face = "bold"),
        axis.ticks = element_blank(),axis.title.x = element_blank(),
        axis.title.y = element_text(hjust=0.5,size=12,face = "bold"),
        axis.text.y = element_blank(),axis.text.x = element_text(size=10),
        plot.caption = element_text(size=12,face="italic",hjust = 1),
        legend.title = element_text(size=15,face = "bold"),legend.text = element_text(size=12))+
  scale_fill_manual(values = palette_status)+
  labs(x="ASVs richness")+coord_flip()
Fig_box_ASV_S


design_redun_richness="
AAAC
AAAC
AAAC
AAAC
BBB#"

Fig_KO_ASV_richness+Fig_box_ASV_S+Fig_box_KO_S+
  plot_layout(design=design_redun_richness,guides='collect')+
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face="bold",size=15))

####____ Evenness ####
colnames(FUNC.hellinger)
fun_alphadiv.df <- data.frame(sample=rownames(FUNC.hellinger),
                              chao1=as.vector(estimateR(FUNC.hellinger)),
                              shannon=as.vector(vegan::diversity(FUNC.hellinger)),
                              richness=as.vector(vegan::specnumber(FUNC.hellinger)),
                              inv_simpson=as.vector(vegan::diversity(FUNC.hellinger,index = "simpson"))) %>%
  dplyr::mutate(.,evenness=shannon/(log(richness)))

View(fun_alphadiv.df)

fun_alphadiv.df %>%
  dplyr::summarise_at(vars("shannon"),list(median=median, sd=sd))


list_K0<-colnames(MDS_FUNC.df)
FUNC.hellinger<- MDS_FUNC.df %>% labdsv::hellinger(.) %>% as.data.frame(.) %>% magrittr::set_colnames(list_K0)

####____ KO Turnover ####

Summer_Ko_codyn.df <- all_fun.df %>% subset(month %in% c("A","B","C","D")) %>% .[,c(1,2,9)] %>%
  group_by(KEGG_ko,sample) %>% dplyr::summarise(CPM=sum(CPM)) %>%
  cbind(.,
        month_full=sample_metadata$month_full[match(.$sample,sample_metadata$sample)],
        lake=sample_metadata$lake[match(.$sample,sample_metadata$sample)]) %>%
  cbind(.,month_number=A_D_rarW@sam_data$month_number[match(.$month_full,A_D_rarW@sam_data$month_full)]) %>%
  dplyr::mutate(.,
        sample_codyn=paste0(lake,"_",month_full,"_",KEGG_ko)) #View(Summer_Ko_codyn.df)

KO_turnover<-turnover(Summer_Ko_codyn.df,time.var = "month_number",species.var = "KEGG_ko",
                       abundance.var = "CPM",replicate.var = "lake",metric = "total") %>%
  dplyr::rename("KO_turnover"="total") %>%
  dplyr::mutate(.,lake_month=paste0(lake,"_",month_number),
                lake=factor(.$lake,levels=c("CHA","VER","GDP","BOI","CTL","CRJ1","CRJ2","VAI","JAB"))) #View(KO_turnover)

ASV_KO_turnover<-ASV_turnover %>%
  dplyr::mutate(.,KO_turnover=KO_turnover$KO_turnover[match(.$lake_month,KO_turnover$lake_month)],.after=1) #View(ASV_KO_turnover)

turnover_centroids.df <- ASV_KO_turnover %>% dplyr::group_by(lake) %>% dplyr::summarise_at(vars(1,2),median)

Fig_KO_ASV_turnover<- ASV_KO_turnover %>%
  ggplot(.,aes(ASV_turnover,KO_turnover,fill=lake,color=lake,group=lake))+theme_classic()+
  geom_convexhull(alpha=0.1,show.legend = F,size=0.6)+
  geom_point(size=3,shape=16,show.legend = F)+
  geom_label(turnover_centroids.df,mapping=aes(ASV_turnover,KO_turnover,label = lake, fill = lake),color = "white",size = 4.5,
             fontface = "bold",show.legend = F,inherit.aes = F)+
  geom_point(size=2,shape=16,alpha=0.5,show.legend = F)+
  scale_color_manual(values=palette_lake)+
  scale_fill_manual(values=palette_lake)+
  theme(axis.ticks = element_blank(),axis.text = element_text(size=12),
        axis.title=element_text(size=15,face="bold"),aspect.ratio = 1,
        plot.caption = element_text(size=12,face="italic",hjust = 1))+
  labs(x="Prokaryotic ASV turnover",y="Prokaryotic KO turnover")
Fig_KO_ASV_turnover

#Stats KO
compare_turnover_KO.stats<-ASV_KO_turnover %>% wilcox_test(KO_turnover ~ Status,p.adjust.method = "BH") %>%
  add_significance("p.adj")  %>% add_xy_position(x = "Status", dodge = 0.8) %>%
  mutate(g1=gsub("_","",as.character(.$group1)),g2=gsub("_","",as.character(.$group2))) %>%
  unite("comp",g1:g2,sep=" - ",remove = FALSE)

letters_box_KO_turnover<-tibble(Status=unique(sort(ASV_KO_turnover$Status)),
                                 letters=cldList(p.adj ~ as.character(comp),compare_turnover_KO.stats,threshold  = 0.05)$Letter)

Fig_box_KO_turnover<- ASV_KO_turnover %>%
  dplyr::group_by(Status) %>% dplyr::mutate(.,mean_turnover=mean(KO_turnover)) %>%
  cbind(.,letter=letters_box_KO_turnover$letters[match(.$Status,letters_box_KO_turnover$Status)]) %>%
  ggplot(.,aes(x=fct_reorder(Status,mean_turnover),y=KO_turnover,fill=Status))+theme_bw()+
  geom_boxplot(show.legend = T,width=0.6,color="black")+
  geom_text(mapping=aes(y=1.01*max(KO_turnover),
                        label=letter),vjust=0,size=4,color="darkred",fontface="bold")+
  theme(panel.grid = element_blank(),plot.title = element_text(hjust=0.5,size=15,face = "bold"),
        axis.ticks = element_blank(),axis.title.x = element_blank(),
        axis.title.y = element_text(hjust=0.5,size=12,face = "bold"),
        axis.text.x = element_blank(),axis.text.y = element_text(size=10),
        plot.caption = element_text(size=12,face="italic",hjust = 1),
        legend.title = element_text(size=15,face = "bold"),legend.text = element_text(size=12))+
  scale_fill_manual(values = palette_status)+
  labs(y="KO turnover")
Fig_box_KO_turnover

#Stats ASV
ASV_KO_turnover
compare_turnover.stats<-ASV_KO_turnover %>% wilcox_test(ASV_turnover ~ Status,p.adjust.method = "BH") %>%
  add_significance("p.adj")  %>% add_xy_position(x = "Status", dodge = 0.8) %>%
  mutate(g1=gsub("_","",as.character(.$group1)),g2=gsub("_","",as.character(.$group2))) %>%
  unite("comp",g1:g2,sep=" - ",remove = FALSE)

letters_box_ASV_turnover<-tibble(Status=unique(sort(ASV_KO_turnover$Status)),
                        letters=cldList(p.adj ~ as.character(comp),compare_turnover.stats,threshold  = 0.05)$Letter)

Fig_box_ASV_turnover<- ASV_KO_turnover %>%
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


design_redun_turnover="
AAAC
AAAC
AAAC
AAAC
BBBD"

Fig_KO_ASV_turnover+Fig_box_ASV_turnover+Fig_box_KO_turnover+guide_area()+
  plot_layout(design=design_redun_turnover,guides='collect')+
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face="bold",size=15))

#### BGCs KO ####

####____ Dissimilarity ####
BGC_matrix_ML<- MDS_BGC.hellinger %>% as.matrix() %>%
  as.data.frame() %>%  rownames_to_column(.,var="sample1_hellinger") %>%
  pivot_longer(-sample1_hellinger,names_to = "sample2_hellinger",values_to = "hellinger_dist") %>%
  dplyr::mutate(.,tmp1=sample1_hellinger,tmp2=sample2_hellinger,
                sample1_hellinger=sample_metadata$lake_month[match(tmp1,
                                                                   sample_metadata$sample)],
                sample2_hellinger=sample_metadata$lake_month[match(tmp2,
                                                                   sample_metadata$sample)],
                comp_hellinger=paste0(sample1_hellinger,"vs",sample2_hellinger)) #

BC_BGC.df<- A_D_rarW %>% subset_samples(.,column=="W2") %>%
  phyloseq::distance(.,method = "wunifrac") %>% as.matrix() %>%
  as.data.frame() %>% rownames_to_column(.,var="sample1_BC") %>%
  pivot_longer(-sample1_BC,names_to = "sample2_BC",values_to = "BC_dist") %>%
  dplyr::mutate(.,
                sample1_BC = all_rarW_wout@sam_data$lake_month[match(.$sample1_BC,all_rarW_wout@sam_data$samples)],
                sample2_BC = all_rarW_wout@sam_data$lake_month[match(.$sample2_BC,all_rarW_wout@sam_data$samples)],
                comp_BC = paste0(sample1_BC,"vs",sample2_BC),
                lake1 = all_rarW_wout@sam_data$lake_order[match(.$sample1_BC,all_rarW_wout@sam_data$samples)],
                lake2 = all_rarW_wout@sam_data$lake_order[match(.$sample2_BC,all_rarW_wout@sam_data$samples)],
                comp_lake = paste0(lake1,"vs",lake2),
                trophicST1=if_else(lake1 %in% c("GDP","VER","CHA"),"Hypereutrophic",
                                   if_else(lake1 %in% c("VAI","JAB"),"Oligotrophic",
                                           if_else(lake1 %in% c("CRJ1","CRJ2"),"Mesotrophic",
                                                   "Eutrophic"))),
                trophicST2=if_else(lake2 %in% c("GDP","VER","CHA"),"Hypereutrophic",
                                   if_else(lake2 %in% c("VAI","JAB"),"Oligotrophic",
                                           if_else(lake2 %in% c("CRJ1","CRJ2"),"Mesotrophic",
                                                   "Eutrophic"))),
                comp_ST = paste0(trophicST1,"vs",trophicST2),
                month_number1 = all_rarW_wout@sam_data$month_number[match(.$sample1_BC,all_rarW_wout@sam_data$samples)],
                month_number2 = all_rarW_wout@sam_data$month_number[match(.$sample2_BC,all_rarW_wout@sam_data$samples)],
                comp_month = paste0(month_number1,"vs",month_number2),
                month_gap=abs(month_number1-month_number2)) %>%
  dplyr::mutate(.,
                hellinger_dist=BGC_matrix_ML$hellinger_dist[match(.$comp_BC,BGC_matrix_ML$comp_hellinger)]) %>%
  subset(.,sample1_BC != sample2_BC) %>%  subset(.,lake1 == lake2) %>%
  dplyr::filter(duplicated(BC_dist) == FALSE) %>%
  dplyr::mutate(.,
                comp_ST=if_else(comp_ST=="HypereutrophicvsHypereutrophic","Hypereutrophic",
                                if_else(comp_ST== "EutrophicvsEutrophic","Eutrophic",
                                        if_else(comp_ST=="OligotrophicvsOligotrophic","Oligotrophic",
                                                "Mesotrophic"))),
                comp_ST=factor(comp_ST,levels=c("Hypereutrophic","Eutrophic","Mesotrophic","Oligotrophic")),
                lake1=factor(lake1,levels=c("CHA","VER","GDP","BOI","CTL","CRJ1","CRJ2","VAI","JAB")),
                comp_lake=factor(comp_lake,levels=c("CHAvsCHA","VERvsVER","GDPvsGDP",
                                                    "BOIvsBOI","CTLvsCTL",
                                                    "CRJ1vsCRJ1","CRJ2vsCRJ2"))) #View(BC_BGC.df)

spearman_cor<-cor(BC_BGC.df$hellinger_dist,BC_BGC.df$BC_dist,method = "spearman")

model <- lm(BC_dist ~ hellinger_dist, data = BC_BGC.df)
summary(model)

mixed.lmer <- lmer(BC_dist ~ hellinger_dist + (1|comp_month), data = BC_BGC.df)
summary(mixed.lmer)


BC_BGC_centroids.df <- BC_BGC.df %>% dplyr::group_by(lake1) %>% dplyr::summarise_at(vars(3,14),median) %>%
  cbind(.,comp_ST=BC_BGC.df$comp_ST[match(.$lake1,BC_BGC.df$lake1)])

Fig_BC_BGC<-ggplot(BC_BGC.df,aes(BC_dist,hellinger_dist,fill=lake1,color=lake1))+
  theme_classic()+
  stat_smooth(method ='lm',mapping=aes(group=comp_ST,color=comp_ST,fill=comp_ST))+
  #geom_convexhull(alpha=0.1,show.legend = F,size=0.6)+
  #geom_point(size=3,shape=16,show.legend = F)+
  #geom_label(BC_BGC_centroids.df,mapping=aes(BC_dist,hellinger_dist,label = lake1, fill = lake1),color = "white",size = 4.5,
             #fontface = "bold",show.legend = F,inherit.aes = F)+
  #scale_color_manual(values=palette_lake)+
  #scale_fill_manual(values=palette_lake)+ 
  theme(axis.ticks = element_blank(),axis.text = element_text(size=12),
        axis.title=element_text(size=15,face="bold"),aspect.ratio = 1,
        plot.caption = element_text(size=12,face="italic",hjust = 1))+
  scale_y_continuous(labels = scales::label_number(accuracy = 0.01))+
  expand_limits(y = c(min(BC_BGC.df$hellinger_dist),max(BC_BGC.df$hellinger_dist)))+
  expand_limits(x = c(min(BC_BGC.df$BC_dist),max(BC_BGC.df$BC_dist)))+
  labs(x="Taxonomical dissimilarity",y="BGCs functions dissimilarity",title="\n\n",
       caption=paste0("Hellinger transformed, Method: Bray-Curtis and W.Unifrac\n"))
Fig_BC_BGC

####____ Richness ####
list_BGC<-colnames(MDS_BGC.df)
BGC.hellinger<- MDS_BGC.df %>% labdsv::hellinger(.) %>% as.data.frame(.) %>% magrittr::set_colnames(list_BGC)

BGC_richness.df <- BGC_fun.df %>%  group_by(sample) %>%
  dplyr::summarise(ko_richness=length(unique(KEGG_ko))) %>%
  dplyr::mutate(.,
                month = sample_metadata$month[match(sample,sample_metadata$sample)],
                month_full = sample_metadata$month_full[match(sample,sample_metadata$sample)],
                lake_month = sample_metadata$lake_month[match(sample,sample_metadata$sample)],
                lake = sample_metadata$lake[match(sample,sample_metadata$sample)],
                lake = factor(lake,levels = c("CHA","GDP","VER","CTL","BOI","CRJ1","CRJ2","JAB","VAI")))
#View(BGC_richness.df)
spearman_cor<-cor(Fig_BGC_ASV_richness$data$ASV_richness,Fig_BGC_ASV_richness$data$ko_richness,method = "spearman")

model <- lm(ASV_richness ~ ko_richness, data = Fig_BGC_ASV_richness$data)
summary(model)

mixed.lmer <- lmer(ASV_richness ~ ko_richness + (1|month_full), data = Fig_BGC_ASV_richness$data)
summary(mixed.lmer)

var_sampling_month=paste0("Spearman=",round(spearman_cor,2),", R2=",round(summary(model)[["adj.r.squared"]],2),", p<0.01, ",
                          round((summary(mixed.lmer)[["varcor"]][["month_full"]]/
                                   (summary(mixed.lmer)[["varcor"]][["month_full"]]+1234)*100),2),
                          "% of the variance explained by the sampling month factor")

Fig_BGC_ASV_richness<- alphaDiv_rarW %>% subset(.,samples !="VAI-B-W1-ADN1") %>%
  subset(.,month %in% c("A","B","C","D")) %>% group_by(lake_month) %>%
  dplyr::summarise_at(.,.var="richness",.funs = median) %>% rename("ASV_richness"="richness") %>%
  dplyr::mutate(.,
                ko_richness=BGC_richness.df$ko_richness[match(.$lake_month,BGC_richness.df$lake_month)],
                lake=BGC_richness.df$lake[match(.$lake_month,BGC_richness.df$lake_month)],
                month_full=BGC_richness.df$month_full[match(.$lake_month,BGC_richness.df$lake_month)],
                Status=sample_metadata$Status[match(.$lake_month,BGC_richness.df$lake_month)]) %>%
  ggplot(.,aes(ASV_richness,ko_richness,fill=Status,color=lake))+theme_classic()+
  geom_point(size=2,alpha=0.5,shape=21,show.legend = F)+
  stat_smooth(method = lm,show.legend = F,se = F)+
  theme(axis.ticks = element_blank(),axis.text = element_text(size=12),
        axis.title=element_text(size=12,face="bold"),aspect.ratio = 1,
        plot.caption = element_text(size=12,face="italic",hjust = 1))+
  scale_color_manual(values=palette_slope)+
  scale_fill_manual(values=palette_slope)+
  labs(x="Prokaryotic ASV richness",y="Prokaryotic  BGCs KO richness",
       caption=var_sampling_month)
Fig_BGC_ASV_richness

####__ Panel  ####
design_redun<-"
AB
CD"

Fig_BC_Hellinger+Fig_BC_WU_Hellinger+Fig_KO_ASV_richness+Fig_BC_BGC+
  plot_layout(design = design_redun)+plot_annotation(tag_levels = "A") &
  theme(aspect.ratio = NULL,plot.tag=element_text(face="bold",size=15))

####. ####

#### Distances all prokaryotes KOs  ####

####____ Status ####
#merge the factor to the metabolites table
jaccard_factor<-  data.frame(MDS_FUNC.df) %>%
  dplyr::mutate(.,Status=Ko_MDS_jaccard.df$Status[match(rownames(.),rownames(Ko_MDS_jaccard.df))]) #View(hellinger_factor)
####______ permanova ####
permanova_ko<-adonis2(MDS_FUNC.jaccard~ jaccard_factor$Status,data=jaccard_factor,permutations=999,)
permanova_ko$`Pr(>F)`[[1]] #General P-value (if ANY sig. diff.)
permanova_ko
####______ Betadisp ####
#Betadisp GENERAL (if NO sig., InterVar >> IntraVar)
betadisper_ko<-as.data.frame(anova(betadisper(d=MDS_FUNC.hellinger,type="centroid",group=hellinger_factor$Status))) #View(betadisper_ASV)
betadisper_ko$`Pr(>F)`[[1]]
####______ Pairwise ####
# sorted Pairwise (adonis between group) 
pairwise_perm_ko_sort<-pairwise.adonis(MDS_FUNC.jaccard,jaccard_factor$Status) %>% .[order(.$pairs),] #View(pairwise_perm_ko_sort)
# sorted Pairwise_Betadisp IF NEDDED (betadisp between group)
pairwise_betadisper_ko_sort<-permutest(betadisper(MDS_FUNC.hellinger,hellinger_factor$Status), pairwise = TRUE) %>% .[[2]] %>% data.frame(.$permuted)
####______ df final ####
final_stats_ko<-tibble(pairwise_perm_ko_sort,betadisp=pairwise_betadisper_ko_sort$permuted)#
View(final_stats_ko)


####__ Hellinger ####
# Hellinger abudance matrix with vegan
?set_colnames
# load matrices
library(magrittr)
View(MDS_FUNC.df)
colnames(MDS_FUNC.df)
list_K0<-colnames(MDS_FUNC.df)
MDS_FUNC.hellinger<- MDS_FUNC.df %>% labdsv::hellinger(.) %>% as.data.frame(.) %>% set_colnames(list_K0) %>%
  vegdist(., method="bray")

(MDS_FUNC.hellinger)
??hellinger
# Compute MDS
MDS_FUNC.hellinger.mds<-cmdscale(MDS_FUNC.hellinger,eig=TRUE, k=2)
library(magrittr)
# Store the axis coordinates with metadat for ggplot
Ko_MDS_hellinger.df <- MDS_FUNC.hellinger.mds$points %>% as.data.frame(.) %>%
  rename(.,"dim1"="V1","dim2"="V2") %>%
  dplyr::mutate(.,
                month=sample_metadata$month[match(rownames(.),sample_metadata$sample)],
                month_full=sample_metadata$month_full[match(rownames(.),sample_metadata$sample)],
                lake=sample_metadata$lake[match(rownames(.),sample_metadata$sample)],
                Status=sample_metadata$Status[match(rownames(.),sample_metadata$sample)],
                Distance=sample_metadata$Distance[match(rownames(.),sample_metadata$sample)]) %>%
  dplyr::mutate(.,
                Status=if_else(lake %in% c("CHA","VER"),"Hypereutrophic",
               if_else(lake %in% c("JAB","VAI"),"Mesotrophic",
                       "Eutrophic")),
               Status=factor(Status,levels=c("Hypereutrophic","Eutrophic","Mesotrophic")))

## For lake factor
# Compute the centroids coordinate
Ko_centroids_hellinger.df <- Ko_MDS_hellinger.df %>% dplyr::group_by(lake) %>% dplyr::summarise_at(vars(1:2),median)
# GGplot !
Fig_hellinger_ko <-
  ggplot(Ko_MDS_hellinger.df,aes(dim1,dim2,color = lake,shape=month_full)) +
  geom_convexhull(aes(color = lake,fill=factor(lake),group = lake),alpha=0.1,show.legend = F,size=0.6)+
  geom_point(size=4.5) + theme_bw() +
  geom_label(Ko_centroids_hellinger.df,mapping=aes(dim1,dim2,label = lake, fill = lake),
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
  labs(color="Lake",shape="Month-2021",
       x=paste0("PCoA Axis 1 [",round(MDS_FUNC.hellinger.mds$eig[1]*100/sum(MDS_FUNC.hellinger.mds$eig),1),"%]"),
       y=paste0("PCoA Axis 2 [",round(MDS_FUNC.hellinger.mds$eig[2]*100/sum(MDS_FUNC.hellinger.mds$eig),1),"%]"))
Fig_hellinger_ko

Fig_hellinger_ko<-Fig_hellinger_ko+theme(legend.position = "bottom",legend.direction = "vertical")+
  guides(color = guide_legend(override.aes=list(size=6),nrow=3))
Fig_hellinger_ko

design_F3="
AAAABBBBC
AAAABBBBC
AAAABBBBC
AAAABBBBC
AAAABBBBC
AAAABBBBC
DDDDEEEE#"

final_F3<-Fig_hellinger_ko+Fig_BC_Hellinger+Fig_box_Hellinger+guide_area()+Fig_box_BC+
  plot_layout(design=design_F3,guides='collect')+
  plot_annotation(tag_levels = "A") & theme(plot.tag = element_text(size=15,face="bold"),
                                            legend.direction = "horizontal")
final_F3




Fig_BC_Hellinger+Fig_box_Hellinger+Fig_box_BC+guide_area()


# Merge them together with patchwork

Fig_AD_wout.BC_polygon+Fig_hellinger_ko+
  plot_layout(guides = 'collect')+plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size=15,face="bold"))

Fig_AD_wout.WU_polygon+Fig_jaccard_lake_ko+
  plot_layout(guides = 'collect')+plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size=15,face="bold"))

## For Status factor
# Compute the centroids coordinate
Ko_centroids_hellinger.df <- Ko_MDS_hellinger.df %>% dplyr::group_by(Status) %>% dplyr::summarise_at(vars(1:2),median)
# GGplot !
Fig_hellinger_ko <- ggplot(Ko_MDS_hellinger.df,aes(dim1,dim2,color = Status,shape=month_full)) +
  geom_convexhull(aes(color = Status,fill=factor(Status),group = Status),alpha=0.1,show.legend = F,size=0.6)+
  geom_point(size=4.5) + theme_bw() +
  geom_label(Ko_centroids_hellinger.df,mapping=aes(dim1,dim2,label = Status, fill = Status),color = "white",size = 4.5,
             fontface = "bold",show.legend = F,inherit.aes = F)+
  theme(plot.title = element_text(hjust=0.5,size=15,face = "bold"),
        axis.title = element_text(size=15,face = "bold"),
        axis.text = element_text(size=12),
        legend.title = element_text(size=15,face = "bold"),legend.text = element_text(size=12),
        panel.grid = element_blank(),
        plot.caption = element_text(size=12,face = "italic",color="black"),axis.ticks = element_blank())+
  scale_color_manual(values=palette_lake_type)+scale_fill_manual(values=palette_lake_type)+scale_shape_manual(values=Summer_2021)+
  guides(color = guide_legend(override.aes=list(size=6)),fill = guide_legend(override.aes=list(size=6)))+
  labs(color="Status",shape="Month (2021)",title ="All prokaryotic annotated KO",
       x=paste0("PCoA Axis 1 [",round(MDS_FUNC.hellinger.mds$eig[1]*100/sum(MDS_FUNC.hellinger.mds$eig),1),"%]"),
       y=paste0("PCoA Axis 2 [",round(MDS_FUNC.hellinger.mds$eig[2]*100/sum(MDS_FUNC.hellinger.mds$eig),1),"%]"),
       caption = paste0("hellinger transformed, ",round(MDS_FUNC.hellinger.mds$GOF[[2]],3)*100," % of total expl. var."))
Fig_hellinger_ko

####__ Stats ####
hellinger_factor
####____ Status * Time ####
#merge the factor to the metabolites table
colnames(hellinger_factor)
hellinger_factor<-  data.frame(FUNC.hellinger) %>%
  dplyr::mutate(.,
                lake=Ko_MDS_hellinger.df$lake[match(rownames(.),rownames(Ko_MDS_hellinger.df))],
                Status=Ko_MDS_hellinger.df$Status[match(rownames(.),rownames(Ko_MDS_hellinger.df))],
                Month=Ko_MDS_hellinger.df$month[match(rownames(.),rownames(Ko_MDS_hellinger.df))]) #View(hellinger_factor)
####______ permanova ####
permanova_ko<-adonis2(MDS_FUNC.hellinger~ Status + Month,data=hellinger_factor,permutations=999,method="bray")
summary(permanova_ko)
permanova_ko$`Pr(>F)`[[1]] #General P-value (if ANY sig. diff.)
permanova_ko
####______ Betadisp ####
#Betadisp GENERAL (if NO sig., InterVar >> IntraVar)
betadisper_ko<-as.data.frame(anova(betadisper(d=MDS_FUNC.hellinger,type="centroid",group=hellinger_factor$Status))) #View(betadisper_ASV)
betadisper_ko$`Pr(>F)`[[1]]
####______ Pairwise ####
# sorted Pairwise (adonis between group) 
pairwise_perm_ko_sort<-pairwise.adonis(MDS_FUNC.hellinger,hellinger_factor$Status) %>% .[order(.$pairs),] #View(pairwise_perm_ko_sort)
# sorted Pairwise_Betadisp IF NEDDED (betadisp between group)
pairwise_betadisper_ko_sort<-permutest(betadisper(MDS_FUNC.hellinger,hellinger_factor$Status), pairwise = TRUE) %>% .[[2]] %>% data.frame(.$permuted)
####______ df final ####
final_stats_ko<-tibble(pairwise_perm_ko_sort,betadisp=pairwise_betadisper_ko_sort$permuted)#
View(final_stats_ko)

####__ vs Fun_Status ####
colnames(Hellinger_status.df)
Hellinger_status.df<-Hellinger_matrix_ML %>%
  dplyr::mutate(.,
                Chla1=param_table_AD$Chla_median[match(.$sample1_hellinger,rownames(param_table_AD))],
                Chla2=param_table_AD$Chla_median[match(.$sample2_hellinger,rownames(param_table_AD))],
                Gradient1=param_MDS_eucli.df$dim1[match(.$sample1_hellinger,rownames(param_MDS_eucli.df))],
                Gradient2=param_MDS_eucli.df$dim1[match(.$sample2_hellinger,rownames(param_MDS_eucli.df))],
                lake_1=sample_metadata$lake[match(sample1_hellinger,sample_metadata$lake_month)],
                lake_2=sample_metadata$lake[match(sample2_hellinger,sample_metadata$lake_month)],
                month_1=A_D_rarW@sam_data$month_number[match(sample1_hellinger,A_D_rarW@sam_data$lake_month)],
                month_2=A_D_rarW@sam_data$month_number[match(sample2_hellinger,A_D_rarW@sam_data$lake_month)],
                month_gap=abs(month_1-month_2),
                Status=sample_metadata$Status[match(sample1_hellinger,sample_metadata$lake_month)]) %>%
  dplyr::mutate(.,
                Status_dist = rowMeans(.[,c(9,10)]),
                Chla_dist= rowMeans(.[,c(7,8)]),
                Havers_dist=Fig_BC_Geo$data$Havers_dist[match(.$lake_1,Fig_BC_Geo$data$G1_lake)]) %>%
  subset(.,sample1_hellinger != sample2_hellinger) %>% #subset(.,lake_1 == lake_2) %>%
  dplyr::filter(duplicated(hellinger_dist) == FALSE)


spearman_KO_status<-cor(Hellinger_status.df$hellinger_dist,Hellinger_status.df$Chla_dist,method = "spearman")

spearman_KO_status<-cor.test(Hellinger_status.df$hellinger_dist,Hellinger_status.df$Chla_dist,method = "spearman")

lm_KO_status_model <- lm(hellinger_dist~Chla_dist,data=Hellinger_status.df)

var_KO_status=paste0("p <0.01, R2 ",round(summary(lm_KO_status_model)[["adj.r.squared"]],2))

Fig_Hellinger_Status<-Hellinger_status.df %>% 
  ggplot(.,aes(Chla_dist,hellinger_dist,group=1))+
  geom_point(size=2,shape=21,fill="grey",color="black")+
  stat_smooth(method = "lm",color="#4166F5",fill="lightblue",show.legend = F)+
  theme_classic()+
  theme(axis.ticks = element_blank(),
        axis.text = element_text(size=10),
        axis.title=element_text(size=12,face="bold"),aspect.ratio = 1,
        plot.caption = element_text(size=12,face="italic"))+
  scale_x_continuous(trans="log10",breaks = c(0,5,10,40,100,200))+
  #scale_color_manual(values=palette_slope)+
  #scale_fill_manual(values=palette_lake_type)+
  labs(x=expression(paste(bold("Chlorophyll "),bold(italic("a")),bold(" ["),bold(µg.L^"-1"))),
       caption=var_KO_status,
       y="Gene-content dissimilarity")
Fig_Hellinger_Status

spearman_KO_geo<-cor(Hellinger_status.df$hellinger_dist,Hellinger_status.df$Havers_dist,method = "spearman")

spearman_KO_geo<-cor.test(Hellinger_status.df$hellinger_dist,Hellinger_status.df$Havers_dist,method = "spearman")

lm_KO_status_geo <- lm(hellinger_dist~Havers_dist,data=Hellinger_status.df)
summary(lm_KO_status_geo)
var_KO_geo=paste0("p>0.05")

Fig_Hellinger_Geo<-Hellinger_status.df %>%
  dplyr::mutate(.,Havers_dist=(Havers_dist/1000)) %>%
  ggplot(.,aes(Havers_dist,hellinger_dist,group=1))+
  geom_point(size=2,shape=21,fill="grey",color="black")+
  stat_smooth(method = "lm",color="#4166F5",fill="lightblue",show.legend = F)+
  theme_classic()+
  theme(axis.ticks = element_blank(),
        axis.text = element_text(size=12),
        axis.title=element_text(size=10,face="bold"),aspect.ratio = 1,
        plot.caption = element_text(size=12,face="italic"))+
  #scale_color_manual(values=palette_slope)+
  #scale_fill_manual(values=palette_lake_type)+
  annotate(geom="text", x=1, y=4,color="black",size=5,hjust=0,label="paste(italic(p), \"<0.01 \", italic(R) ^ 2, \" 0.74\")", parse = TRUE)+
  labs(x="Between-lake distance (km)",caption=var_KO_geo,
       y="Gene-content dissimilarity")
Fig_Hellinger_Geo

spearman_KO_month<-cor(Hellinger_status.df$hellinger_dist,Hellinger_status.df$month_gap,method = "spearman")

lm_KO_month <- lm(hellinger_dist~month_gap,data=Hellinger_status.df)

var_KO_month=paste0("Spearman ",round(spearman_KO_month,2),", R2 ",round(summary(lm_KO_month)[["adj.r.squared"]],2),", p >0.05")

Fig_Hellinger_month<-Hellinger_status.df %>%
  filter(month_gap !=0) %>%
  ggplot(.,aes(month_gap,hellinger_dist,group=1))+
  geom_point(size=2,shape=21,fill="grey",color="black")+
  stat_smooth(method = "lm",color="#4166F5",fill="lightblue",show.legend = F)+
  theme_classic()+
  theme(axis.ticks = element_blank(),
        axis.text = element_text(size=12),
        axis.title=element_text(size=12,face="bold"),aspect.ratio = 1,
        plot.caption = element_text(size=12,face="italic"))+
  scale_color_manual(values=palette_lake)+
  scale_fill_manual(values=palette_lake)+
  labs(x="Number of month between samples",caption=var_KO_month,
       y= "Functional dissimilarity")
Fig_Hellinger_month

####. ####

#### BGC Barplot ####

####__ KO MDS ####

# Transfrom the .df from long to wide
MDS_BGC.df <- all_fun.df %>% subset(KEGG_ko %in% BGC_short_list) %>%
  subset(month %in% c("A","B","C","D")) %>% .[,c(1,2,9)] %>%
  group_by(KEGG_ko,sample) %>% dplyr::summarise(CPM=sum(CPM)) %>%
  pivot_wider(names_from = sample, values_from = CPM,values_fill = 0) %>% remove_rownames() %>%
  column_to_rownames(.,var = "KEGG_ko") %>% t()

####__ Bray-Curtis ####

# Hellinger abudnace matrix with vegan
list_K0<-colnames(MDS_BGC.df)
MDS_BGC.hellinger<- MDS_BGC.df %>% labdsv::hellinger(.) %>% as.data.frame(.) %>% magrittr::set_colnames(list_K0) %>%
  vegdist(., method="bray")
length(list_K0)
# Compute MDS
MDS_BGC.hellinger.mds<- cmdscale(MDS_BGC.hellinger,eig=TRUE, k=2)

# Store the axis coordinates with metadat for ggplot
BGC_MDS_hellinger.df <- MDS_BGC.hellinger.mds$points %>% as.data.frame(.) %>%
  rename(.,"dim1"="V1","dim2"="V2") %>%
  dplyr::mutate(.,
                month=sample_metadata$month[match(rownames(.),sample_metadata$sample)],
                month_full=sample_metadata$month_full[match(rownames(.),sample_metadata$sample)],
                lake=sample_metadata$lake[match(rownames(.),sample_metadata$sample)],
                Status=sample_metadata$Status[match(rownames(.),sample_metadata$sample)])

####____ Status * Time ####
#merge the factor to the metabolites table

BGC_factor<-  data.frame(BGC.hellinger) %>%
  dplyr::mutate(.,
                lake=BGC_MDS_hellinger.df$lake[match(rownames(.),rownames(BGC_MDS_hellinger.df))],
                Status=BGC_MDS_hellinger.df$Status[match(rownames(.),rownames(BGC_MDS_hellinger.df))],
                Month=BGC_MDS_hellinger.df$month[match(rownames(.),rownames(BGC_MDS_hellinger.df))]) #View(hellinger_factor)
####______ permanova ####
permanova_ko<-adonis2(BGC.hellinger~ Status + Month,data=BGC_factor,permutations=999,method="bray")
permanova_ko$`Pr(>F)`[[1]] #General P-value (if ANY sig. diff.)
permanova_ko

BGC.hellinger<- MDS_BGC.df %>% labdsv::hellinger(.) %>% as.data.frame(.) %>% magrittr::set_colnames(list_K0)

####|####

####__ BGC short ####
BGC_short<-read_delim("BGC_PF.csv",delim = ";",show_col_types = F)
View(BGC_short)
BGC_short_list<-BGC_short$Kegg_ko

BGC_short$Level_2

####__ BGC boxplot ####

BGC_short.df<- all_fun.df %>%
  subset(month %in% c("A","B","C","D")) %>%
  dplyr::group_by(sample) %>% dplyr::mutate(.,CPM_tot = sum(CPM)) %>%
  subset(KEGG_ko %in% BGC_short_list) %>%
  cbind(.,
        gene=BGC_short$Level_2[match(.$KEGG_ko,BGC_short$Kegg_ko)],
        gene_metabolism=BGC_short$Level_1[match(.$KEGG_ko,BGC_short$Kegg_ko)]) %>%
  dplyr::group_by(sample,gene_metabolism,gene) %>% dplyr::mutate(.,CPM_rel = (CPM/CPM_tot)*100) %>%
  dplyr::group_by(sample,gene_metabolism,gene) %>% dplyr::summarise(CPM_sum=sum(CPM_rel)) %>%
  cbind(.,
        month_full=sample_metadata$month_full[match(.$sample,sample_metadata$sample)],
        Status=sample_metadata$Status[match(.$sample,sample_metadata$sample)],
        lake=sample_metadata$lake[match(.$sample,sample_metadata$sample)]) %>%
  dplyr::mutate(.,
                gene_metabolism=factor(gene_metabolism,
                               levels = c("Carbon_fixation","N_metabolism",
                               "P_metabolism","S_metabolism","Iron-related")),
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
                                          "Vaires","Jablines")))


BGC_BOI.boxplot<- BGC_short.df %>% subset(.,gene_metabolism %in% c("Carbon_fixation")) %>%
  ggplot(.,aes(lake,CPM_sum,fill=gene))+
  geom_boxplot(width=0.15,color="black",show.legend = F,outlier.size = 0.5)+theme_bw()+
  geom_point(show.legend = F,color="black",size=0.5)+
  theme(panel.grid = element_blank(),axis.ticks = element_blank(),
        axis.title.x = element_blank(),axis.title.y = element_blank(),
        axis.text.y = element_text(size=10),
        axis.text.x = element_blank())+
  #scale_fill_manual(values=c("#7BA8CF","#3E842E"))+
  scale_y_continuous(limits = c(0,NA),name = NULL, sec.axis = sec_axis(~.))+
  guides(y = "none")+
  facet_grid2(gene ~ lake,scales = "free")
              # strip =
              #   strip_color_BGC_CP<- strip_themed(background_x = elem_list_rect(fill = c("#2F6B9D","#E45000","#28A448"),color = "black"),
              #                                     text_x = elem_list_text(colour = "white",face = "bold",size=10),
              #                                     background_y = elem_list_rect(fill = c("#7BA8CF","#3E842E"),color = "black"),
              #                                     text_y = elem_list_text(colour = "white",face = "bold",size=9)),switch = "y")

BGC_BOI.boxplot


BGC_short.boxplot<- BGC_short.df %>% subset(.,gene_metabolism %in% c("Carbon_fixation","P_metabolism")) %>%
  dplyr::mutate(.,gene_metabolism=if_else(gene_metabolism=="Carbon_fixation","C fixation","P metabolism")) %>%
  ggplot(.,aes(Status,CPM_sum,fill=gene_metabolism))+
  geom_boxplot(width=0.15,color="black",show.legend = F,outlier.size = 0.5)+theme_bw()+
  geom_point(show.legend = F,color="black",size=0.5)+
  theme(panel.grid = element_blank(),axis.ticks = element_blank(),
        axis.title.x = element_blank(),axis.title.y = element_blank(),
        axis.text.y = element_text(size=10),
        axis.text.x = element_blank())+
  scale_fill_manual(values=c("#7BA8CF","#3E842E"))+
  scale_y_continuous(limits = c(0,NA),name = NULL, sec.axis = sec_axis(~.))+
  guides(y = "none")+
  facet_grid2(gene_metabolism ~ Status,scales = "free",
              strip =
                strip_color_BGC_CP<- strip_themed(background_x = elem_list_rect(fill = c("#2F6B9D","#E45000","#28A448"),color = "black"),
                                                  text_x = elem_list_text(colour = "white",face = "bold",size=10),
                                                  background_y = elem_list_rect(fill = c("#7BA8CF","#3E842E"),color = "black"),
                                                  text_y = elem_list_text(colour = "white",face = "bold",size=9)),switch = "y")

BGC_short.boxplot

View(all_fun.df %>% subset(KEGG_ko=="ko:K04643"))


#### here ####
all_fun.df %>%
  subset(month %in% c("A","B","C","D")) %>%
  dplyr::group_by(lake, month) %>% 
  dplyr::mutate(.,CPM_tot = sum(CPM)) %>%
  subset(KEGG_ko %in% BGC_short_list) %>%
  dplyr::group_by(lake,month) %>%
  dplyr::mutate(.,CPM_sum = (sum(CPM)/CPM_tot)*100,
                summer="yes") %>%
  group_by(summer) %>%
  dplyr::summarise_at(vars("CPM_sum"),list(median=median, sd=sd))
  


BGC_short_tax.df<- all_fun.df %>% subset(KEGG_ko %in% BGC_short_list) %>%
  subset(month %in% c("A","B","C","D")) %>%
  dplyr::mutate(.,
                gene=BGC_short$Level_2[match(.$KEGG_ko,BGC_short$Kegg_ko)],
                gene_metabolism=BGC_short$Level_1[match(.$KEGG_ko,BGC_short$Kegg_ko)]) %>%
  dplyr::group_by(sample,gene_metabolism,gene,phylum) %>% dplyr::summarise(.,CPM_sum = sum(CPM)) %>%
  cbind(.,
        month_full=sample_metadata$month_full[match(.$sample,sample_metadata$sample)],
        Status=sample_metadata$Status[match(.$sample,sample_metadata$sample)],
        lake=sample_metadata$lake[match(.$sample,sample_metadata$sample)]) %>%
  dplyr::mutate(.,
                gene_metabolism=factor(gene_metabolism,
                                       levels = c("Carbon_fixation","N_metabolism",
                                                  "P_metabolism","S_metabolism","Iron-related")),
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
                                          "Vaires","Jablines")))

#### SIMPER ####

####__ BLR ####

simper_BGC_com<- all_fun.df %>% subset(KEGG_ko %in% BGC_short_list) %>%
  subset(month %in% c("A","B","C","D")) %>%
  dplyr::mutate(.,
                gene=BGC_short$Level_2[match(.$KEGG_ko,BGC_short$Kegg_ko)],
                gene_metabolism=BGC_short$Level_1[match(.$KEGG_ko,BGC_short$Kegg_ko)]) %>%
  dplyr::group_by(sample,gene) %>% dplyr::summarise(.,CPM_sum = sum(CPM)) %>%
  pivot_wider(values_from =CPM_sum,values_fill = 0,names_from = gene) %>%
  cbind(.,Status=BGC_short.df$Status[match(.$sample,BGC_short.df$sample)]) %>%
  column_to_rownames(var = "sample")

simper_BGC_BOI<- all_fun.df %>% subset(KEGG_ko %in% BGC_short_list) %>%
  subset(month %in% c("A","B","C","D")) %>%
  dplyr::mutate(.,
                gene=BGC_short$Level_2[match(.$KEGG_ko,BGC_short$Kegg_ko)],
                gene_metabolism=BGC_short$Level_1[match(.$KEGG_ko,BGC_short$Kegg_ko)]) %>%
  dplyr::group_by(sample,gene) %>% dplyr::summarise(.,CPM_sum = sum(CPM)) %>%
  pivot_wider(values_from =CPM_sum,values_fill = 0,names_from = gene) %>%
  cbind(.,BOI=BGC_short.df$lake[match(.$sample,BGC_short.df$sample)]) %>%
  dplyr::mutate(.,BOI=ifelse(BOI=="BOI","BOI","others")) %>%
  column_to_rownames(var = "sample")
ncol(simper_BGC_BOI)

simper_BGC_BOI.df<-simper(simper_BGC_BOI[,c(-29)],simper_BGC_BOI$BOI,permutations = 999,parallel = 5)
summary(simper_BGC_BOI.df)

simper_per_BOI<- summary(simper_BGC_BOI.df)[[1]] %>% as.data.frame()
percent_BGC_BOI<-simper_per_BOI$cumsum[1]*100

for (i in 2:nrow(simper_per_BOI)) {
  percent<-(simper_per_BOI$cumsum[i]*100-simper_per_BOI$cumsum[i-1]*100)
  percent_BGC_BOI<-c(percent_BGC_BOI,percent)
}
simper_per_BOI$percent<-percent_BGC_BOI

Fig_simper_BOI<-simper_per_BOI %>% rownames_to_column(var="gene") %>%
  cbind(.,gene_metabolism=BGC_short$Level_1[match(.$gene,BGC_short$Level_2)]) %>%
  group_by(gene_metabolism) %>% dplyr::summarise(sum=round(sum(percent),1)) %>%
  dplyr::mutate(.,
                gene_metabolism=
                  if_else(gene_metabolism=="Carbon_fixation","C metabolism",
                          if_else(gene_metabolism=="P_metabolism","P metabolism",
                                  if_else(gene_metabolism=="N_metabolism","N metabolism",
                                          if_else(gene_metabolism=="S_metabolism","S metabolism",
                                                  if_else(gene_metabolism=="Methane_metabolism","Methane metabolism",
                                                          if_else(gene_metabolism=="Methane_metabolism","C metabolism","Iron metabolism")))))),
                gene_metabolism=factor(gene_metabolism,
                                       levels = c("C metabolism","N metabolism","P metabolism",
                                                  "S metabolism","Iron metabolism"))) %>%
  arrange(desc(gene_metabolism)) %>%
  dplyr::mutate(.,
                prop = sum / sum(sum)*100,
                ypos = cumsum(prop)- 0.5*prop) %>%
  ggplot(., aes(x=0, y=sum, fill=gene_metabolism)) +
  geom_bar(stat="identity", width=1, color="black") +
  geom_text(aes(y = ypos,
                label = sum), color = "black", size=6) +
  theme_void()+
  theme(legend.title = element_text(size=15,face='bold'),
        legend.text = element_text(size=12),
        legend.position = "bottom",legend.direction = "vertical")+
  scale_fill_manual(values = palette_BGC_short) +
  guides(fill=guide_legend(override.aes=list(size=5,shape=15),nrow=3))+
  coord_polar("y", start=0) +
  labs(fill="BGC metabolism") # remove background, grid, numeric labels
Fig_simper_BOI

####__ Trophic Status ####
ncol(simper_BGC_com)
simper_BGC.df<-simper(simper_BGC_com[,c(-29)],simper_BGC_com$Status,permutations = 999,parallel = 5)
summary(simper_BGC.df)

# Contraste H _ M
simper_H_M <- summary(simper_BGC.df, ordered =T)[[3]] %>% as.data.frame() %>%
  subset(.,p <=0.05)
list_simper_H_M<-rownames(simper_H_M)

# per metabolism
#simper_per_metabolism<- simper_H_M

simper_per_metabolism<- summary(simper_BGC.df)[[3]]
#simper_per_metabolism<- summary(simper_BGC.df)[[3]] %>% as.data.frame()
percent_BGC_H_M<-simper_per_metabolism$cumsum[1]*100

for (i in 2:nrow(simper_per_metabolism)) {
  percent<-(simper_per_metabolism$cumsum[i]*100-simper_per_metabolism$cumsum[i-1]*100)
  percent_BGC_H_M<-c(percent_BGC_H_M,percent)
}
simper_per_metabolism$percent<-percent_BGC_H_M

Fig_simper_H_M<-simper_per_metabolism %>% rownames_to_column(var="gene") %>%
  cbind(.,gene_metabolism=BGC_short$Level_1[match(.$gene,BGC_short$Level_2)]) %>%
  group_by(gene_metabolism) %>% dplyr::summarise(sum=round(sum(percent),1)) %>%
  dplyr::mutate(.,
                gene_metabolism=
                if_else(gene_metabolism=="Carbon_fixation","C metabolism",
                        if_else(gene_metabolism=="P_metabolism","P metabolism",
                                if_else(gene_metabolism=="N_metabolism","N metabolism",
                                        if_else(gene_metabolism=="S_metabolism","S metabolism",
                                          if_else(gene_metabolism=="Methane_metabolism","Methane metabolism",
                                                  if_else(gene_metabolism=="Methane_metabolism","C metabolism","Iron metabolism")))))),
                gene_metabolism=factor(gene_metabolism,
                                       levels = c("C metabolism","N metabolism","P metabolism",
                                                  "S metabolism","Iron metabolism"))) %>%
  arrange(desc(gene_metabolism)) %>%
  dplyr::mutate(.,
                prop = sum / sum(sum)*100,
                ypos = cumsum(prop)- 0.5*prop) %>%
  ggplot(., aes(x=0, y=sum, fill=gene_metabolism)) +
  geom_bar(stat="identity", width=1, color="white") +
  geom_text(aes(y = ypos,
                label = sum), color = "white", size=6) +
  theme_void()+
  theme(legend.title = element_text(size=15,face="bold"),
        legend.text = element_text(size=12),
        legend.position = "bottom",legend.direction = "vertical")+
  scale_fill_manual(values = palette_BGC_short) +
  guides(fill=guide_legend(override.aes=list(size=5,shape=15),nrow=3))+
  coord_polar("y", start=0) +
  labs(fill="BGC metabolism")  # remove background, grid, numeric labels
Fig_simper_H_M

#per significant genes 

simper_per_gene<- summary(simper_BGC.df, ordered =T)[[3]]
percent_gene_H_M<-simper_per_gene$cumsum[1]*100

for (i in 2:nrow(simper_per_gene)) {
  percent<-(simper_per_gene$cumsum[i]*100-simper_per_gene$cumsum[i-1]*100)
  percent_gene_H_M<-c(percent_gene_H_M,percent)
}
simper_per_gene$percent<-percent_gene_H_M 

simper_per_gene<- simper_per_gene %>% subset(.,p <=0.05)

Fig_simper_H_M_gene<-simper_per_gene %>% rownames_to_column(var="gene") %>%
  cbind(.,gene_metabolism=BGC_short$Level_1[match(.$gene,BGC_short$Level_2)])  %>%
  #group_by(Phylum) %>% dplyr::summarise(sum_var=round(sum(percent),1)) %>%
  add_row(.,gene_metabolism = 'n.s.',gene="n.s.",percent = 100-sum(.$percent)) %>%
  dplyr::mutate(.,
                ypos = cumsum(percent)- 0.5*percent,
                xpos= 0.5,
                gene_percent=paste0(gene,"_",round(percent,1))) %>%
  ggplot(., aes(x=0.5, y=percent, fill=fct_reorder(gene,percent)))+
  geom_bar(stat="identity", width=1, color="black",show.legend = F) +
  geom_text(aes(y = ypos,x = xpos,
                label = gene_percent),fontface = "bold",color = "black", size=4) +
  theme_void()+
  theme(legend.title = element_text(size=15,face="bold"),
        legend.text = element_text(size=12),
        legend.position = "bottom",legend.direction = "vertical")+
  scale_fill_manual(values = rev(c("#5AAA46","gray38","#5AAA46","#C84D4C","#5AAA46","#9B5C97"))) +
  guides(fill=guide_legend(override.aes=list(size=5,shape=15),nrow=3))+
  coord_polar("y", start=0) +
  labs(fill="Microbial Phyla")  # remove background, grid, numeric labels
Fig_simper_H_M_gene

palette_BGC_short<-c("#317EC2","#9B5C97","#5AAA46","#C84D4C","grey","black")

# Contraste E _ M
simper_E_M <- summary(simper_BGC.df, ordered =T)[[2]] %>% as.data.frame() %>%
  subset(., p <=0.05)
list_simper_E_M<-rownames(simper_E_M)

#simper_per_metabolism<- simper_E_M

simper_per_metabolism<- summary(simper_BGC.df)[[2]]

#simper_per_metabolism<- summary(simper_BGC.df)[[3]] %>% as.data.frame()
percent_BGC_E_M<-simper_per_metabolism$cumsum[1]*100

for (i in 2:nrow(simper_per_metabolism)) {
  percent<-(simper_per_metabolism$cumsum[i]*100-simper_per_metabolism$cumsum[i-1]*100)
  percent_BGC_E_M<-c(percent_BGC_E_M,percent)
}
simper_per_metabolism$percent<-percent_BGC_E_M


# Contraste E _ H
simper_E_H <- summary(simper_BGC.df, ordered =T)[[1]] %>% as.data.frame() %>%
  subset(.,p <=0.05)
list_simper_E_H<-rownames(simper_E_H)

#simper_per_metabolism<- simper_E_H

simper_per_metabolism<- summary(simper_BGC.df)[[1]]

#simper_per_metabolism<- summary(simper_BGC.df)[[3]] %>% as.data.frame()
percent_BGC_E_H<-simper_per_metabolism$cumsum[1]*100

for (i in 2:nrow(simper_per_metabolism)) {
  percent<-(simper_per_metabolism$cumsum[i]*100-simper_per_metabolism$cumsum[i-1]*100)
  percent_BGC_E_H<-c(percent_BGC_E_H,percent)
}
simper_per_metabolism$percent<-percent_BGC_E_H

# simper BGCs list p 0.05
list_simper<-unique(c(list_simper_H_M, list_simper_E_M,list_simper_E_H)) %>% as.data.frame() %>%
  rename("BGC_simper"=".")
length(list_simper)

####__ Carbon & Methane BGC tax ####

C_tax<-
  all_fun.df %>% subset(month %in% c("A","B","C","D")) %>%
  cbind(.,gene=BGC_short$Level_2[match(.$KEGG_ko,BGC_short$Kegg_ko)]) %>%
  dplyr::group_by(lake,month) %>% dplyr::mutate(.,CPM_tot = sum(CPM)) %>%  # obtain total CPM
  dplyr::group_by(lake,month,gene,phylum) %>% dplyr::summarise(.,CPM_sum = sum((CPM/CPM_tot)*100)) %>% # merge phyla by gene and sample
  subset(gene %in% c("rbcL","psbA","pufM","fmoA","coxL","aclA","tauA","mcrA","mmoB"))  %>%
  dplyr::mutate(.,gene=factor(gene,levels =c("rbcL","psbA","pufM","fmoA","coxL","aclA","tauA","mcrA","mmoB")),
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
                                              "VSM","JAB"))),
                month_full=if_else(month=="A","Jun.",
                                   if_else(month=="B","Jul.",
                                           if_else(month=="C","Aug.","Sep."))),
                month_full=factor(month_full,
                                  levels=c("Jun.","Jul.","Aug.","Sep."))) %>%
  # dplyr::group_by(sample) %>% dplyr::mutate(.,CPM_tot = sum(CPM)) %>%
  # dplyr::group_by(sample,gene,phylum) %>% dplyr::mutate(.,CPM_sum = sum(CPM)) %>%
  # dplyr::group_by(lake,gene,phylum) %>% dplyr::summarise(CPM_median=median(CPM_sum)) %>%
  # dplyr::group_by(lake,gene) %>% dplyr::mutate(CPM_rel=((CPM_median/sum(CPM_median))*100)) %>%
  dplyr::mutate(.,
                phylum = if_else(phylum %in% c("Cyanobacteria","Chlorobi",
                                               "Actinobacteria","Proteobacteria","Candidatus Gracilibacteria","Bacteria",
                                               "Euryarchaeota"),phylum,"Others")) %>%
  dplyr::group_by(lake,month_full,gene,phylum) %>%  dplyr::summarise(CPM_taxa_rel = sum(CPM_sum)) %>%
  ggplot(.,aes(month_full,CPM_taxa_rel,fill=fct_rev(fct_reorder(phylum,CPM_taxa_rel))))+
  geom_col(width=0.8,color="black",size=0.1)+theme_bw()+
  theme(panel.grid = element_blank(),
        legend.title = element_text(face="bold",size=15),
        legend.text = element_text(size=12),
        axis.title = element_blank(),axis.ticks = element_blank(),
        axis.text.x = element_text(size=9),axis.text.y = element_text(size=9))+
  scale_fill_manual(values=c("#CCEBC5","#FED9A6","#B3CDE3","firebrick","gray38",
                             "lightgrey","black","salmon"))+
                             #  FED9A6
                             # "#F2F2F2","#000000","#E5D8BD","firebrick",
                             # 
                             # "#70F4BB","#5D9BEF","#F391CC","#EEF4A3","#CAAB7D","#8AF3D4","#EAF7B9","#4263B2","#E0B333","#2CB550",
                             # "#BE6AE7","#8B7664","#7161EC","#CACCA8","#F138F4","#8EAC73"))+
  guides(fill=guide_legend(ncol = 1,override.aes = list(shape=22,size=5)))+
  #scale_x_discrete()+
  scale_y_continuous(position = 'right')+
  labs(fill="Prokaryotic Phyla",y="Relative abundance (CPM)")+
  facet_grid2(gene~lake,scales = "free",switch = "y",
              strip=
                strip_themed(background_x = elem_list_rect(fill = rev(palette_lake_3T),color = "black"),
                             text_x = elem_list_text(colour = "white",face = "bold",size=10),
                             background_y = elem_list_rect(fill = c("#317EC2","#317EC2","#317EC2","#317EC2",
                                                                    "#317EC2","#317EC2","#317EC2",
                                                                    "grey","grey"),color = "black"),
                             text_y = elem_list_text(colour = "white",face = "bold.italic",size=14)))
C_tax

C_tax$data %>% subset(gene =="mcrA") %>%
  dplyr::group_by(lake,month_full) %>%
  dplyr::mutate(.,CPM_tot=sum(CPM_taxa_rel)) %>%
  dplyr::mutate(.,Status=if_else(lake %in% c("VSS","CSM"),"Hypereutrophic",
                                 if_else(lake %in% c("JAB","VSM"),"Mesotrophic",
                                         "Eutrophic"))) %>%
  dplyr::group_by(gene) %>%
  dplyr::summarise_at(vars("CPM_tot"),list(median=median, sd=sd))

C_tax$data %>% subset(gene %in% c("rbcL")) %>%
  dplyr::group_by(lake,month_full) %>%
  dplyr::mutate(.,CPM_tot=sum(CPM_taxa_rel)) %>%
  dplyr::mutate(.,Status=if_else(lake =="BLR","BLR","other"),
                summer="yes") %>%
  group_by(Status) %>%
  dplyr::summarise_at(vars("CPM_tot"),list(mean=mean, sd=sd))

C_tax$data %>% subset(gene %in% c("mcrA")) %>%
  cbind(.,gene_metabolism=BGC_short$Level_1[match(.$gene,BGC_short$Level_2)]) %>%
  dplyr::group_by(lake,month_full,gene_metabolism) %>%
  dplyr::mutate(.,CPM_tot=sum(CPM_taxa_rel)) %>%
  dplyr::mutate(.,Status=if_else(lake %in% c("VSS","CSM"),"Hypereutrophic",
                                 if_else(lake %in% c("JAB","VSM"),"Mesotrophic",
                                         "Eutrophic"))) %>%
  dplyr::group_by(gene_metabolism) %>%
  dplyr::summarise_at(vars("CPM_tot"),list(median=median, sd=sd))

C_tax$data %>% subset(gene =="rbcL") %>%
  dplyr::group_by(lake,month_full) %>%
  dplyr::mutate(.,CPM_tot=sum(CPM_taxa_rel)) %>%
  dplyr::mutate(.,Status=if_else(lake %in% c("BLR"),"BLR","others")) %>%
  dplyr::group_by(Status) %>%
  dplyr::summarise_at(vars("CPM_tot"),list(mean=mean, sd=sd))

C_tax$data %>% subset(gene =="rbcL") %>%
  dplyr::group_by(lake,month_full) %>%
  dplyr::mutate(.,CPM_tot=sum(CPM_taxa_rel)) %>%
  dplyr::group_by(lake,month_full,phylum) %>%
  dplyr::mutate(.,percent_CPM=CPM_taxa_rel/CPM_tot*100) %>%
  dplyr::group_by(phylum) %>%
  dplyr::summarise_at(vars("percent_CPM"),list(mean=mean, sd=sd))

C_tax$data %>% subset(gene =="rbcL") %>%
  dplyr::mutate(.,CPM_tot=sum(CPM_taxa_rel)) %>%
  dplyr::mutate(.,Status=if_else(lake %in% c("BLR"),"BLR","others")) %>%
  dplyr::group_by(Status,phylum) %>%
  #dplyr::group_by(lake,phylum) %>%
  dplyr::mutate(.,percent_CPM=CPM_taxa_rel/CPM_tot) %>%
  dplyr::group_by(phylum) %>%
  dplyr::summarise_at(vars("percent_CPM"),list(mean=mean, sd=sd))

####__ Nitrogen BGC tax ####

N_tax<-
  all_fun.df %>% subset(month %in% c("A","B","C","D")) %>%
  cbind(.,gene=BGC_short$Level_2[match(.$KEGG_ko,BGC_short$Kegg_ko)]) %>%
  dplyr::group_by(lake,month) %>% dplyr::mutate(.,CPM_tot = sum(CPM)) %>%  # obtain total CPM
  dplyr::group_by(lake,month,gene,phylum) %>% dplyr::summarise(.,CPM_sum = sum((CPM/CPM_tot)*100)) %>% # merge phyla by gene and sample
  subset(gene %in% c("nifH","anfG","narB","nasA","ureC","amoA"))  %>%
  dplyr::mutate(.,gene=factor(gene,levels=c("nifH","anfG","narB","nasA","ureC","amoA")),
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
                                              "VSM","JAB"))),
                month_full=if_else(month=="A","Jun.",
                                   if_else(month=="B","Jul.",
                                           if_else(month=="C","Aug.","Sep."))),
                month_full=factor(month_full,
                                  levels=c("Jun.","Jul.","Aug.","Sep."))) %>%
  dplyr::mutate(.,
               phylum = if_else(phylum %in% c("Cyanobacteria","Chlorobi","Planctomycetes",
                                              "Actinobacteria","Proteobacteria"),phylum,"Others")) %>%
  dplyr::group_by(lake,month_full,gene,phylum) %>%  dplyr::summarise(CPM_taxa_rel = sum(CPM_sum)) %>%
  ggplot(.,aes(month_full,CPM_taxa_rel,fill=fct_rev(fct_reorder(phylum,CPM_taxa_rel))))+
  geom_col(width=0.8,color="black",size=0.1)+theme_bw()+
  theme(panel.grid = element_blank(),
        legend.title = element_text(face="bold",size=14),legend.text = element_text(size=12),
        axis.title = element_blank(),axis.ticks = element_blank(),
        axis.text.x = element_text(size=9),axis.text.y = element_text(size=9))+
  scale_fill_manual(values=c("#FED9A6","#CCEBC5","#B3CDE3","firebrick","lightgrey","#E5D8BD"))+
  guides(fill=guide_legend(ncol = 1,override.aes = list(shape=22,size=5)))+
  scale_y_continuous(position = 'right')+
  labs(fill="Prokaryotic Phyla",y="Relative abundance (CPM)")+
  facet_grid2(gene~lake,scales = "free",switch = "y",
              strip=
                strip_themed(background_x = elem_list_rect(fill = rev(palette_lake_3T),color = "black"),
                             text_x = elem_list_text(colour = "white",face = "bold",size=10),
                             background_y = elem_list_rect(fill = c("#9B5C97","#9B5C97","#9B5C97","#9B5C97",
                                                                    "#9B5C97","#9B5C97"),color = "black"),
                             text_y = elem_list_text(colour = "white",face = "bold.italic",size=14)))
N_tax

tmp<-
  all_fun.df %>% subset(month %in% c("A","B","C","D")) %>%
  cbind(.,gene=BGC_short$Level_2[match(.$KEGG_ko,BGC_short$Kegg_ko)]) %>%
  dplyr::group_by(lake,month) %>% dplyr::mutate(.,CPM_tot = sum(CPM)) %>%  # obtain total CPM
  dplyr::group_by(lake,month,gene) %>% dplyr::summarise(.,CPM_sum = sum((CPM/CPM_tot)*100)) %>% # merge phyla by gene and sample
  subset(gene %in% c("asrB","aprA","dsrB","soxB","dmdA","fecA"))  %>%
  dplyr::mutate(.,lake_month=paste0(lake,"_",month)) %>%
  .[,c(3,4,5)] %>% pivot_wider(names_from = "lake_month",values_from = "CPM_sum",values_fill = NA)

View(tmp)
write.csv(tmp,"tmp.csv")
  
N_tax$data %>% subset(gene =="narB") %>%
  dplyr::group_by(lake,month_full) %>%
  dplyr::mutate(.,CPM_tot=sum(CPM_taxa_rel)) %>%
  dplyr::mutate(.,Status=if_else(lake %in% c("VSS","CSM"),"Hypereutrophic",
                                 if_else(lake %in% c("JAB","VSM"),"Mesotrophic",
                                         "Eutrophic")),
                summer="yes") %>%
  dplyr::group_by(Status) %>%
  dplyr::summarise_at(vars("CPM_tot"),list(mean=mean, sd=sd))

P_tax$data %>% subset(gene =="ppk1") %>%
  dplyr::group_by(lake,month_full) %>%
  dplyr::mutate(.,CPM_tot=sum(CPM_taxa_rel)) %>%
  dplyr::mutate(.,Status=if_else(lake %in% c("JAB","VSM"),"Hypereutrophic","others")) %>%
  dplyr::group_by(Status) %>%
  dplyr::summarise_at(vars("CPM_tot"),list(median=median, sd=sd))

N_tax$data %>% subset(gene =="ureC") %>%
  dplyr::group_by(lake,month_full) %>%
  dplyr::mutate(.,CPM_tot=sum(CPM_taxa_rel)) %>%
  dplyr::mutate(.,Status=if_else(lake %in% c("BLR"),"BLR","others")) %>%
  dplyr::group_by(Status) %>%
  dplyr::summarise_at(vars("CPM_tot"),list(mean=mean, sd=sd))

N_tax$data %>% subset(gene =="ureC" & lake !="BLR") %>%
  dplyr::group_by(lake,month_full) %>%
  dplyr::mutate(.,CPM_tot=sum(CPM_taxa_rel)) %>%
  dplyr::group_by(lake,month_full,phylum) %>%
  dplyr::mutate(.,percent_CPM=CPM_taxa_rel/CPM_tot*100) %>%
  dplyr::group_by(phylum) %>%
  dplyr::summarise_at(vars("percent_CPM"),list(mean=mean, sd=sd))


N_tax$data %>% subset(gene =="ureC") %>%
  dplyr::group_by(lake,month_full) %>%
  dplyr::mutate(.,CPM_tot=sum(CPM_taxa_rel)) %>%
  dplyr::group_by(lake,month_full,phylum) %>%
  dplyr::mutate(.,percent_CPM=CPM_taxa_rel/CPM_tot*100) %>%
  dplyr::group_by(phylum) %>%
  dplyr::summarise_at(vars("percent_CPM"),list(mean=mean, sd=sd))

####__ Phosphorus BGC tax ####

P_tax<-
  all_fun.df %>% subset(month %in% c("A","B","C","D")) %>%
  cbind(.,gene=BGC_short$Level_2[match(.$KEGG_ko,BGC_short$Kegg_ko)]) %>%
  dplyr::group_by(lake,month) %>% dplyr::mutate(.,CPM_tot = sum(CPM)) %>%  # obtain total CPM
  dplyr::group_by(lake,month,gene,phylum) %>% dplyr::summarise(.,CPM_sum = sum((CPM/CPM_tot)*100)) %>% # merge phyla by gene and sample
  #dplyr::group_by(lake,month,gene,phylum) %>% dplyr::summarise(CPM_rel=(CPM_sum/CPM_tot)*100) %>% # merge phyla by gene and sample
  subset(gene %in% c("phnD","phnM","phoD","phoX","ppk1","ppx","pstS"))  %>%
  dplyr::mutate(.,gene=factor(gene,levels = c("ppk1","ppx","pstS","phnD","phnM","phoD","phoX")),
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
                                              "VSM","JAB"))),
                month_full=if_else(month=="A","Jun.",
                                    if_else(month=="B","Jul.",
                                            if_else(month=="C","Aug.","Sep."))),
                month_full=factor(month_full,
                                    levels=c("Jun.","Jul.","Aug.","Sep."))) %>%
  # dplyr::group_by(sample) %>% dplyr::mutate(.,CPM_tot = sum(CPM)) %>%
  # dplyr::group_by(sample,gene,phylum) %>% dplyr::mutate(.,CPM_sum = sum(CPM)) %>%
  # dplyr::group_by(lake,gene,phylum) %>% dplyr::summarise(CPM_median=median(CPM_sum)) %>%
  # dplyr::group_by(lake,gene) %>% dplyr::mutate(CPM_rel=((CPM_median/sum(CPM_median))*100)) %>%
  dplyr::mutate(.,
                phylum = if_else(phylum %in% c("Cyanobacteria","Bacteroidetes",
                                               "Actinobacteria","Proteobacteria"),phylum,"Others")) %>%
  dplyr::group_by(lake,month_full,gene,phylum) %>%  dplyr::summarise(CPM_taxa_rel = sum(CPM_sum)) %>%
  ggplot(.,aes(month_full,CPM_taxa_rel,fill=fct_rev(fct_reorder(phylum,CPM_taxa_rel))))+
  geom_col(width=0.8,color="black",size=0.1)+theme_bw()+
  theme(panel.grid = element_blank(),
        legend.title = element_text(face="bold",size=15),
        legend.text = element_text(size=12),
        axis.title = element_blank(),axis.ticks = element_blank(),
        axis.text.x = element_text(size=9), axis.text.y = element_text(size=9))+
  scale_fill_manual(values=c("#FED9A6","#DECBE4","firebrick","#B3CDE3","#CCEBC5"))+
  # "#F2F2F2","#000000","#E5D8BD","firebrick",
  # 
  # "#70F4BB","#5D9BEF","#F391CC","#EEF4A3","#CAAB7D","#8AF3D4","#EAF7B9","#4263B2","#E0B333","#2CB550",
  # "#BE6AE7","#8B7664","#7161EC","#CACCA8","#F138F4","#8EAC73"))+
  guides(fill=guide_legend(ncol = 1,override.aes = list(shape=22,size=5)))+
  scale_y_continuous(expand = c(0,0),position = 'right')+
  labs(fill="Prokaryotic Phyla",y="Relative abundance (CPM)")+
  facet_grid2(gene~lake,scales = "free",switch = "y",
              strip=
                strip_themed(background_x = elem_list_rect(fill = rev(palette_lake_3T),color = "black"),
                             text_x = elem_list_text(colour = "white",face = "bold",size=10),
                             background_y = elem_list_rect(fill = c("#5AAA46","#5AAA46","#5AAA46","#5AAA46",
                                                                    "#5AAA46","#5AAA46","#5AAA46"),color = "black"),
                             text_y = elem_list_text(colour = "white",face="bold.italic",size=14)))
P_tax

P_tax$data %>% subset(gene =="phnM") %>%
  dplyr::group_by(lake,month_full) %>%
  dplyr::mutate(.,CPM_tot=sum(CPM_taxa_rel)) %>%
  dplyr::mutate(.,Status=if_else(lake %in% c("VSS","CSM"),"Hypereutrophic",
                                 if_else(lake %in% c("JAB","VSM"),"Mesotrophic",
                                         "Eutrophic"))) %>%
  dplyr::group_by(Status) %>%
  dplyr::summarise_at(vars("CPM_tot"),list(mean=mean, sd=sd))

P_tax$data %>% subset(gene =="phnM") %>%
  dplyr::group_by(lake,month_full) %>%
  dplyr::mutate(.,CPM_tot=sum(CPM_taxa_rel)) %>%
  dplyr::mutate(.,Status=if_else(lake %in% c("JAB","VSM"),"Mesotrophic","others")) %>%
  dplyr::group_by(Status) %>%
  dplyr::summarise_at(vars("CPM_tot"),list(mean=mean, sd=sd))

P_tax$data %>% subset(gene =="phnM") %>%
  dplyr::group_by(lake,month_full) %>%
  dplyr::mutate(.,CPM_tot=sum(CPM_taxa_rel)) %>%
  dplyr::group_by(lake,month_full,phylum) %>%
  dplyr::mutate(.,percent_CPM=CPM_taxa_rel/CPM_tot*100) %>%
  dplyr::group_by(phylum) %>%
  dplyr::summarise_at(vars("percent_CPM"),list(mean=mean, sd=sd))

####__ Sulfur & Iron BGC tax ####

S_I_tax<-
  all_fun.df %>% subset(month %in% c("A","B","C","D")) %>%
  cbind(.,gene=BGC_short$Level_2[match(.$KEGG_ko,BGC_short$Kegg_ko)]) %>%
  dplyr::group_by(lake,month) %>% dplyr::mutate(.,CPM_tot = sum(CPM)) %>%  # obtain total CPM
  dplyr::group_by(lake,month,gene,phylum) %>% dplyr::summarise(.,CPM_sum = sum((CPM/CPM_tot)*100)) %>% # merge phyla by gene and sample
  subset(gene %in% c("asrB","aprA","dsrB","soxB","dmdA","fecA"))  %>%
  dplyr::mutate(.,gene=factor(gene,levels = c("asrB","aprA","dsrB","soxB","dmdA","fecA")) ,
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
                                              "VSM","JAB"))),
                month_full=if_else(month=="A","Jun.",
                                   if_else(month=="B","Jul.",
                                           if_else(month=="C","Aug.","Sep."))),
                month_full=factor(month_full,
                                  levels=c("Jun.","Jul.","Aug.","Sep."))) %>%
  # dplyr::group_by(sample) %>% dplyr::mutate(.,CPM_tot = sum(CPM)) %>%
  # dplyr::group_by(sample,gene,phylum) %>% dplyr::mutate(.,CPM_sum = sum(CPM)) %>%
  # dplyr::group_by(lake,gene,phylum) %>% dplyr::summarise(CPM_median=median(CPM_sum)) %>%
  # dplyr::group_by(lake,gene) %>% dplyr::mutate(CPM_rel=((CPM_median/sum(CPM_median))*100)) %>%
  dplyr::mutate(.,
                phylum = if_else(phylum %in% c("Bacteroidetes",
                                               "Actinobacteria","Proteobacteria"),phylum,"Others")) %>%
  dplyr::group_by(lake,month_full,gene,phylum) %>%  dplyr::summarise(CPM_taxa_rel = sum(CPM_sum)) %>%
  ggplot(.,aes(month_full,CPM_taxa_rel,fill=fct_rev(fct_reorder(phylum,CPM_taxa_rel))))+
  geom_col(width=0.8,color="black",size=0.1)+theme_bw()+
  theme(panel.grid = element_blank(),
        legend.title = element_text(face="bold",size=15),
        legend.text = element_text(size=12),
        axis.title = element_blank(),axis.ticks = element_blank(),
        axis.text.x = element_text(size=9),axis.text.y = element_text(size=9))+
  scale_fill_manual(values=c("#DECBE4","#B3CDE3","#FED9A6","firebrick"))+
  # "#F2F2F2","#000000","#E5D8BD","firebrick",
  # 
  # "#70F4BB","#5D9BEF","#F391CC","#EEF4A3","#CAAB7D","#8AF3D4","#EAF7B9","#4263B2","#E0B333","#2CB550",
  # "#BE6AE7","#8B7664","#7161EC","#CACCA8","#F138F4","#8EAC73"))+
  guides(fill=guide_legend(ncol = 1,override.aes = list(shape=22,size=5)))+
  scale_y_continuous(expand = c(0,0),position = 'right')+
  labs(fill="Microbial Phyla",y="Relative abundance (CPM)")+
  facet_grid2(gene~lake,scales = "free",switch = "y",
              strip=
                strip_themed(background_x = elem_list_rect(fill = rev(palette_lake_3T),color = "black"),
                             text_x = elem_list_text(colour = "white",face = "bold",size=10),
                             background_y = elem_list_rect(fill = c("#C84D4C","#C84D4C","#C84D4C","#C84D4C",
                                                                    "#C84D4C","black"),color = "black"),
                             text_y = elem_list_text(colour = "white",face="bold.italic",size=14)))
S_I_tax

S_I_tax$data %>% subset(gene =="soxB") %>%
  dplyr::group_by(lake,month_full) %>%
  dplyr::mutate(.,CPM_tot=sum(CPM_taxa_rel)) %>%
  dplyr::mutate(.,Status=if_else(lake %in% c("VSS","CSM"),"Hypereutrophic",
                                 if_else(lake %in% c("JAB","VSM"),"Mesotrophic",
                                         "Eutrophic"))) %>%
  dplyr::group_by(Status) %>%
  dplyr::summarise_at(vars("CPM_tot"),list(median=median, sd=sd))

S_I_tax$data %>% subset(gene =="dmdA") %>%
  dplyr::group_by(lake,month_full) %>%
  dplyr::mutate(.,CPM_tot=sum(CPM_taxa_rel)) %>%
  dplyr::mutate(.,Status=if_else(lake %in% c("JAB"),"Mesotrophic","other")) %>%
  dplyr::group_by(Status) %>%
  dplyr::summarise_at(vars("CPM_tot"),list(mean=mean, sd=sd))

S_I_tax$data %>%subset(gene %in% c("dsrB")) %>%
  dplyr::group_by(lake,month_full) %>%
  dplyr::mutate(.,CPM_tot=sum(CPM_taxa_rel)) %>%
  dplyr::group_by(lake,month_full,phylum) %>%
  dplyr::mutate(.,percent_CPM=CPM_taxa_rel/CPM_tot*100) %>%
  dplyr::group_by(phylum) %>%
  dplyr::summarise_at(vars("percent_CPM"),list(mean=mean, sd=sd))

S_I_tax$data %>% subset(gene %in% c("soxB","aprA","asrB","dsrB")) %>%
  dplyr::group_by(lake,month_full) %>%
  dplyr::mutate(.,CPM_tot=sum(CPM_taxa_rel)) %>%
  dplyr::group_by(lake,month_full,gene,phylum) %>%
  dplyr::mutate(.,percent_CPM=CPM_taxa_rel/CPM_tot*100) %>%
  dplyr::group_by(lake,gene,phylum) %>%
  dplyr::summarise_at(vars("percent_CPM"),list(mean=mean, sd=sd)) %>%
  cbind(.,gene_metabolism=BGC_short$Level_1[match(.$gene,BGC_short$Level_2)]) %>%
  subset(gene_metabolism =="S_metabolism") %>%
  dplyr::group_by(gene_metabolism,phylum) %>%
  dplyr::summarise_at(vars("mean"),list(mean1=mean, sd=sd))

####|####
