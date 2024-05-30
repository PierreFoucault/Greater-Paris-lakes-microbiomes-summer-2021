#### Import data ####

sample_metadata<-read_delim("bin_metadata.csv",";", escape_double = FALSE,trim_ws = TRUE,show_col_types = FALSE) %>%
  dplyr::mutate(.,lake_order=factor(lake,levels=c("CHA","GDP","VER","CTL","BOI","CRJ1","CRJ2","JAB","VAI")))

all_fun.df<-rbind(
  M_A_1_fun.df,M_A_2_fun.df,M_A_3_fun.df,M_A_4_fun.df,M_A_5_fun.df,M_A_7_fun.df,M_A_8_fun.df,M_A_9_fun.df,
  M_B_1_fun.df,M_B_2_fun.df,M_B_3_fun.df,M_B_4_fun.df,M_B_5_fun.df,M_B_6_fun.df,M_B_7_fun.df,M_B_8_fun.df,M_B_9_fun.df,
  M_C_1_fun.df,M_C_2_fun.df,M_C_3_fun.df,M_C_4_fun.df,M_C_5_fun.df,M_C_6_fun.df,M_C_7_fun.df,M_C_8_fun.df,M_C_9_fun.df,
  M_D_1_fun.df,M_D_2_fun.df,M_D_3_fun.df,M_D_4_fun.df,M_D_5_fun.df,M_D_6_fun.df,M_D_7_fun.df,M_D_8_fun.df,M_D_9_fun.df) %>%
  cbind(.,
        lake=(sample_metadata$lake_order[match(.$sample,sample_metadata$sample)]),
        month=(sample_metadata$month[match(.$sample,sample_metadata$sample)]),
        lake_month=(sample_metadata$lake_month[match(.$sample,sample_metadata$sample)])) %>%
as.data.frame()

write.csv(all_fun.df,"all_fun.csv")
