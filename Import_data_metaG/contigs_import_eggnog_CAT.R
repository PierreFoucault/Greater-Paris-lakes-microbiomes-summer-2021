#### Import packages ####

library(BiocManager)
library(tidyverse)
library(ggpubr) 
library(tidytext)

####. ####

#### Data table format ####

n_threads=8
import_eggnog <- function (file,col.names,col.type) {
  col_names = col.names
  read_tsv(file, col_names = col_names, comment = "#",show_col_types = F,num_threads = n_threads,col_types = col.type,)
}

# col names and types for emapper seed_orthologs file
col.ortholog<-c("qseqid","egg_sid","egg_evalue","egg_bitscore","egg_qstart","egg_qend","egg_sstart","egg_send","egg_pident","egg_qcov","egg_scov")
col.type.ortholog<-cols(.default = "d",qseqid = "character",egg_sid = "character")

# col names and types for emapper diamond hits file
col.diamond <- c("qseqid","dmd_sid","dmd_pident","dmd_qlength","dmd_TBD","dmd_mismatch","dmd_qstart","dmd_qend",
                 "dmd_sstart","dmd_send","dmd_evalue","dmd_bitscore","dmd_qcov","dmd_scov") 
col.type.diamond<-cols(.default = "d",qseqid = "character",dmd_sid = "character")

# col names and types for emapper annotation file
col.annotation<-c("qseqid","seed_ortholog","evalue","score","eggNOG_OGs","max_annot_lvl","COG_category","Description","Preferred_name",
                  "GOs","EC","KEGG_ko","KEGG_Pathway","KEGG_Module","KEGG_Reaction","KEGG_rclass","BRITE","KEGG_TC","CAZy","BiGG_Reaction","PFAMs")
col.type.annotation<-cols(.default = "character",evalue = "d",score = "d")

####. ####

#### Import data ####

#### M_A ####

####__ M_A_1 ####
sample="M_A_1"

M_A_1.cat<-read_delim("/media/duperron/Pierre_PhD/CAT/M_A_1.cat.names.txt","\t",escape_double = F,trim_ws = T,show_col_types = F,num_) %>%
  dplyr::rename("contig" = "# contig") %>% drop_na(lineage) %>% .[,c(1,6:ncol(.))] %>%
  dplyr::filter(., superkingdom %in% c("Bacteria: 1.00","Archaea: 1.00")) #View(M_A_1.cat)

M_A_1.contigs.counts<-read_tsv("/media/duperron/Pierre_PhD/CAT_quant/M_A_1.quant.counts",num_threads = 8,show_col_types = F)

M_A_1_fun.df <- list(M_A_1.annotation = import_eggnog("/media/duperron/Pierre_PhD/emapper_contigs_data/M_A_1.emapper.annotations",col.annotation,col.type.annotation) %>%
                       dplyr::mutate(.,contig_info = gsub('_[^_]*$', '', qseqid),coverage = gsub(".*\\_", "",contig_info),coverage=as.numeric(coverage),.after = 1) %>%
                       separate(KEGG_ko,c("KEGG_ko", "KEGG_ko_grouped"), ",",fill = "right",extra = "merge") %>%
                       separate(eggNOG_OGs,c("root_OG","kingdom_OG","phylum_OG","class_OG","tax_OG"), ",",fill = "right",extra = "merge") %>%
                       dplyr::mutate(.,kingdom_OG = gsub(".*\\|", '',kingdom_OG),phylum_OG = gsub(".*\\|", '',phylum_OG),class_OG = gsub(".*\\|", '',class_OG)) %>%
                       .[,c(1:6,8:10,14:15,18:19)],
                     M_A_1.diamond.hits = import_eggnog("/media/duperron/Pierre_PhD/emapper_contigs_data/M_A_1.emapper.hits",col.diamond,col.type.diamond) %>% .[,c(1:3,11:14)],
                     M_A_1.ortholog = import_eggnog("/media/duperron/Pierre_PhD/emapper_contigs_data/M_A_1.emapper.seed_orthologs",col.ortholog,col.type.ortholog) %>% .[,c(1:4,9:11)]) %>%
  reduce(left_join,"qseqid") %>% dplyr::mutate(.,sample = sample,.before = 1) %>% dplyr::filter(., dmd_qcov >= 80 | egg_qcov >= 80) %>%
  dplyr::mutate(.,
                CPM = M_A_1.contigs.counts$count[match(.$contig_info,M_A_1.contigs.counts$transcript)],
                level_1 = kegg_id$Level_1[match(.$KEGG_ko,kegg_id$Kegg_ko)],level_2 = kegg_id$Level_2[match(.$KEGG_ko,kegg_id$Kegg_ko)],
                level_3 = kegg_id$Level_3[match(.$KEGG_ko,kegg_id$Kegg_ko)],
                domain_score = M_A_1.cat$superkingdom[match(.$contig_info,M_A_1.cat$contig)],domain = if_else(domain_score %in% c("no support",NA),kingdom_OG,gsub(':[^:]*$','',domain_score)),
                phylum_score = M_A_1.cat$phylum[match(.$contig_info,M_A_1.cat$contig)],phylum = if_else(phylum_score %in% c("no support",NA),phylum_OG,gsub(':[^:]*$', '',phylum_score)),
                class_score = M_A_1.cat$class[match(.$contig_info,M_A_1.cat$contig)],class = if_else(class_score %in% c("no support",NA),class_OG,gsub(':[^:]*$', '',class_score))) %>%
  subset(.,KEGG_ko != "-") %>% group_by(KEGG_ko,sample,domain,phylum,class) %>%
  summarize(nb_sample_tax =  n(),cov_sample_tax = sum(coverage),CPM = sum(CPM))
#View(M_A_1_fun.df)
