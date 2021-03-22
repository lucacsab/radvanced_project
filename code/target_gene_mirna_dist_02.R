##############################################################################
# Check the connection between the target gene and regulatory miRNA distance #
##############################################################################

# 2021-02-22
# Benedek Dankó


# Load libraries ----
library("tidyverse")
library("data.table")
library("viridis")


# Read data ----
regulatory_data <- read_tsv("/home/dbenedek/uni/msc/Semester_4/Advanced_R_programming_for_biologists/project/data/Exp-miBRS_track_information_hg38.tsv")

# Load Gencode annotation, format:
gencode <- "/home/dbenedek/uni/msc/Semester_4/Advanced_R_programming_for_biologists/project/data/gencode.v37.annotation.gtf.gz"
gencode_annot <- fread(gencode,
                       col.names = c("chr", "ensemble", "type",
                                     "start", "end", "V6", "strand",
                                     "V8", "info")) %>% 
  mutate(gene_id = str_extract(info, "ENSG\\d+")) %>% 
  mutate(gene_name = str_extract(info, 'gene_name\\s.\\w+')) %>% 
  mutate(gene_name = str_replace(gene_name, 'gene_name \"', ""))

# Load chromosome info:
chr_info <- fread("/home/dbenedek/uni/msc/Semester_4/Advanced_R_programming_for_biologists/project/data/chromInfo.txt") %>% 
  dplyr::select(-V3)


# Functions ----
create_pairs_1 <- function(x){
  mirnas <- strsplit(x[1], ",") %>% 
    `[[`(1)
  target_gene <- x[2]
  df <- data.frame(mirnas=mirnas,
                   target_gene=rep(target_gene, length(mirnas)))
  return(df)
}

create_pairs_2 <- function(x){
  genes <- strsplit(x[1], ",") %>% 
    `[[`(1)
  mirna <- x[2]
  df <- data.frame(target_gene=genes,
                   mirna=rep(mirna, length(genes)))
  return(df)
}


# Process data ----
# Select the miRNA(s) and target gene(s) interactions:
regulatory_data_reduced <- regulatory_data %>% 
  dplyr::select(c(MIRNA, GENES)) %>% 
  filter(!is.na(MIRNA))

# Process data, to have single, pairwise target gene - miRNA interactions:
df_interactions <- apply(regulatory_data_reduced, 1, create_pairs_1)  
df_interactions <- do.call("rbind", df_interactions) %>% 
  dplyr::select(target_gene, mirnas)
df_interactions <- apply(df_interactions, 1, create_pairs_2)  
df_interactions <- do.call("rbind", df_interactions) %>% 
  distinct() %>% 
  mutate(mirna=str_replace(mirna, "miR", "mir")) %>% 
  mutate(mirna=str_replace(mirna, "hsa-", "")) %>% 
  mutate(mirna=str_replace(mirna, "-", "")) %>%
  mutate(mirna=str_replace(mirna, "mir", "MIR")) %>% 
  mutate(mirna=str_replace(mirna, "-\\dp", "")) %>% 
  mutate(mirna=str_replace(mirna, "let", "MIRLET")) %>% 
  mutate(mirna=toupper(mirna)) %>% 
  mutate(mirna=str_replace(mirna, "-", ""))


# Add Gencode annotation:
df_interactions_annot <- df_interactions %>% 
  left_join(gencode_annot, by=c("mirna"="gene_name")) %>% 
  dplyr::select(c(chr, target_gene, mirna, gene_id)) %>% 
  dplyr::rename(mirna_chr=chr,
                minra_gene_id = gene_id)
df_interactions_annot <- df_interactions_annot %>% 
  left_join(gencode_annot, by=c("target_gene"="gene_name")) %>% 
  dplyr::rename(target_gene_chr=chr,
                target_gene_id = gene_id) %>% 
  dplyr::select(-c(info, ensemble, type, start, end, V6, strand, V8)) %>% 
  filter(!is.na(target_gene_chr)) %>% 
  mutate(chr_interaction = paste0(mirna_chr, " - ",target_gene_chr))

chr_interactions <- df_interactions_annot %>% 
  dplyr::select(chr_interaction, target_gene_chr, mirna_chr) %>% 
  group_by(chr_interaction) %>% 
  mutate(count=n()) %>% 
  ungroup() %>% 
  arrange(desc(count)) %>% 
  distinct()

interactions_same_chr <- chr_interactions %>% 
  filter(mirna_chr == target_gene_chr) %>% 
  inner_join(chr_info, by=c("target_gene_chr"="V1")) %>% 
  dplyr::rename(chr_length=V2)


# Plot ----
# Plot the miRNA - taregt gene pairs where they are on the same chromosome:
barplot_same_chr <- ggplot(interactions_same_chr, 
                           aes(x=reorder(target_gene_chr, -count), y=count,
                               fill=chr_length))+
  scale_fill_viridis(discrete = F, name="Chr length",
                     breaks=c(50000000, 100000000, 150000000, 200000000),
                     labels=c("50 Mb", "100 Mb", "150 Mb", "200 Mb"))+
  geom_bar(stat="identity")+
  theme_minimal() +
  theme(panel.border = element_rect(color = "black", fill = NA),
        axis.text.x = element_text(size=12, hjust = 0.5, vjust = 0.5, angle=90),
        axis.title.x = element_text(size = 14, margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size = 14, margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.text.y = element_text(size = 12),
        legend.title = element_text(size=12))+
  xlab("Chromosome")+
  ylab("Interaction count")


# Save plot ----
ggsave("/home/dbenedek/uni/msc/Semester_4/Advanced_R_programming_for_biologists/project/docs/intrachr_mirna_tagetgene_counts.png",
       barplot_same_chr, units = "in", width = 12, height = 8)