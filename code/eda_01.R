# Load, format gencode annotation
# EDA


library("tidyverse")
library("data.table")
library("ggplot2")


# Load Gencode annotation, format:
gencode <- "/home/dbenedek/uni/msc/Semester_4/Advanced_R_programming_for_biologists/project/data/gencode.v37.annotation.gtf.gz"
gencode_annot <- fread(gencode,
                       col.names = c("chr", "ensemble", "type",
                                     "start", "end", "V6", "strand",
                                     "V8", "info")) %>% 
  mutate(gene_id = str_extract(info, "ENSG\\d+")) %>% 
  mutate(gene_name = str_extract(info, 'gene_name\\s.\\w+')) %>% 
  mutate(gene_name = str_replace(gene_name, 'gene_name \"', ""))


# Get only the miRNA data:
mirna_annot <- gencode_annot %>% 
  filter(grepl("miRNA", info)) # 5643 rows

# Count the number of miRNAs per chromosome:
mirna_per_chr <- mirna_annot %>% 
  group_by(chr) %>% 
  summarize(count=n()) %>% 
  mutate(chr=factor(chr, levels=c("chr1", "chr2","chr3","chr4","chr5",
                                  "chr6","chr7","chr8","chr9","chr10",
                                  "chr11","chr12","chr13","chr14","chr15",
                                  "chr16","chr17","chr18","chr19","chr20",
                                  "chr21","chr22","chrX","chrY")))


# Load chromosome info:
chr_info <- fread("/home/dbenedek/uni/msc/Semester_4/Advanced_R_programming_for_biologists/project/data/chromInfo.txt") %>%
  dplyr::select(-V3) %>%
  filter(V1 %in% mirna_per_chr$chr) %>%
  left_join(mirna_per_chr, by=c("V1"="chr"))

 
# Plot number of miRNAs per chromosome:
boxplot_mirnacount <- ggplot(mirna_per_chr,
                             aes(x=chr,
                                 y=count)) + 
  geom_bar(position="dodge", stat="identity") +
  scale_y_continuous(trans='log10')+
  theme_minimal() +
  theme(panel.border = element_rect(color = "black", fill = NA),
        axis.text.x = element_text(size=12, hjust = 1, vjust = 0.5, angle=90),
        axis.title.x = element_text(size = 13),
        axis.title.y = element_text(size = 13),
        axis.text.y = element_text(size = 12),
        legend.position="none") +
  xlab("Chromosome") + 
  ylab("miRNA count")


# Correlation test - length of chromsome vs. number of miRNAs:
cor.test(chr_info$V2, chr_info$count)
# Pearson's corr. = 0.52


# Plot the length of the chromosome vs. number of miRNAs:
cor_plot <- ggplot(chr_info, aes(x=V2, y=count,label = V1))+
  geom_point(size=3)+
  ggrepel::geom_text_repel()+
  theme_minimal() +
  theme(panel.border = element_rect(color = "black", fill = NA),
        axis.text.x = element_text(size=12, hjust = 1, vjust = 0.5, angle=0),
        axis.title.x = element_text(size = 14, margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size = 14, margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.text.y = element_text(size = 12),
        legend.position="none")+
  xlab("Chromosome length (bp)")+
  ylab("Number of miRNAs")


# Save plot:
ggsave("/home/dbenedek/uni/msc/Semester_4/Advanced_R_programming_for_biologists/project/docs/corr_plot.png",
       cor_plot, units = "in", width = 14, height = 10)
