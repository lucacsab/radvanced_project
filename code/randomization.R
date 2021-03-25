# Randomization

# Data cleaning

# Loading the source table
raw_data <- read.table("Exp-miBRS_track_information_hg38.tsv")
# Colum names from the first row
names(raw_data) <- raw_data[1,]
# Deleting the first row
raw_data <- raw_data[-1,]
# Deleting NA rows
miRNA_data <- subset(raw_data,subset=!is.na(MIRNA))
# Deleting useless informations
miRNA_data_min <- miRNA_data[,c(1,10,12)]
# Separating the miRNAs to rows for each regulated gene
library(tidyr)
# miRNA_data_2 <- separate(miRNA_data_min,MIRNA,c("MIRNA1","MIRNA2","MIRNA3","MIRNA4","MIRNA5","MIRNA6","MIRNA7","MIRNA8","MIRNA9","MIRNA10",
#                                             "MIRNA11","MIRNA12","MIRNA13","MIRNA14","MIRNA15","MIRNA16","MIRNA17","MIRNA18","MIRNA19","MIRNA20",
#                                             "MIRNA21","MIRNA22","MIRNA23","MIRNA24","MIRNA25","MIRNA26","MIRNA27","MIRNA28","MIRNA29","MIRNA30",
#                                             "MIRNA31","MIRAN32","MIRNA33","MIRNA34","MIRNA35","MIRNA36","MIRNA37","MIRNA38","MIRNA39","MIRNA40",
#                                             "MIRNA41","MIRNA42","MIRNA43","MIRNA44","MIRNA45","MIRNA46","MIRNA47","MIRNA48","MIRNA49","MIRNA50",
#                                             "MIRNA51","MIRAN52","MIRNA53","MIRNA54","MIRNA55","MIRNA56","MIRNA57","MIRNA58","MIRNA59","MIRNA60",
#                                             "MIRNA61","MIRNA62","MIRNA63","MIRNA64","MIRNA65","MIRNA66","MIRNA67","MIRNA68","MIRNA69","MIRNA70",
#                                             "MIRNA71","MIRNA72"),sep=",")
miRNA_data_3 <- separate_rows(miRNA_data_min,MIRNA,sep=",")
miRNA_data_4 <- subset(miRNA_data_3,subset=MIRNA!="NA")

# Randomizing
miRNA_data_randomized <- miRNA_data_4
miRNA_data_randomized$CHROM <- sample(miRNA_data_4$CHROM,
                                size=nrow(miRNA_data_4),
                                replace=FALSE)