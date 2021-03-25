#################
# Randomization #
#################

# Sophie Farkas

# Data cleaning ----

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
miRNA_data_3 <- separate_rows(miRNA_data_min,MIRNA,sep=",")
miRNA_data_4 <- subset(miRNA_data_3,subset=MIRNA!="NA")

# Randomizing ----

miRNA_data_randomized <- miRNA_data_4
# Randomizing the target gene's location
miRNA_data_randomized$CHROM <- sample(miRNA_data_4$CHROM,
                                size=nrow(miRNA_data_4),
                                replace=FALSE)
# Randomizing the target gene
miRNA_data_randomized$GENES <- sample(miRNA_data_4$GENES,
                                size=nrow(miRNA_data_4),
                                replace=FALSE)
