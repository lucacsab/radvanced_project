##############################################################################
#                       Analyze miRNA - target networks                     #
##############################################################################

# 2021-03-24
# Luca Csabai

install.packages("igraph")
if (!requireNamespace("RCy3", quietly = TRUE)) 
  BiocManager::install("RCy3")

# Load libraries ----
library("tidyverse")
library("data.table")
library("igraph")
library("RCy3")



# Read data ----
network_raw <- read_tsv("input/Exp-miBRS_track_information_hg38.tsv")

# Functions from Benedek----
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
# Select the miRNA(s) and target gene(s) interactions only exonic:
network_raw <- network_raw[network_raw$REGION == 'exonic', ] 

int_data <- network_raw %>% 
  dplyr::select(c(MIRNA, GENES)) %>% 
  filter(!is.na(MIRNA))

# Process data, to have single, pairwise target gene - miRNA interactions:
int_df <- apply(int_data, 1, create_pairs_1)  
int_df <- do.call("rbind", int_df) %>% 
  dplyr::select(target_gene, mirnas)
int_df <- apply(int_df, 1, create_pairs_2)  
int_df <- do.call("rbind", int_df)

# Remove NA values that were in a list 
int_df[ int_df == "NA" ] <- NA
int_df <- int_df %>% drop_na()

# Separate nodes
target_gene <- int_df %>%
  distinct(target_gene) %>%
  rename(label = target_gene)

source_mirna <- int_df %>%
  distinct(mirna) %>%
  rename(label = mirna)

node_df <- full_join(source_mirna, target_gene, by = "label")

# Change direction 
int_df <- int_df[,c(2, 1)]

# Create an igraph object
net <- graph.data.frame(int_df, directed=T)

# Visualize with igraph
#V(net)$size <- 10

#l <- layout.kamada.kawai(net)
#l <- layout.norm(l, ymin=-1, ymax=1, xmin=-1, xmax=1)

#plot(net, rescale=F, layout=l*2, vertex.label=NA)

# Calculate node betweenness values
betweenness_res <- betweenness(net, v = V(net), directed = TRUE, weights = NULL,
                               nobigint = TRUE, normalized = FALSE)
# Put into dataframe
betweenness_res2 <- data.frame(betweenness_centrality=sort(betweenness_res, decreasing=TRUE)) %>% tibble::rownames_to_column()

# Add to node table
nodes_2 <- left_join(node_df, betweenness_res2, by = c("label"="rowname"))

##### Import into Cytoscape ##### ----

# Preprocess
net2 <- int_df %>% dplyr::rename("source"= "mirna", "target"="target_gene")
net2['interaction'] = 'inhibits'
nodes2 <- nodes_2 %>% dplyr::rename("id" = "label")

# Determine if betweenness centrality calcualtions have been carried out
if("betweenness_centrality" %in% colnames(nodes_2)){
  bet = TRUE
} else {
  bet=FALSE 
}

# See if cytoscape is open
cyto_error <- FALSE
tryCatch( {msg <- cytoscapePing () } , error = function(e) {cyto_error <<- TRUE})

if (!cyto_error){
  continue = TRUE
  message('Successfully connected to Cytoscape - continuing with visualisation')
} else {
  continue = FALSE
  message('Could not connect to Cytoscape - skipping visualisation')
}

# Run visualisation if cytoscape open
if(continue){
  createNetworkFromDataFrames(nodes2, edges=net2[,1:2], title="exonic network")
  edges <- net2 %>%
    mutate(key=paste(source, "(interacts with)", target))
  print(edges)
  loadTableData(edges[,3:4], data.key.column = 'key', table = 'edge')
  
  # Colour by Betweenness centrality if data exists in network
  if (bet){
    setNodeColorMapping("betweenness_centrality",mapping.type="c", network="exonic network")
  }
  
  # Save cys file
  saveSession(filename = file.path("output/exonic_network.cys"))
}

# Node degrees
# Compute node degrees (#links) and use that to set node size:
deg <- degree(net, mode="all")
V(net)$size <- deg*3
