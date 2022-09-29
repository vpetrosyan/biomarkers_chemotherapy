args = commandArgs(trailingOnly=TRUE)
adj_matrix = args[1]
module_genes = args[2]
ind = as.numeric(args[3])

require(CTD)
require(igraph)

module_df = read.csv(module_genes,header=TRUE)
colnames(module_df) = c("genes","module")
module_df$genes = tolower(module_df$genes)

#Remove those in the grey module - unassigned by WGCNA 
module_df = module_df[module_df$module != "grey",]


adj_df = as.matrix(read.csv(adj_matrix,header=TRUE,row.names = 1))
rownames(adj_df) = tolower(rownames(adj_df))
colnames(adj_df) = tolower(colnames(adj_df))



adj_df = as.matrix(read.csv(adj_matrix,header=TRUE,row.names = 1))
rownames(adj_df) = tolower(rownames(adj_df))
colnames(adj_df) = tolower(colnames(adj_df))

ig = graph.adjacency(adj_df, mode="undirected", weighted=TRUE, add.colnames="name")
V(ig)$name = tolower(V(ig)$name)

G = vector(mode="list", length=length(V(ig)$name))
names(G) = V(ig)$name
adj_mat = as.matrix(get.adjacency(ig, attr="weight"))

p0=0.1
p1=0.9
thresholdDiff=0.01
all_nodes = tolower(V(ig)$name)
S_perturbed_nodes = unlist(module_df$genes)

current_node_set = singleNode.getNodeRanksN(ind,G,p1=0.9,thresholdDiff=0.01,
                                            adj_mat,S=S_perturbed_nodes,
                                            num.misses=log2(length(G)),TRUE)

names = sub(".*\\/", "", adj_matrix)
names = gsub("_5000pruned02.csv","",names)


save(current_node_set, file = sprintf("%s/%d_publish.RData", names, ind))



