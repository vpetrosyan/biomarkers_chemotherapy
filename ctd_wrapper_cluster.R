args = commandArgs(trailingOnly=TRUE)
adj_matrix = args[1]
module_genes = args[2]


require(igraph)
require(dplyr)

module_df = read.csv(module_genes,header=TRUE)
colnames(module_df) = c("genes","module")
module_df$genes = tolower(module_df$genes)


#Remove those in the grey module - unassigned by WGCNA
module_df = module_df[module_df$module != "grey",]
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
print(length(S_perturbed_nodes))

for (s_node in S_perturbed_nodes){
  if (!(s_node %in% rownames(adj_mat))){
    stop(paste('Node ', s_node, ' not in graph. Exiting program.'))
  }}


names = sub(".*\\/", "", adj_matrix)
names = gsub("_5000pruned02.csv","",names)


for (n in 1:length(S_perturbed_nodes)) {
  print(n)
  ind = which(names(G)==S_perturbed_nodes[n])
  print(ind)
  print(sprintf("qsub -v adj_matrix=%s, module_genes=%s,ind=%d getPerm_genes_test.pbs", adj_matrix, module_genes, ind))
  str = sprintf("qsub -v adj_matrix=%s,module_genes=%s,ind=%d getPerm_genes_test.pbs", adj_matrix, module_genes, ind)
  system(str, wait=FALSE)
  system("sleep 0.2")
}




#Count the number of gene perms that have been run
ready=FALSE
last_sum = 0
while (!ready) {
  f = list.files(pattern="publish.RData")
  print(sprintf("#permutation files ready = %d", length(f)))
  if (length(f)==length(S_perturbed_nodes)) {
    ready=TRUE
  } else {
    system("sleep 20")
    system("rm *.pbs.*")
    curr_sum = length(f)
    if (curr_sum > last_sum) {
      last_sum = curr_sum
    }
  }
}
system("rm *.pbs.*")

# collate permutations into one R object


permutationByStartNode = list()
for (n in 1:length(S_perturbed_nodes)) {
  load(sprintf("%s/%d_publish.RData", names, ind))
  permutationByStartNode[[n]] = current_node_set
}
names(permutationByStartNode) = S_perturbed_nodes

save.image(file= sprintf("%s_all_publish.RData", names))
print("collate complete...")

