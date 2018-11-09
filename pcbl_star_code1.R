args <- commandArgs(TRUE)

if (length(args) != 3) stop ("Required command line arugments are: input_tree_file outgroup_name output_name")

library(phybase, warn.conflicts=FALSE)

#define output file names by adding various extensions to the provided base names
spp_tree = paste (args[3], "_species-tree.txt", sep="")
output_nodes = paste (args[3], "_node_defs.txt", sep="")

genetrees<-read.tree.string(file=args[1],format="phylip")

#read.tree.string only extract taxa names from the first tree. Add any additional names from other trees 
for (input_tree in genetrees$tree){
	input_tree_nodes = read.tree.nodes(input_tree)
	for (input_name in input_tree_nodes$names){
		if(!input_name %in% genetrees$names){
			genetrees$names = append(genetrees$names, input_name)
		}
	}
}

#define a species matrix. This assumes that there is one terminal for each species. STAR can accomodate mutliple terminals per species, but the following code would have to be modified to specify that in the matrix.
numtax=length(genetrees$names)
species.structure<-matrix(0,numtax,numtax)
diag(species.structure)<-1

#generate species tree with STAR based on all gene trees. Note this requires specification of a single terminal as an outgroup for rooting the resulting trees. This is specificied by the user on the command line 
full_sptree = star.sptree(genetrees$tree, speciesname=genetrees$names, taxaname=genetrees$names,species.structure=species.structure, outgroup=args[2], method="nj")

#print species tree to output file
sink(spp_tree)
cat (full_sptree)
sink()

#extract node information from full species tree
full_sptree_nodes = read.tree.nodes(full_sptree)

#print definition of nodes to output file
sink(output_nodes)
internal_node_count = 0
#loop through each node, and obtain all taxa names that are descendents of that node
for (node in 1:dim(full_sptree_nodes$nodes)[1]){
	
	full_node_indexes = offspring.species(node,full_sptree_nodes$nodes,numtax)
	if (length(full_node_indexes) == 1 || length(full_node_indexes) >= numtax - 1) next
	internal_node_count = internal_node_count + 1
	full_node_names = c();
	for (index_counter in 1:length(full_node_indexes)){
		full_node_names[index_counter] = full_sptree_nodes$names[full_node_indexes[index_counter]]
	}	
	cat (internal_node_count, sort(full_node_names), "\n")
}
sink()