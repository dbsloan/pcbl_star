args <- commandArgs(TRUE)

if (length(args) != 7) stop ("Required command line arugments are: input_tree_file  outgroup_name  output_name  start_tree_num  end_tree_num  batch_num  spp_tree_file")

library(phybase, warn.conflicts=FALSE)

#define output file names by adding various extensions to the provided base names
output_trees = paste (args[3], "_trees_", args[6], ".txt", sep="")
output_pcbl = paste (args[3], "_pcbl_", args[6], ".txt", sep="")
output_presence_absence = paste (args[3], "_node_presence_", args[6], ".txt", sep="")

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

#read in previously generated species tree
full_sptree_obj = read.tree(args[7])
full_sptree = write.tree(root(full_sptree_obj, args[2], resolve.root=TRUE))

#extract node information from full species tree
full_sptree_nodes = read.tree.nodes(full_sptree)

#define number of trees based on provided start and end indexes
start_tree_num = as.numeric(args[4]);
end_tree_num = as.numeric(args[5]);

num_gene_trees = end_tree_num - start_tree_num + 1;

#define empty matrix to store pcbl scores
pcbl_mat = matrix(, nrow=num_gene_trees, ncol=length(genetrees$names)-3)

#define empty matrix to store node presence/absence info. This will include a 1 if a node is still in a species tree after removing a given gene tree and a 0 if not.
presence_absence_mat = matrix(, nrow=num_gene_trees, ncol=length(genetrees$names)-3)

sink(output_trees)
#loop through each gene tree
for (gt in start_tree_num:end_tree_num){
	#generate new species tree based on all gene trees except the current one
	sub_sptree=star.sptree(genetrees$tree[-gt], speciesname=genetrees$names, taxaname=genetrees$names,species.structure=species.structure, outgroup=args[2], method="nj")
	
	#skipping trees if they phybase fails to produce a species tree. This can occur if two taxa only co-occur in a single tree (so removing that tree results in undefined distances).
	sub_sptree_obj = read.tree(text=sub_sptree)
	if (class(sub_sptree_obj) != "phylo"){
		cat ("Warning: skipping analysis of gene tree ", gt, " because no species tree was produced after removing this tree.\n", sep="");
		cat (sub_sptree, "\n")
		pcbl_mat[gt - start_tree_num + 1, 1] = "Skipped"
		presence_absence_mat[gt - start_tree_num + 1, 1] = "Skipped"
		next
	}

	cat (sub_sptree, "\n")
	sub_sptree_nodes = read.tree.nodes(sub_sptree)
	internal_node_count = 0
	#loop through all nodes in the full tree
	for (node in 1:dim(full_sptree_nodes$nodes)[1]){
		#get a vector of species (by index number) that are descendents of that node
		full_node_indexes = offspring.species(node,full_sptree_nodes$nodes,numtax)
		#skip the root node, the internal node that includes everything except the outgroup, or indvidual terminal nodes (i.e., number of taxa = all species, all species minus 1, or just 1 species)
		if (length(full_node_indexes) == 1 || length(full_node_indexes) >= numtax - 1) next
		internal_node_count = internal_node_count + 1;
		#convert species to actual names
		full_node_names = c();
		for (index_counter in 1:length(full_node_indexes)){
			full_node_names[index_counter] = full_sptree_nodes$names[full_node_indexes[index_counter]]
		}
		#sort species name vector
		sorted_full_node_names = sort (full_node_names)
		#loop through all nodes in the subtree to see if one exists with the same species content
		match_node=0
		for (sub_node in 1:dim(sub_sptree_nodes$nodes)[1]){
			sub_node_indexes = offspring.species(sub_node,sub_sptree_nodes$nodes,numtax)
			sub_node_names = c();
			for (sub_index_counter in 1:length(sub_node_indexes)){
				sub_node_names[sub_index_counter] = sub_sptree_nodes$names[sub_node_indexes[sub_index_counter]]
			}
			sorted_sub_node_names = sort (sub_node_names)
			if (identical(sorted_sub_node_names, sorted_full_node_names)){
				if(match_node) stop ("ERROR: multiple nodes in the sub_tree were identified as a match")
				match_node=sub_node
			}
		}
		#calculate PCBL by substracting branch length from this leave-one-out tree from the corresponding branch length in the "full" tree. If the current tree does not have the node from the full tree, use a branch length of 0 in the calculation. Store PCBL values in marix 
		sub_bl = 0
		present = 0
		if(match_node){ 
			sub_bl = sub_sptree_nodes$nodes[match_node,4]
			present = 1
		}
		pcbl_mat[gt - start_tree_num + 1, internal_node_count] = full_sptree_nodes$nodes[node,4] - sub_bl
		presence_absence_mat[gt - start_tree_num + 1, internal_node_count] = present
		
	}
}
sink()

#print PCBL matrix to output file
write.table(pcbl_mat, file=output_pcbl, sep ="\t", row.names=FALSE)

#print presence/absence matrix to output file
write.table(presence_absence_mat, file=output_presence_absence, sep ="\t", row.names=FALSE)