# pcbl_star
Scripts for automating Partitioned Coalescence Branch Length (PCBL) analysis using STAR

`pcbl_star.pl`

## Overview: 
This script implements the "partitioned coalescence branch length" method by automating the process of removing one tree at a time from a set of gene trees and recalculating a species tree using the STAR method in [Phybase](http://faculty.franklin.uga.edu/lliu/phybase). New species trees are then compared against a species tree inferred from the full set of gene trees to determine changes in branch length that resulted from removing the gene tree in question. Changes in subtending branch lengths (PCBL values) are reported for every clade in the species tree. A summary of which clades are "lost" (i.e., are no longer monophyletic) after removing each gene is also reported.



## Requirements: 

This automation is implemented with Perl and R scripts. They have been designed for a Unix environment (Mac OSX or Linux). They have been tested in Mac OSX 10.11 and Linux CentOS 6, but they should work in most Unix environments.

Perl - The provided Perl script should be called by users (`pcbl_star.pl`). It in turn calls the accompanying R code (`pcbl_star_code1.R` and `pcbl_star_code2.R`). Perl is pre-installed in most Mac OSX and Linux distributions. If using the `--threads` option to process batches of trees in parallel, the [Parallel::ForkManager](https://metacpan.org/pod/Parallel::ForkManager) Perl module must be installed.

R - The STAR method is part of the R package Phybase [Phybase](http://faculty.franklin.uga.edu/lliu/phybase). Therefore, R and the Phybase package must be installed. This also requires the installation of additional R packages: `ape`, `Matrix`, and `methods`. The `Rscript` application that allows running R code must be either in your PATH or you must provide a full path to it when calling `pcbl_star.pl` (see below). The scripts have been tested on Phybase 1.4 and 1.5. **They do not currently work with Phybase v2.0**.



## Running pcbl_star.pl:
The scripts can be called from the command line to analyze a set of gene trees with STAR. There are a mix of required and optional parameters that are specified at the command line as described below. Sample data and expected output files are provided in the sample_data directory.


Usage: `perl pcbl_star.pl [arguments]`
   
   REQUIRED ARGUMENTS
   
   Gene trees must be specified with one of the following options (but
   not both): 

      --gt_file      - a single file containing one or more trees

      --gt_dir       - a directory containing one or more files each containing
                       a single tree 

   Outgroup species:
   
      --outgroup     - name of outgroup taxon (must be a single species, not a
                       clade)

   Base name for output files:
   
      --output       - name for output files


   OPTIONAL ARGUMENTS

   Multithreading:

      --threads      - number of parallel batches of trees to run [default: 1]

      --batch_size   - number of trees per batch if multithreading [default: 10]

      --retain_subs  - Add this flag if individual files from batches should be 
                       kept after merging into final output files.


   Subsetting jobs by tree numbers:
   
      --start_tree   - tree number to end with if only doing PCBL calculations
                       on a subset of trees. Note that reference species tree
                       will still be inferred from all gene trees. This is used
                       to split up an analysis across multiple different
                       machines. [default: 1]
		
      --end_tree     - tree number to end with if only doing PCBL calculations
                       on a subset of trees. Note that reference species tree
                       will still be inferred from all gene trees. This is used
                       to split up an analysis across multiple different
                       machines. [default: MAX]


   Path and filenames for R code to run STAR:
   
      --r_code1      - default: pcbl_star_code1.R

      --r_code2      - default: pcbl_star_code2.R


   Path and filename for the Rscript application installed on local machine:
   
      --Rscript      - default: Rscript

   
   EXAMPLE

      perl pcbl_star.pl
         --gt_file=sample_data/sample_trees_file.tre
         --outgroup=Danio
         --output=star_sample
         --threads=5
         --batch_size=6  	


## Output Files

#### Output_pcbl.txt
A tab-delimited text file, reporting PCBL scores at each node for each tree in the dataset. Positive values indicate support for a node such that excluding that gene tree reduces the branch length supporting the node in the inferred species tree.

#### Output_species-tree.txt
A text file with the inferred species tree based on all gene trees (newick format).

#### Output_trees.txt
A text file with the species trees inferred after removing one gene tree at a time in order (newick format). There is one tree per line. The first line is a note about the file contents. Warnings about excluded trees will also appear in this file (for example, if removing a gene tree prevents building a species tree because it removes the only tree in which a specific pair of taxa co-occur).

#### Output\_node_presence.txt
A tab-delimited text file, indicating which nodes from the full species-tree are retained (1) or lost (0) after removing each gene tree and reanalyzing.

#### Output\_node_defs.txt
A text file defining the node numbering scheme used in `Output_pcbl.txt` and `Output_node_presence.txt`. The taxa descended from each node are listed.