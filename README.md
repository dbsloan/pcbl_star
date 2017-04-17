# pcbl_star
Scripts for automating Partitioned Coalescence Branch Length (PCBL) analysis using STAR

pcbl_star.pl

Overview: 
This script implements the "partitioned coalescence branch length" method by automating the process of removing one tree at a time from a set of gene trees and recalculating a species tree using the STAR method (Phybase: http://faculty.franklin.uga.edu/lliu/content/phybase). New species trees are then compared against a species tree inferred from the full set of gene trees to determine changes in branch length that resulted from removing the gene tree in question. Changes in subtending branch lengths (PCBL values) are reported for every clade in the species tree. A summary of which clades are "lost" (i.e., are no longer monophyletic) after removing each gene is also reported.



Requirements: 

This automation is implemented with Perl and R scripts. They have been designed for a Unix environment (Mac OSX or Linux). They have been tested in Mac OSX 10.11 and Linux CentOS 6, but they should work in most Unix environments.

Perl - The provided Perl script should called by users (pcbl_star.pl). It in turn calls the R code. Perl is pre-installed in most Mac OSX and Linux distributions.

R - The STAR method is part of the R package Phybase (http://faculty.franklin.uga.edu/lliu/content/phybase). Therefore, R and Phybase must be installed. This also requires the installation of additional R packages: ape, Matrix, and methods. The Rscript application that allows running R code must be either in your PATH or you must provide a full path to it when calling pcbl_star.pl (see below). The scripts have been tested on Phybase 1.4 and 1.5, but would likely work on other versions as well.



Running pcbl_star.pl:
The scripts can be called from the command line to analyze a set of gene trees with STAR. There are a mix of required and optional parameters that are specified at the command line as described below. Sample data and expected output files are provided in the sample_data directory.


Usage: perl pcbl_star.pl [arguments]

ARGUMENTS

Gene trees must be specified with one of the following options (but
not both): 

--gt_file      - a single file containing one or more trees

--gt_dir       - a directory containing one or more files each 
containing a single tree 

Path and filename for R code to run STAR:

--r_code        - default: pcbl_star.R

Path and filename for the Rscript application installed on local machine:

--Rscript       - default: Rscript

Outgroup species:

--outgroup        - name of outgroup taxon (must be a single species, not a clade)

Base name for output files:

--output        - name for output files


EXAMPLE
perl pcbl_star.pl --gt_file=sample_data/sample_trees_file.tre --outgroup=Danio --output=star_sample

