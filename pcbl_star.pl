#!/usr/bin/perl

use strict;
use warnings;
use List::Util qw (min);
use Getopt::Long;
use File::Which;

our $gt_file;
our $gt_dir;
our $r_code1 = "pcbl_star_code1.R";
our $r_code2 = "pcbl_star_code2.R";
our $Rscript = "Rscript";
our $outgroup;
our $output;
our $threads = 1;
our $batch_size = 10;
our $retain_subs;
our $start_tree = 1;
our $end_tree;


my $usage = 
"\nUsage: perl $0 [arguments]
   
   REQUIRED ARGUMENTS
   
   Gene trees must be specified with one of the following options (but
   not both): 

      --gt_file      - a single file containing one or more trees

      --gt_dir       - a directory containing one or more files each containing
                       a single tree 

   Outgroup species:
   
      --outgroup     - name of outgroup taxon (must be a single species, not a
                       clade). Rooting of the output species tree is required
                       for data processing, but the position of the root does
                       not affect PCBL calculations for any bipartitions within
                       the tree. Therefore, it is OK to arbitrarily select one
                       species if the actual outgroup is a clade.

   Base name for output files:
   
      --output       - name for output files


   OPTIONAL ARGUMENTS

   Multithreading:

      --threads      - number of parallel batches of trees to run [default: 1]

      --batch_size   - number of trees per batch if multithreading [default: 10]

      --retain_subs  - Add this flag if individual files from batches should be 
                       kept after merging into final output files.


   Subsetting jobs by tree numbers:
   
      --start_tree   - tree number to start with if only doing PCBL calculations
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

      perl $0
         --gt_file=sample_data/sample_trees_file.tre
         --outgroup=Danio
         --output=star_sample
         --threads=5
         --batch_size=6  	
\n\n";

GetOptions(
    'gt_file=s'  => \$gt_file,
    'gt_dir=s'  => \$gt_dir,
    'r_code1=s'  => \$r_code1,
    'r_code1=s'  => \$r_code1,
    'Rscript=s'  => \$Rscript,
    'outgroup=s'  => \$outgroup,
    'output=s'  => \$output,
    'threads=i'  => \$threads,
    'batch_size=i'  => \$batch_size,
    'retain_subs'  => \$retain_subs,
    'start_tree=i'  => \$start_tree,
    'end_tree=i'  => \$end_tree
);


if ($threads > 1){
	use Parallel::ForkManager;
}


#Checked that required command line arguments were provided
$output or die ("\nERROR: Must specify base name for all output files with --output\n\n$usage");
$outgroup or die ("\nERROR: Must specify base name for all output files with --output\n\n$usage");
$gt_file or $gt_dir or die ("\nERROR: Must specify gene trees with either --gt_file or --gt_dir\n\n$usage");
$gt_file and $gt_dir and die ("\nERROR: --gt_file and --gt_dir cannot both be used at the same time\n\n$usage");
if ($gt_file){
	-d $gt_file and die ("\nERROR: the --gt_file option was used but $gt_file is a directory. You may want to use --gt_dir instead.\n\n$usage");
}else{
	-d $gt_dir or die ("\nERROR: the --gt_dir option was used but $gt_dir is NOT a directory. You may want to use --gt_file instead.\n\n$usage");
}

-e $r_code1 or die ("\nERROR: $r_code1 does not exist. Provide correct path to the provided R code for running star with --r_code1. \n\n$usage");
-e $r_code2 or die ("\nERROR: $r_code2 does not exist. Provide correct path to the provided R code for running star with --r_code2. \n\n$usage");

which $Rscript or die ("\nERROR: $Rscript does not exist. Provide correct path to the Rscript application on your system with --Rscript. \n\n$usage");

#Define random number that will be used in saving temp file names to avoid clashing if multiple runs are done in parallel.
my $randnum = int(rand(10000));

print STDERR "\n\n". (localtime) ."\nRunning PCBL using STAR as implemented in phybase\n\n";
my $gt_source;
if ($gt_file){$gt_source = $gt_file;}else{$gt_source = $gt_dir;}
print STDERR (localtime) . "\nAnalyzing gene trees in $gt_source using $outgroup as an outgroup to build full STAR species tree...\n\n";



#concatenate files from $gt_dir if user provided a directory of gene trees
if ($gt_dir){
	substr ($gt_dir, -1) eq '/' or $gt_dir .= '/';
	my @files = get_file_names ($gt_dir);
	my $FH_GT_OUT = open_output (".$randnum\_genetrees.txt");
	$gt_file = ".$randnum\_genetrees.txt";
	foreach (@files){
		my $FH_GT = open_file($gt_dir.$_);
		my $first_line = <$FH_GT>;
		chomp $first_line;
		print $FH_GT_OUT "$first_line\n";
	}
}

system ("$Rscript $r_code1 $gt_file $outgroup $output");

#Unless user has specified a subset of trees to run, get a total number of (non-empty) lines in tree file and set that as the end tree number
unless ($end_tree){
	my @arr = file_to_array($gt_file);
	foreach (@arr){
		$_ =~ /^\s*$/ or ++$end_tree;
	}
}

print STDERR "\n". (localtime) . "\nGenerating PCBL data by rerunning STAR after excluding one gene tree at a time...\n\n\tCheck progress in $output\_trees file(s).\n\n\tNumber of specified threads: $threads\n\n";

#keep track of how many batches were run (will only be 1 if there is no multithreading)
my $subfile_count = 0;

if ($threads == 1){
	$subfile_count = 1;
	system ("$Rscript $r_code2 $gt_file $outgroup $output $start_tree $end_tree $subfile_count $output\_species-tree.txt");
}else{
	#if multithreading is specified, use ForkManager to fork batches of trees to child processes
	my $pm = Parallel::ForkManager->new($threads);
	
	for (my $i = $start_tree; $i <= $end_tree; $i += $batch_size){
		my $batch_end = min ($end_tree, $i + $batch_size - 1);
		++$subfile_count;	
		my $pid = $pm->start and next; #this is the forking process that will start multiple R calls in parallel
		system ("$Rscript $r_code2 $gt_file $outgroup $output $i $batch_end $subfile_count $output\_species-tree.txt");
		$pm->finish; # end of child processes
	}
	$pm->wait_all_children;
}


#if $gt_file was created as temp file from $gt_dir, then delete it (but delete user provided file).
$gt_dir and unlink($gt_file);

my $FH_TREES = open_output ("$output\_trees.txt");
print $FH_TREES "STAR trees based on removal of one gene tree at a time...\n";
my $FH_NP = open_output ("$output\_node_presence.txt");
my $FH_NP_1 = open_file ("$output\_node_presence_1.txt");
my $np_first_line = <$FH_NP_1>;
close $FH_NP_1;
print $FH_NP $np_first_line;
my $FH_PCBL = open_output ("$output\_pcbl.txt");
my $FH_PCBL_1 = open_file ("$output\_pcbl_1.txt");
my $pcbl_first_line = <$FH_PCBL_1>;
close $FH_PCBL_1;
print $FH_PCBL $pcbl_first_line;

#loop through all subfiles for combining
for (my $i = 1; $i <= $subfile_count; ++$i){
	#tree files
	my $tree_file_data = file_to_string ("$output\_trees_$i\.txt");
	print $FH_TREES $tree_file_data;
	$retain_subs or unlink ("$output\_trees_$i\.txt");
	
	#node_presence files
	my @np_file_lines = file_to_array ("$output\_node_presence_$i\.txt");
	shift @np_file_lines;
	foreach (@np_file_lines){
		print $FH_NP $_;
	}
	$retain_subs or unlink ("$output\_node_presence_$i\.txt");
	
	#pcbl files
	my @pcbl_file_lines = file_to_array ("$output\_pcbl_$i\.txt");
	shift @pcbl_file_lines;
	foreach (@pcbl_file_lines){
		print $FH_PCBL $_;
	}
	$retain_subs or unlink ("$output\_pcbl_$i\.txt");
}

close $FH_NP;
close $FH_TREES;
close $FH_PCBL;

#clean-up formating from raw R table output
rewrite_output ("$output\_pcbl.txt", $gt_dir, $start_tree);
rewrite_output ("$output\_node_presence.txt", $gt_dir, $start_tree);

print STDERR "\n" . (localtime) . "\nRun Completed\n\n";


###define subroutines###

sub rewrite_output {

	use strict;
	use warnings;

    my ($file, $dir, $count_offset) = @_;
    
    my @filenames;
    
    $dir and @filenames = get_file_names ($dir);
    
    my $FH_RW = open_file ($file);
    
    my $output_string = "Removed Tree";
    
    my $first_line = <$FH_RW>;
    my @sfl = split (/\t/, $first_line);
    my $node_total = scalar (@sfl);
    
    for (my $i = 1; $i <= $node_total; ++$i){
    	$output_string .= "\tNode $i";
    }
    $output_string .= "\n";
    
    my $gt_count = $count_offset - 1;
    
    while (<$FH_RW>){
    	++$gt_count;
    	if (@filenames){
    		$output_string .= "$filenames[$gt_count-1]\t$_";
    	}else{
    		$output_string .= "Tree $gt_count\t$_";
    	}
    }
    
    close $FH_RW;
    my $FH_RWO = open_output ($file);
    
    print $FH_RWO $output_string;
	close $FH_RWO;

}

sub get_file_names {
	use strict;
	use warnings;

    my ($directory) = @_;
    my @files = (  );
    my @filedata =(  );
    	

    # Open the directory
    unless(opendir(DIRECTORY, $directory)) {
        print "Cannot open directory $directory!\n";
        exit;
    }
    
    # Read the directory, ignoring special entries starting with "."
    @files = grep (!/^\./, readdir(DIRECTORY));
    
    closedir(DIRECTORY);
    
    return (@files);
   
}

sub file_to_string {
	use strict;
	use warnings;

    my($filename) = @_;

    my $filedata;

    unless( open(GET_FILE_DATA, $filename) ) {
        print STDERR "Cannot open file \"$filename\"\n\n";
        exit;
    }
    
    while (<GET_FILE_DATA>){
    	$filedata .= $_;
    }
    
    close GET_FILE_DATA;

    return $filedata;
}

sub file_to_array {
	use strict;
	use warnings;

    my($filename) = @_;

    my @filedata = (  );

    unless( open(GET_FILE_DATA, $filename) ) {
        print STDERR "Cannot open file \"$filename\"\n\n";
        exit;
    }

    @filedata = <GET_FILE_DATA>;

    close GET_FILE_DATA;

    return @filedata;
}


sub open_file {
	use strict;
	use warnings;

    my($filename) = @_;
    my $fh;

    unless(open($fh, $filename)) {
        print "Cannot open file $filename\n";
        exit;
    }
    return $fh;
}

sub open_output {
	use strict;
	use warnings;

    my($filename) = @_;
    my $fh_output;

    unless(open($fh_output, ">$filename")) {
        print "Cannot open file $filename\n";
        exit;
    }
    return $fh_output;
}