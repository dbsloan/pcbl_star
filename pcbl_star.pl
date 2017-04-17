#!/usr/bin/perl

use strict;
use warnings;
use List::Util qw (min);
use Getopt::Long;
use File::Which;

our $gt_file;
our $gt_dir;
our $r_code = "pcbl_star.R";
our $Rscript = "Rscript";
our $outgroup;
our $output;


my $usage = 
"\nUsage: perl $0 [arguments]
   
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
         perl $0 --gt_file=sample_data/sample_trees_file.tre --outgroup=Danio --output=star_sample
         	
\n\n";

GetOptions(
    'gt_file=s'  => \$gt_file,
    'gt_dir=s'  => \$gt_dir,
    'r_code=s'  => \$r_code,
    'Rscript=s'  => \$Rscript,
    'outgroup=s'  => \$outgroup,
    'output=s'  => \$output,
);

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

-e $r_code or die ("\nERROR: $r_code does not exist. Provide correct path to the provided R code for running star with --r_code. \n\n$usage");

which $Rscript or die ("\nERROR: $Rscript does not exist. Provide correct path to the Rscript application on your system with --Rscript. \n\n$usage");

#Define random number that will be used in saving temp file names to avoid clashing if multiple runs are done in parallel.
my $randnum = int(rand(10000));

print STDERR "\n\n". (localtime) ."\nRunning PCBL using STAR as implemented in phybase\n\n";
my $gt_source;
if ($gt_file){$gt_source = $gt_file;}else{$gt_source = $gt_dir;}
print STDERR "Analyzing gene trees in $gt_source using $outgroup as an outgroup.\n\n\tCheck progress in $output\_trees.txt\n\n";



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

system ("$Rscript $r_code $gt_file $outgroup $output");


#if $gt_file was created as temp file from $gt_dir, then delete it (but delete user provided file).
$gt_dir and unlink($gt_file);


rewrite_output ("$output\_pcbl.txt", $gt_dir);
rewrite_output ("$output\_node_presence.txt", $gt_dir);

print STDERR "\n" . (localtime) . "\nRun Completed\n\n";


###define subroutines###

sub rewrite_output {

	use strict;
	use warnings;

    my ($file, $dir) = @_;
    
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
    
    my $gt_count = 0;
    
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