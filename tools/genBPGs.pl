#!/usr/bin/perl -w
#
# Patient Network LPIA-style = PNL (read: panel)
# Usage: perl PNL.pl gene_score.file GO_sets.file PATH_sets.file [up/down/abs] > <output_name>

### UPDATE 06-10-2015 ###
# desc:     create a patient network, through integration of posterior probability of a differential effect for a gene, GO terms and pathway information
# author:   Maciej M Kandula (maciej.kandula@gmail.com)

use strict;
use Data::Dumper;

open( SCORES, "$ARGV[0]" ) || die ( "Cannot open file $ARGV[0] for input\n" );
open( GO_YGENE, "$ARGV[1]" ) || die ( "Cannot open file $ARGV[3]\n" );
open( PATH_YGENE, "$ARGV[2]" ) || die ( "Cannot open $ARGV[4] for output\n" );

######## Set Up Variables, Expression Score Hash, and Gene Set Hashes/Arrays #########
my $direction = "abs"; # default; direction is either up, down, or abs
if (defined($ARGV[3])) {
    $direction =  "$ARGV[3]";
}

# output file names
# MK: redirecting the output to stdout; for automatic naming see the commented out code far below

# score hashes and other global variables
my %score_hash = ();
my $norm_score = 0;
my @row_names = ();
my $s;

# read in gene scores and save in score hash
while ( my $line = <SCORES> ) {
	chomp($line);
	$line =~ tr/[a-z]/[A-Z]/;
	my @tmp = split( /\t/, $line);
	if ( $direction =~ /up/ ){
		$s = exp($tmp[1]);
	}
	elsif ( $direction =~ /down/ ){
		$s = exp(-1 * $tmp[1]);
	}
	else {
		$s = abs($tmp[1]);
	}
	my @ids = split( /\s\/\/\/\s/, $tmp[0] );
	foreach my $i (@ids){
		if ( exists $score_hash{$i} ){
			# if exists in the hash (i.e. if this gene is repeated) we will take the max score
			if ( $score_hash{$i} < $s ){
				$score_hash{$i} = $s;
			}	
		}
		else {
			$score_hash{$i} = $s;
		}
	}

}

# Set up the Pathway:Genes hash 
# Set up an array of Pathway names
my %pathway_hash = ();
my @path_list = ();
while ( my $line = <PATH_YGENE> ) {
	chomp( $line );
	$line =~ tr/[a-z]/[A-Z]/;
	my ( $one, $rest ) = split( /\s+/, $line, 2);
	my @array = split( /\s+/, $rest );
	@{$pathway_hash{ $one }} = split(/\s+/, $rest);
	push(@path_list, $one);
}

# Set up the GO:Genes hash 
# Set up an array of GO names
my @go_list = ();
my %go_hash = ();
while ( my $line = <GO_YGENE> ) {
	chomp( $line );
	$line =~ tr/[a-z]/[A-Z]/;
	my ( $one, $rest ) = split( /\s+/, $line, 2);
	@{$go_hash{ $one }} = split(/\s+/, $rest);
	push(@go_list, $one);
}


# number of go terms
my $go_size = scalar(@go_list);

# number of pathways
my $path_size = scalar(@path_list);

#print STDERR "GO LIST: $go_size\n";
#print STDERR "PATH LIST: $path_size\n";

######## Main Script #########

my $median = 0;

# loop through the GO terms
for( my $g = 0; $g < $go_size; $g++ ){
    
    # loop through the PATH terms
    for( my $p = 0; $p < $path_size; $p++ ){
	
	# Compute the jaccard index and the median gene expression score in the intersection
	if ( exists $go_list[$g] && exists $path_list[$p] ) {
	    if (exists $go_hash{$go_list[$g]} && exists $pathway_hash{$path_list[$p]} ){
		my @go_genes = @{ $go_hash{$go_list[$g]} };
		my @path_genes = @{ $pathway_hash{$path_list[$p]} };
		
		# get intersection of the two arrays
		my %path_genes=map{$_ =>1} @path_genes;
		my %go_genes=map{$_=>1} @go_genes;
		my @int = grep( $path_genes{$_}, @go_genes );
		
		# get the union of the two arrays
		my %union = ();
		foreach(@path_genes,@go_genes){
		    $union{$_}=1;
		}
		my @union2 = keys %union;
		
		# get the sizes of the union and intersection
		my $union_size = scalar(@union2);						
		my $int_size = scalar(@int);
		
		my @score_list = ();
		# check to see if there is an intersection...
		if ( $int_size > 0 ){
		    
		    # loop through the intersecting genes
		    foreach my $i (@int) {
			if ( exists $score_hash{$i} ){
			    push( @score_list, $score_hash{$i} );
			}
		    }
		    my @sorted_int = sort { $a <=> $b } @score_list;
		    # get the number of gene scores, and compute the median!
		    my $num = scalar(@sorted_int);
		    if ( $num > 0 ) {
			if ($num%2==1) {
			    $median = $sorted_int[($#sorted_int / 2)];
			}
			else {
			    $median = $sorted_int[int($#sorted_int / 2)] + (($sorted_int[int($#sorted_int / 2) + 1] - $sorted_int[int($#sorted_int / 2)]) / 2);
			}
			$norm_score = $median * ($int_size/$union_size);
		    }
		    else { $norm_score = 0 };
		}
		else {
		    $norm_score = $int_size;
		}
		if ( $p == 0 ){
		    print STDOUT "$norm_score";
		}
		elsif ($p == $path_size-1) {
		    print STDOUT "\t$norm_score\n";
		}
		else {
		    print STDOUT "\t$norm_score";
		}
	    }
	}
    }
}

########################################################################
### Version of the tool that finds an output file name automatically ###
########################################################################

# comment out the above and uncomment the below to use the version that automatically sets the output file name

## output file names
#my $base =  "$ARGV[0]";
#my $frontbase = $base;
#$frontbase =~ s/.*\///; # this is to deal with input files from a different directory
#my $outputfile = "PNL_$frontbase";
#
## score hashes and other global variables
#my %score_hash = ();
#my $norm_score = 0;
#my @row_names = ();
#my $s;
#
## read in gene scores and save in score hash
#while ( my $line = <SCORES> ) {
#	chomp($line);
#	$line =~ tr/[a-z]/[A-Z]/;
#	my @tmp = split( /\t/, $line);
#	if ( $direction =~ /up/ ){
#		$s = exp($tmp[1]);
#	}
#	elsif ( $direction =~ /down/ ){
#		$s = exp(-1 * $tmp[1]);
#	}
#	else {
#		$s = abs($tmp[1]);
#	}
#	my @ids = split( /\s\/\/\/\s/, $tmp[0] );
#	foreach my $i (@ids){
#		if ( exists $score_hash{$i} ){
#			# if exists in the hash (i.e. if this gene is repeated) we will take the max score
#			if ( $score_hash{$i} < $s ){
#				$score_hash{$i} = $s;
#			}	
#		}
#		else {
#			$score_hash{$i} = $s;
#		}
#	}
#
#}
#
## Set up the Pathway:Genes hash 
## Set up an array of Pathway names
#my %pathway_hash = ();
#my @path_list = ();
#while ( my $line = <PATH_YGENE> ) {
#	chomp( $line );
#	$line =~ tr/[a-z]/[A-Z]/;
#	my ( $one, $rest ) = split( /\s+/, $line, 2);
#	my @array = split( /\s+/, $rest );
#	@{$pathway_hash{ $one }} = split(/\s+/, $rest);
#	push(@path_list, $one);
#}
#
## Set up the GO:Genes hash 
## Set up an array of GO names
#my @go_list = ();
#my %go_hash = ();
#while ( my $line = <GO_YGENE> ) {
#	chomp( $line );
#	$line =~ tr/[a-z]/[A-Z]/;
#	my ( $one, $rest ) = split( /\s+/, $line, 2);
#	@{$go_hash{ $one }} = split(/\s+/, $rest);
#	push(@go_list, $one);
#}
#
#
## number of go terms
#my $go_size = scalar(@go_list);
#
## number of pathways
#my $path_size = scalar(@path_list);
#
##print STDERR "GO LIST: $go_size\n";
##print STDERR "PATH LIST: $path_size\n";
#
######### Main Script #########
#
#my $median = 0;
#my $weighted_matrix = $outputfile;
#open ( ORIG, ">$weighted_matrix" ) || die ( "Could not open $weighted_matrix for output \n" );
#
## loop through the GO terms
#for( my $g = 0; $g < $go_size; $g++ ){
#    
#    # loop through the PATH terms
#    for( my $p = 0; $p < $path_size; $p++ ){
#	
#	# Compute the jaccard index and the median gene expression score in the intersection
#	if ( exists $go_list[$g] && exists $path_list[$p] ) {
#	    if (exists $go_hash{$go_list[$g]} && exists $pathway_hash{$path_list[$p]} ){
#		my @go_genes = @{ $go_hash{$go_list[$g]} };
#		my @path_genes = @{ $pathway_hash{$path_list[$p]} };
#		
#		# get intersection of the two arrays
#		my %path_genes=map{$_ =>1} @path_genes;
#		my %go_genes=map{$_=>1} @go_genes;
#		my @int = grep( $path_genes{$_}, @go_genes );
#		
#		# get the union of the two arrays
#		my %union = ();
#		foreach(@path_genes,@go_genes){
#		    $union{$_}=1;
#		}
#		my @union2 = keys %union;
#		
#		# get the sizes of the union and intersection
#		my $union_size = scalar(@union2);						
#		my $int_size = scalar(@int);
#		
#		my @score_list = ();
#		# check to see if there is an intersection...
#		if ( $int_size > 0 ){
#		    
#		    # loop through the intersecting genes
#		    foreach my $i (@int) {
#			if ( exists $score_hash{$i} ){
#			    push( @score_list, $score_hash{$i} );
#			}
#		    }
#		    my @sorted_int = sort { $a <=> $b } @score_list;
#		    # get the number of gene scores, and compute the median!
#		    my $num = scalar(@sorted_int);
#		    if ( $num > 0 ) {
#			if ($num%2==1) {
#			    $median = $sorted_int[($#sorted_int / 2)];
#			}
#			else {
#			    $median = $sorted_int[int($#sorted_int / 2)] + (($sorted_int[int($#sorted_int / 2) + 1] - $sorted_int[int($#sorted_int / 2)]) / 2);
#			}
#			$norm_score = $median * ($int_size/$union_size);
#		    }
#		    else { $norm_score = 0 };
#		}
#		else {
#		    $norm_score = $int_size;
#		}
#		if ( $p == 0 ){
#		    print ORIG "$norm_score";
#		}
#		elsif ($p == $path_size-1) {
#		    print ORIG "\t$norm_score\n";
#		}
#		else {
#		    print ORIG "\t$norm_score";
#		}
#	    }
#	}
#    }
#}
#
#close(ORIG);
