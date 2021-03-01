#!/usr/bin/perl -w

use strict;
use Getopt::Long qw(GetOptions);

#Brooke Wolford
#Feb 10, 2014
#
#Vivek Rai
#Oct 23, 2017
#CHANGELOG
# - Don't append `chr` when CHRM column contains chr already

# get command line options:
my ($h, $vcf, $het, $id,$snp,$indel,$keepDupHet);
GetOptions(
    'h'   => \$h,
    'v=s' => \$vcf,
    't'   => \$het,
    'i=s' => \$id,
    's' => \$snp,
    'd' => \$indel,
    'k' => \$keepDupHet,
);


my $usage = <<USAGE;

    USAGE: perl vcf2bed.pl -v <vcf> -i <ID> -t -s -d -k 

    options:
    -h display this help message
    -v vcf file name 
    -t only print data for het sites (optional)
    -i id you want (optional) (default is all)
    -s if you want to exclude indels (optional) (default is SNPs and indels)
    -d if you want to exclude SNPs (optional) (default is SNPs and indels)
    -k if the vcf has two entries for the same position, the heterozgyous entry
       will be kept (if both are het, neither will be printed) (if both are
       homozygous, neither will be printed)

    NOTE*
    There are rows at beginning starting with pound
    columns are chr pos(1 based) ID Ref Alt Qual Filter Info Format IDs
    genotype is formatted like 1|1:100:44 ie only phased genotypes
    output is formatted chr, posS, posE, alleleA/hap1, alleleB/hap2, refAllele, ID
USAGE

# check command line arguments and display help if needed
if ($h) { die "$usage"; }
unless ($vcf) { die "$usage"; }

# open vcf and cycle through line by line:
if ( $vcf =~ /\.gz$/ ) {
    open( VCF, "gunzip -c $vcf |" );
}
else {
    open( VCF, $vcf );
}

#retrieve sample IDs from file 
my @ID_array = &fetchID($vcf);


#give specific ID to print out genotypes for (default is to print all)
if ($id){
    my $column = &idColumn($id,@ID_array);
    if ($column == 0) {
        if ($ID_array[0] ne $id) {
            die "This ID is not in VCF\n";
        }    
    }
    &readVCF_col($column);
}
else{
    &readVCF_all;
}

# #####################
# #### SUBROUTINES ####
# #####################

sub fetchID {
	my @ID_array=();
	while(my $line = <VCF>) {
	chomp $line;
	my ($chr, $pos, $ID, $ref, $alt, $qual, $filter, $info, $format, @alleles)=split /\s+/, $line;
		if ($chr eq '#CHROM'){
			@ID_array=@alleles;
			return @ID_array;
		}
	}
}
#what column number has the info for the id given
sub idColumn {
    my $id_name = $_[0];
    my @array = @_[1..$#_];
    for (my $i=0; $i <= $#array; $i++){
	if ($array[$i] eq $id_name) {
	    return $i;
	}
    }
}

sub readVCF_col {
    my $col_num = shift;
    my $lastPos;
    my $lastGenotype;
    my $outputLine;
    while(my $line = <VCF>) {
	chomp $line;
	my ($chr, $pos, $ID, $ref, $alt, $qual, $filter, $info, $format, @alleles)=split /\s+/, $line;
	    if ($chr !~ /\#/){
	    my ($alleleChunk, @junk) = split (":", $alleles[$col_num]);
	    if ($alleleChunk =~ /\|/){
		my ($alleleA, $alleleB) = split ('\|',$alleleChunk);
		$alleleA =~ s/1/$alt/g;
		$alleleA  =~ s/0/$ref/g;
		$alleleB =~ s/1/$alt/g;
		$alleleB =~ s/0/$ref/g;
		if ($lastPos && $lastPos==$pos){
		    if ($lastGenotype ne $alleleA."\t".$alleleB){
			warn "Duplicate genotypes at $chr:$pos $lastGenotype\t$alleleA\t$alleleB\n";
			if ($keepDupHet) {
			    if ($lastGenotype =~ /^(.)\t\1/) {
				$outputLine = "";
				warn "Both duplicate genotypes are homozygous\n";
			    }
			    if ($alleleA eq $alleleB) {
				next;
			    }
			    if ($outputLine){
				$outputLine = "";
				warn "Both duplicate genotypes are heterozygous\n";
				next;
			    }
			}
			else {
			    $outputLine="";
			    next;
			}
		    }
		    else {
			next;
		    }
		}
		if ($outputLine) {
		    print $outputLine."\n";
		}
		$outputLine="";
		$lastPos=$pos;
		$lastGenotype=$alleleA."\t".$alleleB;
		if ($het && $snp){
		    if ($alleleA ne $alleleB && length($alleleA)==1 && length($alleleB)==1 && length($ref)==1) {
			if ($chr =~ /chr/){
			    $outputLine=join("\t", $chr, $pos-1, $pos, $alleleA, $alleleB, $ref, $ID_array[$col_num]);
			}
			else { 
			    $outputLine=join("\t", "chr".$chr, $pos-1, $pos, $alleleA, $alleleB, $ref, $ID_array[$col_num]);
			}
			
		    }
		}
		elsif ($het && $indel) {
		    if ($alleleA ne $alleleB && (length($alleleA)>1 || length($alleleB)>1 || length($ref)>1)) {
			if ($chr =~ /chr/){
                            $outputLine=join("\t", $chr, $pos-1, $pos, $alleleA, $alleleB, $ref, $ID_array[$col_num]);
			}
			else {
			    $outputLine=join("\t", "chr".$chr, $pos-1, $pos, $alleleA, $alleleB, $ref, $ID_array[$col_num]);
			}
                    }
                }
		elsif ($het) {
		    if ($alleleA ne $alleleB) {
			if ($chr =~ /chr/){
                            $outputLine=join("\t", $chr, $pos-1, $pos, $alleleA, $alleleB, $ref, $ID_array[$col_num]);
			}
			else {
			    $outputLine=join("\t", "chr".$chr, $pos-1, $pos, $alleleA, $alleleB, $ref, $ID_array[$col_num]);
			}
		    }
		}
		elsif ($indel) {
		    if (length($alleleA)>1 || length($alleleB)>1 || length($ref)>1){
			if ($chr =~ /chr/){
                            $outputLine=join("\t", $chr, $pos-1, $pos, $alleleA, $alleleB, $ref, $ID_array[$col_num]);
			}
			else {
			    $outputLine=join("\t", "chr".$chr, $pos-1, $pos, $alleleA, $alleleB, $ref, $ID_array[$col_num]);
			}
		    }
		}
		elsif ($snp) {
		    if (length($alleleA)==1 && length($alleleB)==1 && length($ref)==1){
			if ($chr =~ /chr/){
                            $outputLine=join("\t", $chr, $pos-1, $pos, $alleleA, $alleleB, $ref, $ID_array[$col_num]);
			}
			else {
			    $outputLine=join("\t", "chr".$chr, $pos-1, $pos, $alleleA, $alleleB, $ref, $ID_array[$col_num]);
			}
		    }
		}
		else {
		    if ($chr =~ /chr/){
			$outputLine=join("\t", $chr, $pos-1, $pos, $alleleA, $alleleB, $ref, $ID_array[$col_num]);
		    }
		    else {
			$outputLine=join("\t", "chr".$chr, $pos-1, $pos, $alleleA, $alleleB, $ref, $ID_array[$col_num]);
		    }
		}
	    }
	}
    }
    if ($outputLine){
	print $outputLine."\n";
    }
}

sub readVCF_all {
    while(my $line = <VCF>) {
        chomp $line;
        my ($chr, $pos, $ID, $ref, $alt, $qual, $filter, $info, $format, @alleles)=split /\s+/, $line;
        if ($chr !~ /\#/){
	    for (my $i=0; $i <= $#alleles; $i++){
		    #break apart 0|1 format and replace numbers with letters
		my ($alleleChunk, $sep_two, $dosageChunk)=split(":", $alleles[$i]);
		if ($alleleChunk =~ /\|/){
		    my ($alleleA, $sep, $alleleB)=split("|", $alleleChunk);
		    $alleleA =~ s/1/$alt/g;
		    $alleleA  =~ s/0/$ref/g;
		    $alleleB =~ s/1/$alt/g;
		    $alleleB =~ s/0/$ref/g;
                    if ($chr !~ /chr/) {
                      $chr = "chr".$chr;
                    }
		    if ($het && $snp){
			if ($alleleA ne $alleleB && length($alleleA)==1 && length($alleleB)==1 && length($ref)==1) {
			    print(join("\t", $chr, $pos-1, $pos, $alleleA, $alleleB, $ref, $ID_array[$i]) . "\n");
			}
		    }
		    elsif ($het && $indel){
			if ($alleleA ne $alleleB && (length($alleleA)>1 || length($alleleB)>1 || length($ref)>1)) {
			    print(join("\t", $chr, $pos-1, $pos, $alleleA, $alleleB, $ref, $ID_array[$i]) . "\n");
                        }
                    }
		    elsif ($indel){
			if (length($alleleA)>1 || length($alleleB)>1 || length($ref)>1){
			    print(join("\t", $chr, $pos-1, $pos, $alleleA, $alleleB, $ref, $ID_array[$i]) . "\n");
			}
		    }
		    elsif ($het) {
			if ($alleleA ne $alleleB) {
			    print(join("\t", $chr, $pos-1, $pos, $alleleA, $alleleB, $ref, $ID_array[$i]) . "\n");
			}
		    }
		    elsif ($snp) {
			if (length($alleleA)==1 && length($alleleB)==1 && length($ref)==1){
			    print(join("\t", $chr, $pos-1, $pos, $alleleA, $alleleB, $ref, $ID_array[$i]) . "\n");
			}
		    }
		    else {
			print(join("\t", $chr, $pos-1, $pos, $alleleA, $alleleB, $ref, $ID_array[$i]) . "\n");
		    }
		}
	    }
	}
    }
}
