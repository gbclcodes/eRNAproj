#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use File::Path;

my ($motiffile,$outdir,$help);

GetOptions(
	"motiffile=s" => \$motiffile,
	"outdir=s" => \$outdir,
	"help!" => \$help,
);

open(MOT,"<$motiffile") or die "$!\n";
my $header="MEME version 4

ALPHABET= ACGT

strands: + -

Background letter frequencies
A 0.25 C 0.25 G 0.25 T 0.25

";

my @array = <MOT>;
foreach my $element (@array){
	chomp $element;
	if($element =~ /MOTIF/){
		my @fieldValues = split /\s+/,$element;
		open(OUT,">$outdir/$fieldValues[1].meme") or die "$!\n";
		print OUT "$header\n";
		print OUT "$element\n";
	}else{
		print OUT "$element\n";
	}
}
close MOT;
close OUT;

# perl step2.pl --motiffile /media/yuhua/yuhua_projects/enhProj/annodata/JASPAR_CORE_non-redundant_pfms_meme.txt --outdir /media/yuhua/yuhua_projects/enhProj/annodata/JASPAR_CORE

