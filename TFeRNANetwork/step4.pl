#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use File::Path;

my ($indir,$outfile,$help);
GetOptions(
	"indir=s" => \$indir,
    "outfile=s" => \$outfile,
	"help!" => \$help,
);

open(OUT,">$outfile") or die "$!\n";
my @memefiles = `find $indir -name "*meme"`;
foreach my $memefile (@memefiles){
	open(IN,"<$memefile") or die "$!\n";
	while(<IN>){
		my $line = $_;
		chomp $line;
		if($line =~ /^MOTIF/){
			my @fieldValues = split /\s+/,$line;
			my @parts = split /\./,$fieldValues[2];
			print OUT "$fieldValues[1]\t$parts[2]\n";
		}
	}
	close IN;
}
close OUT;

# perl step4.pl --indir /media/yuhua/yuhua_projects/enhProj/annodata/JASPAR_CORE/ --outfile /media/yuhua/yuhua_projects/enhProj/annodata/JASPAR_CORE_ID_to_NAME.txt
