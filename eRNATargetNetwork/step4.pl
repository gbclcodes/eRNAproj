#!/usr/bin/perl -w

use strict;
use Getopt::Long;
use File::Path;

my ($coexprfile,$posfile,$outfile,$help);
GetOptions(
	"coexprfile=s" => \$coexprfile,
	"posfile=s" => \$posfile,
	"outfile=s" => \$outfile,
	"help!" => \$help,
);

my (%coHash,%posHash);
open(CO,"<$coexprfile") or die "$!\n";
while(<CO>){
	my $line = $_;
	chomp $line;
	my @fieldValues = split /\s+/,$line;
	$coHash{"$fieldValues[0]\t$fieldValues[1]"} = "$fieldValues[2]\t$fieldValues[3]\t$fieldValues[4]\t$fieldValues[5]"; 
}
close CO;

open(POS,"<$posfile") or die "$!\n";
while(<POS>){
	my $line = $_;
	chomp $line;
	my @fieldValues = split /\s+/,$line;
	$fieldValues[1] =~ s/\.\d+//g;
	$posHash{"$fieldValues[0]\t$fieldValues[1]"} = "$fieldValues[2]"; 
}
close POS;

open(OUT,">$outfile") or die "$!\n";
foreach my $key (keys %coHash){
	if(exists $posHash{$key}){
		print OUT "$key\t$coHash{$key}\t$posHash{$key}\n";
	}
}
close OUT;

# perl fig3a_step4.pl --coexprfile /media/yuhua/yuhua_projects/enhProj/ENHData/enh_gene_coexpr_pairs_GSR_all.txt --posfile /media/yuhua/yuhua_projects/enhProj/ENHData/enh_gene_within_5mb_GSR.txt --outfile /media/yuhua/yuhua_projects/enhProj/ENHData/enh_target_pairs_GSR_all.txt

# perl fig3a_step4.pl --coexprfile /media/yuhua/yuhua_projects/enhProj/ENHData/enh_gene_coexpr_pairs_XW_all.txt --posfile /media/yuhua/yuhua_projects/enhProj/ENHData/enh_gene_within_5mb_XW.txt --outfile /media/yuhua/yuhua_projects/enhProj/ENHData/enh_target_pairs_XW_all.txt