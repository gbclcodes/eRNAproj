#!/usr/bin/perl -w

use strict;
use Getopt::Long;
use File::Path;

my ($enhfile,$genefile,$outfile,$help);
GetOptions(
	"enhfile=s" => \$enhfile,
	"genefile=s" => \$genefile,
	"outfile=s" => \$outfile,
	"help!" => \$help,
);

sub min{
    my $mn = $_[0];
    for my $e(@_) {$mn = $e if ($e < $mn);}
    return $mn;
}

my %enhHash;
open(ENH,"<$enhfile") or die "$!\n";
while(<ENH>){
    my $line = $_;
    chomp $line;
    my @fieldValues = split /\t/,$line;
    $enhHash{$fieldValues[3]}{"chrom"} = $fieldValues[0];
    $enhHash{$fieldValues[3]}{"start"} = $fieldValues[1];
    $enhHash{$fieldValues[3]}{"end"} = $fieldValues[2];
}
close ENH;

my (%geneposHash,%genetypeHash);
open(GENE,"<$genefile") or die "$!\n";
while(<GENE>){
    my $line = $_;
    chomp $line;
    next if($line =~ /^#/);
    my @fieldValues = split /\t/,$line;
    if($fieldValues[2] eq "gene"){
        $fieldValues[8] =~ /gene_id=(.*?);gene_type=(.*?);/;
        my $geneid = $1;
		$geneid =~ s/\.\d+//g;
        my $genetype = $2;
        $genetypeHash{$geneid} = $genetype;
        $geneposHash{$fieldValues[0]}{$geneid}{"start"} = $fieldValues[3];
        $geneposHash{$fieldValues[0]}{$geneid}{"end"} = $fieldValues[4];
    }
}
close GENE;

open(OUT,">$outfile") or die "$!\n";
foreach my $enhid (keys %enhHash){
    my $enhchrom = $enhHash{$enhid}{"chrom"};
    my $enhstart = $enhHash{$enhid}{"start"};
    my $enhend = $enhHash{$enhid}{"end"};
    my @geneids = keys %{$geneposHash{$enhchrom}};
    foreach my $geneid (@geneids[0..$#geneids]){
        my $genestart = $geneposHash{$enhchrom}{$geneid}{"start"};
        my $geneend = $geneposHash{$enhchrom}{$geneid}{"end"};
        my $min = min(abs($geneend-$enhstart),abs($genestart-$enhend));
        if($min < 5000000){
            print OUT "$enhid\t$geneid\t$min\n";
        }
    }
}
close OUT;

# perl step2.pl --enhfile /media/yuhua/yuhua_projects/enhProj/ENHData/enhancer_annotation_step6_GSR_6column.bed --genefile /media/yuhua/yuhua_projects/genomeanno/gencode.vM17.annotation.gff3 --outfile /media/yuhua/yuhua_projects/enhProj/ENHData/enh_gene_within_5mb_GSR.txt
# perl step2.pl --enhfile /media/yuhua/yuhua_projects/enhProj/ENHData/enhancer_annotation_step6_XW_6column.bed --genefile /media/yuhua/yuhua_projects/genomeanno/gencode.vM17.annotation.gff3 --outfile /media/yuhua/yuhua_projects/enhProj/ENHData/enh_gene_within_5mb_XW.txt
