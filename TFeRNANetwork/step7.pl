#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use File::Path;

my ($motifdir,$outfile,$help);

GetOptions(
	"motifdir=s" => \$motifdir,
	"outfile=s" => \$outfile,
	"help!" => \$help,
);

my @motifseqfiles = `find $motifdir -name "sequences.tsv"`;

open(OUT,">$outfile") or die "$!\n";
foreach my $motifseqfile (@motifseqfiles){
	chomp $motifseqfile;
	my $enfile = $motifseqfile;
	$enfile =~ s/sequences/ame/;
	my $stat ="No";
	open(EN,"<$enfile") or die "$!\n";
	while(<EN>){
		my $line = $_;
		chomp $line;
		next if($line =~ /^rank/);
		my @fieldValues = split /\t/,$line;
		if(defined $fieldValues[6] && $fieldValues[6] < 0.01){
			$stat="Yes";
		}
		last;
	}
	close EN;
	$motifseqfile=~/mememotiffile\/(.*?)\/.*/;
	my $stageName = $1;
	open(MS,"<$motifseqfile") or die "$!\n";
	while(<MS>){
		my $line = $_;
		chomp $line;
		my @fieldValues = split /\s+/,$line;
		if(defined $fieldValues[6] && $fieldValues[6] eq "tp"){
			if($fieldValues[3] !~ /shuf/){
				print OUT "$fieldValues[1]\t$fieldValues[3]\t$fieldValues[5]\t$stageName\t$stat\n";
			}
		}
	}
	close MS;
}
close OUT;

# perl step7.pl --motifdir /media/yuhua/yuhua_projects/enhProj/ENHData/enh_stage_group_files/GSR/mememotiffile --outfile /media/yuhua/yuhua_projects/enhProj/ENHData/enh_stage_group_files/GSR/mememotiffile/TF_enhancer_network.txt
# perl step7.pl --motifdir /media/yuhua/yuhua_projects/enhProj/ENHData/enh_stage_group_files/XW/mememotiffile --outfile /media/yuhua/yuhua_projects/enhProj/ENHData/enh_stage_group_files/XW/mememotiffile/TF_enhancer_network.txt
