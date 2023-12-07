#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use File::Path;

my ($inputdir,$outputdir,$motifdir,$help);
GetOptions(
	"inputdir=s" => \$inputdir,
    "outputdir=s" => \$outputdir,
	"motifdir=s" => \$motifdir,
	"help!" => \$help,
);

my @fastafiles = `find $inputdir -name "*fasta"`;
foreach my $fastafile (@fastafiles){
    chomp $fastafile;
    $fastafile =~ /.*\/(.*)\_enh.fasta/;
    my $stagename = $1;
    my @motiffiles = `find $motifdir -name "*meme"`;
    foreach my $motiffile (@motiffiles){
		$motiffile =~ /.*\/(.*)\.meme/;
		my $motifname = $1;
		if(!-e "$outputdir/mememotiffile/$stagename/$motifname"){
			mkpath("$outputdir/mememotiffile/$stagename/$motifname",0644);
			if($@){
				print "Make path $outputdir/mememotiffile/$stagename/$motifname failed:$@\n";
				exit(1);
			}

			my $out = system("ame --control --shuffle-- --oc $outputdir/mememotiffile/$stagename/$motifname $fastafile $motiffile");
			if($out == 0){
				print "The task of $stagename/$motifname is successfully submitted\n";
			}
		}
	}
}

# perl step3.pl --inputdir /media/yuhua/yuhua_projects/enhProj/ENHData/enh_stage_group_files/GSR --outputdir /media/yuhua/yuhua_projects/enhProj/ENHData/enh_stage_group_files/GSR --motifdir /media/yuhua/yuhua_projects/enhProj/annodata/JASPAR_CORE
# perl step3.pl --inputdir /media/yuhua/yuhua_projects/enhProj/ENHData/enh_stage_group_files/XW --outputdir /media/yuhua/yuhua_projects/enhProj/ENHData/enh_stage_group_files/XW --motifdir /media/yuhua/yuhua_projects/enhProj/annodata/JASPAR_CORE
