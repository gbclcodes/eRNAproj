#!/usr/bin/perl -w
# split -l 2000 enhancer_annotation_step2.bed
# mv x* splitfiles/
use strict;
use Getopt::Long;

my ($indir,$help);

GetOptions(
	"indir=s" => \$indir,
	"help!" => \$help,
);

my @infiles = `find $indir -name "x*"`;
foreach my $infile (@infiles){
	chomp $infile;
	my $taskNum = `ps -aux | grep identifyEnhancer | wc -l`;
	while($taskNum > 10){
		print "The num of task remaining $taskNum\n";
		$taskNum = `ps -aux | grep identifyEnhancer | wc -l`;
	}
	print "$infile\n";
	my $out = system("Rscript --vanilla /public/ZhangJinLab/project_enhancer/identifyEnhancerTwoDirectionTranscription.R $infile &");
	if($out==0){
		print "The task of $infile is successfully submitted\n";
	}
}

# perl step3_identify_two_direction_transcription_main.pl --indir /public/ZhangJinLab/project_enhancer/splitfiles

