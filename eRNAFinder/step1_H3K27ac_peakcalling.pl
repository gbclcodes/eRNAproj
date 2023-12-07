#!/usr/bin/perl -w

use strict;
use Getopt::Long;
use File::Path;

my ($inputdir,$outputdir,$help);

GetOptions(
	"inputdir|i=s" => \$inputdir,
	"outputdir|o=s" => \$outputdir,
	"help!" => \$help,
);

my @samples = `find $inputdir -name "*.sorted.unique.bam"`;
foreach my $sample (@samples){
    chomp $sample;
    if($sample !~ /accepted/){
    	$sample =~ /.*\/(.*)_NHi1.sorted.unique.bam/;
    	my $sample_id = $1;
    	if(!-e "$outputdir/$sample_id"){
        	mkpath("$outputdir/$sample_id",0644);
        	if($@){
            		print "Make path $outputdir/$sample_id failed:\n";
            		exit(1);
       		}
   	}
	open(SH,">$outputdir/$sample_id/${sample_id}_macs2.sh") or die "$!\n";
	print SH "macs2 callpeak -t $sample -c $oocyte_control -g mm -n $sample_id -B --nomodel -q 0.05 --nomodel --nolambda --size=37 --extsize=73 --nolambda --outdir $outputdir/$sample_id\n";
	close SH;

	my $out = system("sh $outputdir/$sample_id/${sample_id}_macs2.sh 1>>$outputdir/$sample_id/std.log 2>>$outputdir/$sample_id//error.log &");
	if($out==0){
		print "The task of $sample_id is successfully submitted\n";
	}
}

# perl fig1a_step2_H3K27ac_peakcalling.pl --inputdir /public/ZhangJinLab/project_enhancer/H3K27acData/bowtie2file --outputdir /public/ZhangJinLab/project_enhancer/H3K27acData/macs2file
