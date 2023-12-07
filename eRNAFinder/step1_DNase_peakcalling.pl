#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use File::Path;

my ($inputdir,$outputdir,$mapfile,$help);

GetOptions(
	"inputdir|i=s" => \$inputdir,
	"outputdir|o=s" => \$outputdir,
    "mapfile=s" => \$mapfile,
	"help!" => \$help,
);

my @samples = `find $inputdir -name "accepted_hits_NHi1.sorted.unique.sam"`;
foreach my $sample (@samples){
    chomp $sample;
    $sample =~ /.*\/(.*)\/accepted_hits_NHi1.sorted.unique.sam/;
    my $sample_id = $1;
    if(!-e "$outputdir/$sample_id"){
        mkpath("$outputdir/$sample_id",0644);
        if($@){
            print "Make path $outputdir/$sample_id failed:\n";
            exit(1);
        }
    }
    open(SH,">$outputdir/$sample_id/${sample_id}_macs2.sh") or die "$!\n";
    print SH "macs2 callpeak -t $sample -g mm -n $sample_id -B --nomodel -q 0.05 --shift -100 --extsize 200 --outdir $outputdir/$sample_id\n";
    close SH;
  
    my $out = system("sh $outputdir/$sample_id/${sample_id}_macs2.sh 1>>$outputdir/$sample_id/std.log 2>>$outputdir/$sample_id//error.log &");
    if($out==0){
        print "The task of $sample_id is successfully submitted\n";
    }
}

# perl step2_DNase_peakcalling.pl --inputdir /public/ZhangJinLab/project_enhancer/DNasefile/bowtie2file --outputdir /public/ZhangJinLab/project_enhancer/DNasefile/macs2file
