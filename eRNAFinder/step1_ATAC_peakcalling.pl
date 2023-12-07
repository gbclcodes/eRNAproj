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

my %mapHash;
open(MAP,"<$mapfile") or die "$!\n";
while(<MAP>){
    my $line = $_;
    chomp $line;
    my @fieldValues = split /\s+/,$line;
    $mapHash{$fieldValues[0]} = 1;
}
close MAP;

my @samples = `find $inputdir -name "accepted_hits_NHi1.sorted.unique.bam"`;
foreach my $sample (@samples){
    chomp $sample;
    $sample =~ /.*\/(.*)\/accepted_hits_NHi1.sorted.unique.bam/;
    my $sample_id = $1;
    if(exists $mapHash{$sample_id}){
        if(!-e "$outputdir/$sample_id"){
            mkpath("$outputdir/$sample_id",0644);
            if($@){
                print "Make path $outputdir/$sample_id failed:\n";
                exit(1);
            }
        }
        open(SH,">$outputdir/$sample_id/${sample_id}_macs2.sh") or die "$!\n";
        print SH "macs2 callpeak -t $sample -g mm -n $sample_id -B --nomodel -q 0.05 --nolambda --outdir $outputdir/$sample_id\n";
        close SH;
      
        my $out = system("sh $outputdir/$sample_id/${sample_id}_macs2.sh 1>>$outputdir/$sample_id/std.log 2>>$outputdir/$sample_id//error.log &");
        if($out==0){
            print "The task of $sample_id is successfully submitted\n";
        }
    }
}

# perl step2_ATAC_peakcalling.pl --inputdir /public/ZhangJinLab/project_enhancer/ATACfile/bowtie2file --outputdir /public/ZhangJinLab/project_enhancer/ATACfile/macs2file --mapfile /public/ZhangJinLab/project_enhancer/ATACfile/stringtiefile/ATAC_run2samplename.txt
