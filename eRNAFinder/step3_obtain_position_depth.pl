#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use File::Path;

my ($infilepath,$bedfile,$outfilepath,$help);

GetOptions(
	"infilepath=s" => \$infilepath,
	"bedfile=s" => \$bedfile,
	"outfilepath=s" => \$outfilepath,
	"help!" => \$help,
);

my @bamfiles = `find $infilepath -name "*bam"`;
foreach my $bamfile (@bamfiles){
	chomp $bamfile;
	$bamfile =~ /.*\/(.*).bam/;
    my $sample_id = $1;
	
	if(!-e "$outfilepath/$sample_id"){
		mkpath("$outfilepath/$sample_id",0644);
		if($@){
			print "Make path $outfilepath/$sample_id failed:\n";
			exit(1);
		}
    }
	open(SH,">$outfilepath/$sample_id/${sample_id}_depth.sh") or die "$!\n";
	print SH "samtools depth -a -b $bedfile $bamfile > $outfilepath/$sample_id/depth.txt\n";
	close SH;
  
	my $out = system("sh $outfilepath/$sample_id/depth.sh 1>>$outfilepath/$sample_id/std.log 2>>$outfilepath/$sample_id//error.log &");
	if($out==0){
		print "The task of $sample_id is successfully submitted\n";
	}
}

# perl step3_obtain_position_depth.pl --infilepath /media/yuhua/yuhua_projects/enhProj/GSRData/RNAseq_GSR --bedfile /media/yuhua/yuhua_projects/enhProj/ENHData/enhancer_annotation_step2.bed --outfilepath /media/yuhua/yuhua_projects/enhProj/ENHData/depthfiles_GSR 
# perl step3_obtain_position_depth.pl --infilepath /media/yuhua/yuhua_projects/enhProj/GSRData/RNAseq_GSR --bedfile /media/yuhua/yuhua_projects/enhProj/GENEData/known_genes.bed --outfilepath /media/yuhua/yuhua_projects/enhProj/GENEData/depthfiles_GSR 
