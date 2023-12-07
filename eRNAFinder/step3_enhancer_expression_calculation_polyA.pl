#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use File::Path;

my ($inputdir,$outputdir,$htseqdir,$gfffile,$threads,$help);

GetOptions(
	"inputdir|i=s" => \$inputdir,
	"outputdir|o=s" => \$outputdir,
	"htseqdir=s" => \$htseqdir,
	"gfffile|gff=s" => \$gfffile,
	"threads=s" => \$threads,
	"help!" => \$help,
);

my @samples = `find $inputdir -name "*.bam"`;
print join("\n",@samples)."\n";
foreach my $sample (@samples){
	chomp $sample;
	$sample =~ /.*\/(.*).bam/;
	my $sample_id = $1;
    
	if(!-e "$outputdir/$htseqdir/$sample_id"){
		mkpath("$outputdir/$htseqdir/$sample_id",0644);
		if($@){
			print "Make path $outputdir/$htseqdir/$sample_id failed:\n";
			exit(1);
		}
	}
    
	open(SH,">$outputdir/$htseqdir/$sample_id/${sample_id}_expcal.sh") or die "$!\n";
    if(!-e "$outputdir/$htseqdir/$sample_id/htseq_count.txt" || -z "$outputdir/$htseqdir/$sample_id/htseq_count.txt"){
		print SH "htseq-count -f bam -s no -i ID -t gene --nonunique all -q $sample $gfffile > $outputdir/$htseqdir/$sample_id/htseq_count.txt\n";
	}
	close SH;
    
    my $taskNum =`ps -aux | grep expcal.sh | wc -l`; 
    while($taskNum > 5){
        print "The num of task remaining $taskNum\n";
        sleep 30;
        print `date`;
        $taskNum = `ps -aux | grep expcal.sh | wc -l`;
    }
    
	my $out = system("sh $outputdir/$htseqdir/$sample_id/${sample_id}_expcal.sh 1>>$outputdir/$htseqdir/$sample_id/std.log 2>>$outputdir/$htseqdir/$sample_id/error.log &");
    
	if($out==0){
		print "The task of $sample_id is successfully submitted\n";
	}
}

# perl enhancer_expression_calculation_polyA.pl --inputdir /media/yuhua/yuhua_projects/enhProj/XWData/RNAseqData_XW --outputdir /media/yuhua/yuhua_projects/enhProj/ENHData --htseqdir htseqcountfile --gfffile /media/yuhua/yuhua_projects/enhProj/ENHData/enhancer_annotation_step4.gff3 --threads 2 
