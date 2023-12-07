#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use File::Path;

my ($inputdir,$outputdir,$indexdir,$bowtie2dir,$picarddir,$threads,$help);

GetOptions(
	"inputdir|i=s" => \$inputdir,
	"outputdir|o=s" => \$outputdir,
    "indexdir=s" => \$indexdir,
	"bowtie2dir=s" => \$bowtie2dir,
	"picarddir=s" => \$picarddir,
	"threads=s" => \$threads,
	"help!" => \$help,
);

my @samples = `find $inputdir -name "*_1.fastq"`;
print join("\n",@samples)."\n";
foreach my $sample_p1 (@samples){
    if($sample_p1 !~ /fastq/){
        chomp $sample_p1;
        $sample_p1 =~ /.*\/(.*)\_1.fastq/;
        my $sample_id = $1;
        if(!-e "$outputdir/$bowtie2dir/$sample_id"){
            mkpath("$outputdir/$bowtie2dir/$sample_id",0644);
            if($@){
                print "Make path $outputdir/$bowtie2dir/$sample_id failed:\n$@";
                exit(1);
            }
        }
        
        open(SH,">$outputdir/$bowtie2dir/$sample_id/${sample_id}_expcal.sh") or die "$!\n";
        if(!-e "$outputdir/$bowtie2dir/$sample_id/accepted_hits.sam" || -z "$outputdir/$bowtie2dir/$sample_id/accepted_hits.sam"){
			if(-e "$inputdir/$sample_id/$sample_id\_2.fastq"){
                print SH "bowtie2 -x $indexdir -p $threads -t -q -N 1 -L 25 -X 2000 --no-mixed --no-discordant --rg-id $sample_id --rg SM:$sample_id -1 $inputdir/$sample_id/$sample_id\_1.fastq -2 $inputdir/$sample_id/$sample_id\_2.fastq -S $outputdir/$bowtie2dir/$sample_id/accepted_hits.sam\n";
                print SH "samtools view -bS $outputdir/$bowtie2dir/$sample_id/accepted_hits.sam -o $outputdir/$bowtie2dir/$sample_id/accepted_hits.bam\n";
                print SH "samtools sort $outputdir/$bowtie2dir/$sample_id/accepted_hits.bam -o $outputdir/$bowtie2dir/$sample_id/accepted_hits.sorted.bam\n";
                print SH "samtools index $outputdir/$bowtie2dir/$sample_id/accepted_hits.sorted.bam\n";
            }else{
                print SH "bowtie2 -x $indexdir -p $threads -t -q -N 1 -L 25 --no-mixed --no-discordant --rg-id $sample_id --rg SM:$sample_id -U $inputdir/$sample_id/$sample_id\.sra_1.fastq -S $outputdir/$bowtie2dir/$sample_id/accepted_hits.sam\n";
                print SH "samtools view -bS $outputdir/$bowtie2dir/$sample_id/accepted_hits.sam -o $outputdir/$bowtie2dir/$sample_id/accepted_hits.bam\n";
                print SH "samtools sort $outputdir/$bowtie2dir/$sample_id/accepted_hits.bam -o $outputdir/$bowtie2dir/$sample_id/accepted_hits.sorted.bam\n";
                print SH "samtools index $outputdir/$bowtie2dir/$sample_id/accepted_hits.sorted.bam\n";
            }
        }
        
        if(!-e "$outputdir/$bowtie2dir/$sample_id/accepted_hits_NHi1.sam" || -z "$outputdir/$bowtie2dir/$sample_id/accepted_hits_NHi1.sam"){
            print SH "grep -v -E -w 'NH:i:2|NH:i:3|NH:i:4|NH:i:5|NH:i:6|NH:i:7|NH:i:8|NH:i:9|NH:i:10|NH:i:11|NH:i:12|NH:i:13|NH:i:14|NH:i:15|NH:i:16|NH:i:17|NH:i:18|NH:i:19|NH:i:20' $outputdir/$bowtie2dir/$sample_id/accepted_hits.sam > $outputdir/$bowtie2dir/$sample_id/accepted_hits_NHi1.sam\n";
        }
        
        if(!-e "$outputdir/$bowtie2dir/$sample_id/accepted_hits_NHi1.sorted.sam" || -z "$outputdir/$bowtie2dir/$sample_id/accepted_hits_NHi1.sorted.sam"){
            print SH "samtools sort -@ $threads -o $outputdir/$bowtie2dir/$sample_id/accepted_hits_NHi1.sorted.sam $outputdir/$bowtie2dir/$sample_id/accepted_hits_NHi1.sam\n";
        }
        
        if(!-e "$outputdir/$bowtie2dir/$sample_id/accepted_hits_NHi1.sorted.unique.sam" || -z "$outputdir/$bowtie2dir/$sample_id/accepted_hits_NHi1.sorted.unique.sam"){
            print SH "java -Xmx15g -jar $picarddir/picard.jar MarkDuplicates I=$outputdir/$bowtie2dir/$sample_id/accepted_hits_NHi1.sorted.sam O=$outputdir/$bowtie2dir/$sample_id/accepted_hits_NHi1.sorted.unique.sam METRICS_FILE=$outputdir/$bowtie2dir/$sample_id/${sample_id}.metricsFile VALIDATION_STRINGENCY=LENIENT REMOVE_DUPLICATES=true ASSUME_SORT_ORDER=coordinate\n";
            print SH "samtools view -bS $outputdir/$bowtie2dir/$sample_id/accepted_hits_NHi1.sorted.unique.sam -o $outputdir/$bowtie2dir/$sample_id/accepted_hits_NHi1.sorted.unique.bam\n";
            print SH "samtools index $outputdir/$bowtie2dir/$sample_id/accepted_hits_NHi1.sorted.unique.bam\n";
            print SH "rm $outputdir/$bowtie2dir/$sample_id/accepted_hits.sam\n";
            print SH "rm $outputdir/$bowtie2dir/$sample_id/accepted_hits_NHi1.sam\n";
            print SH "rm $outputdir/$bowtie2dir/$sample_id/accepted_hits_NHi1.sorted.sam\n";
            print SH "rm $outputdir/$bowtie2dir/$sample_id/accepted_hits_NHi1.sorted.unique.sam\n";
        }
        
        my $taskNum =`ps -aux | grep sh | grep bowtie2 | wc -l`; 
        while($taskNum > 20){
            print "The num of task remaining $taskNum\n";
            sleep 30;
            print `date`;
            $taskNum = `ps -aux | grep sh | grep bowtie2 | wc -l`;
        }
        
        my $out = system("nohup sh $outputdir/$bowtie2dir/$sample_id/${sample_id}_expcal.sh 1>>$outputdir/$bowtie2dir/$sample_id/std.log 2>>$outputdir/$bowtie2dir/$sample_id/error.log &");
        if($out==0){
            print "The task of $sample_id is successfully submitted\n";
        }
    }
}

# perl step1_H3K27ac_aligment.pl --inputdir /public/ZhangJinLab/project_enhancer/H3K27acData --outputdir /public/ZhangJinLab/project_enhancer/H3K27acData --indexdir /media/data/TFChIPseqData/bowtie2index/GRCm38 --picarddir /software/picard-2.18.2 --bowtie2dir bowtie2file --threads 2
