#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use File::Path;

my ($genomefile,$enhposfile,$highenhjsfile,$medianenhjsfile,$lowenhjsfile,$outputdir,$help);
GetOptions(
    "genomefile=s" => \$genomefile,
	"enhposfile=s" => \$enhposfile,
	"highenhjsfile=s" => \$highenhjsfile,
	"medianenhjsfile=s" => \$medianenhjsfile,
	"lowenhjsfile=s" => \$lowenhjsfile,
    "outputdir=s" => \$outputdir,
	"help!" => \$help,
);

open(GENOME,"<$genomefile") or die "Could not open the file $genomefile:$!\n";
my (%chromosome,$chromname);
while(<GENOME>){
	my $line = $_;
	chomp $line;
	if($line =~ />(chr.*)/){
		my @lineContent = split /\s+/,$1;
		$chromname = $lineContent[0];
		print "$chromname\n";
		$chromosome{$chromname} = "";
	}else{
        $chromosome{$chromname} .= $line;
    }
}
close GENOME;

my (%stageHash,%jsHash);
open(HI,"<$highenhjsfile") or die "$!\n";
while(<HI>){
	my $line = $_;
	chomp $line;
	my @fieldValues = split /\s+/,$line;
	push @{$stageHash{$fieldValues[2]}}, $fieldValues[0];
	$jsHash{$fieldValues[0]} = $fieldValues[1];
}
close HI;

open(ME,"<$medianenhjsfile") or die "$!\n";
while(<ME>){
	my $line = $_;
	chomp $line;
	my @fieldValues = split /\s+/,$line;
	push @{$stageHash{$fieldValues[2]}}, $fieldValues[0];
	$jsHash{$fieldValues[0]} = $fieldValues[1];
}
close ME;

open(LW,"<$lowenhjsfile") or die "$!\n";
while(<LW>){
	my $line = $_;
	chomp $line;
	my @fieldValues = split /\s+/,$line;
	push @{$stageHash{$fieldValues[2]}}, $fieldValues[0];
	$jsHash{$fieldValues[0]} = $fieldValues[1];
}
close LW;

my %posHash;
open(ENH,"<$enhposfile") or die "Could not open the file $enhposfile:$!\n";
while(<ENH>){
	my $line = $_;
	chomp $line;
	my @fieldValues = split /\s+/,$line;
	$posHash{$fieldValues[0]} = "$fieldValues[1]\t$fieldValues[2]\t$fieldValues[3]";
}
close ENH;

foreach my $stagename (keys %stageHash){
	if(!-e "$outputdir/$stagename"){
		mkpath("$outputdir/$stagename",0644);
		if($@){
			print "Make path $outputdir/$stagename failed:$@\n";
			exit(1);
		}
	}
	open(FA,">$outputdir/$stagename/$stagename\_enh.fasta") or die "Could not open the file $outputdir/$stagename/$stagename.fasta:$!\n";
	foreach my $enhid (@{$stageHash{$stagename}}){
		my @fieldValues = split /\t/,$posHash{$enhid};
		if(exists $jsHash{$enhid}){
			my $geneseq = substr($chromosome{$fieldValues[0]},$fieldValues[1],($fieldValues[2]-$fieldValues[1]+1));
			print ">$enhid\t$jsHash{$enhid}\n";
			print FA ">$enhid\t$jsHash{$enhid}\n";
			my $strnum = length($geneseq);
			my $iternum = int($strnum/50);
			for(my $iter=1;$iter<=$iternum;$iter++){
				print FA substr($geneseq,50*($iter-1),50)."\n";
			}
			print FA substr($geneseq,50*$iternum,($strnum-50*$iternum)-1)."\n";
		}
	}
	close FA;
}

# perl step1.pl --genomefile /media/yuhua/yuhua_projects/enhProj/annodata/GRCm38.p6.genome.fa --enhposfile /media/yuhua/yuhua_projects/enhProj/ENHData/enh_GSR_clean_td_0.01EPM.txt --highenhjsfile /media/yuhua/yuhua_projects/enhProj/ENHData/enhancer_stage_specific_highscores_GSR.txt --medianenhjsfile /media/yuhua/yuhua_projects/enhProj/ENHData/enhancer_stage_specific_medianscores_GSR.txt --lowenhjsfile /media/yuhua/yuhua_projects/enhProj/ENHData/enhancer_stage_specific_lowscores_GSR.txt --outputdir /media/yuhua/yuhua_projects/enhProj/ENHData/enh_stage_group_files/GSR 

# perl step1.pl --genomefile /media/yuhua/yuhua_projects/enhProj/annodata/GRCm38.p6.genome.fa --enhposfile /media/yuhua/yuhua_projects/enhProj/ENHData/enh_XW_clean_td_0.01EPM.txt --highenhjsfile /media/yuhua/yuhua_projects/enhProj/ENHData/enhancer_stage_specific_highscores_XW.txt --medianenhjsfile /media/yuhua/yuhua_projects/enhProj/ENHData/enhancer_stage_specific_medianscores_XW.txt --lowenhjsfile /media/yuhua/yuhua_projects/enhProj/ENHData/enhancer_stage_specific_lowscores_XW.txt --outputdir /media/yuhua/yuhua_projects/enhProj/ENHData/enh_stage_group_files/XW 



