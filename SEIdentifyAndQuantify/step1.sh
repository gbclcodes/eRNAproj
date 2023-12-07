#1.get all super enhancers together
cat list|while read i g
do
python2.7 ROSE_main.py -g MM10 -i gff/$g.gff -r bam/$i.sorted.bam -o rose_results/$i -s 12500 -t 2500
cat rose_results/$i/${i}_Gateway_SuperEnhancers.bed |sed "s/Enh/${i}_Enh/g" >> merge01.bed
done
#2 merge with previous enhancers
sort -k1,1 -k2,2n merge01 |sed 's/_lociStitched//g' > merge02.bed
bedtools merge -i merge02.bed -c 4 -o collapse > merge03.bed
#3 name merged super enhancers
#if the number of superenhancers is greater than 9999, the script should be change to'printf("%05d\n")'
n=1
cat merge03.bed |while read i
do
c=`echo $n|awk '{printf("%04d\n",$0)}'`
let n++
c2="SpEnh$c"
echo -e "$i\t$c2" >> merge04.bed
done
#4 get the source of merged super enhancers
cat list|while read i
do
grep "PolII_${i}" merge04.bed|cut -f 1,2,3,4,5 > tmp/$i.tmp
bedtools intersect -a tmp/$i.tmp -b gff/$i.bed -wao >> merge.s1.out
done
sort -k 5 merge.s1.out > merge.s2.out
bedtools groupby -i merge.s2.out -g 1,4,5 -c 9 -o distinct > merge.s3.out
sed -i 's/,/;/g' merge.s3.out
awk -v OFS=',' '{print $1,$2,$3,$4}' merge.s3.out > SuperE3.info.csv