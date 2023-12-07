awk -v OFS='\t' '{print $1,"EnhDataset","enh",$2,$3,".",".",".","enh_id \""$5"\"; enh_source \""$4"\""}' merge04.bed > SuperE3.merge.gtf
cat GSR.list|while read i
do
htseq-count -f bam -r pos -s no -a 10 -t enh -i enh_id -m union --nonunique=none GSR_bam/$i.sorted.bam SuperE3.merge.gtf > htseq_count/GSR.$i.count
done

cat XW.list|while read i
do
htseq-count -f bam -r pos -s no -a 10 -t enh -i enh_id -m union --nonunique=none XW_bam/$i.sorted.bam SuperE3.merge.gtf > htseq_count/XW.$i.count
done

#generate count csv
cd htseq_count
cat ../GSR.list |while read i; do sed -i "1 i Enh\t$i" GSR.$i.count; awk '{print $2}' GSR.$i.count > GSR.$i.count.new; done
cat ../XW.list |while read i; do sed -i "1 i Enh\t$i" XW.$i.count;awk '{print $2}' XW.$i.count > XW.$i.count.new; done
cat ../SuperE3.merge.gtf|cut -f1,4,5|sed "1 i Enh\tstart\tend" > location
paste GSR.MIIOocyte.count GSR.2cell.count.new GSR.4cell.count.new GSR.8cell.count.new  GSR.morula.count.new GSR.ICM.count.new GSR.TE.count.new GSR.EpiE65.count.new GSR.ExeE65.count.new location > GSR.count
awk -v FS='\t' -v OFS=',' '{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12}' GSR.count> GSR.count.csv
paste XW.MIIOocyte.count XW.Zygote.count.new XW.E2C.count.new XW.L2C.count.new XW.M4C.count.new XW.M8C.count.new XW.ICM.count.new location > XW.count
awk -v FS='\t' -v OFS=',' '{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10}' XW.count> XW.count.csv