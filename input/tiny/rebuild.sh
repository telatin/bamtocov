export PATH=/Users/telatina/git/bamtocov/scripts/benchmarking/:"$PATH"
if [ ! -e target.bed ] || [ ! -e reads.bed ];then
 echo "Missing target.bed, reads.bed"
 exit
fi

for i in target*bed;
do
   echo "Converting $i to ${i/bed/gtf}"
   bed-to-gtf.py -i $i -o ${i/bed/gtf} --min-len 1000
done
echo "Making ref.fa"
bed-to-fasta.py -i target.bed --min-len 1000  > ref.fa
samtools faidx ref.fa

for i in reads*bed;
do
   echo "Making ${i/bed/bam}"
   make-bam-from-bed.py -i $i -o ${i/bed/bam} --min-len 1000
   samtools index  ${i/bed/bam} 
done

featureCounts -a target.gtf -o fc reads.bam

