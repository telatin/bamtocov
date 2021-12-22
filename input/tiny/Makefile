

all: target.gtf reads.bam.bai ref.fa.fai  fc

target.gtf: target.bed
	bed-to-gtf.py -i $< -o $@  --min-len 1000

reads.bam.bai: reads.bam
	samtools index $<

reads.bam: reads.bed
	make-bam-from-bed.py -i $< -o $@ --min-len 1000

ref.fa.fai: ref.fa
	samtools faidx $<

ref.fa: target.bed
	bed-to-fasta.py -i $< --min-len 1000  > $@

mos-target.regions.bed: target.bed reads.bam.bai
	mosdepth --by target.bed mos-target reads.bam
	gunzip mos*gz

fc: target.gtf reads.bam
	featureCounts -a target.gtf -o fc reads.bam

clean:
	rm target.gtf reads.bam* ref.fa* fc*