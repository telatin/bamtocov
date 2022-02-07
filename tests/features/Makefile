BIN=../../bin
SCRIPTS=../../scripts

all: clean compare.txt ref.fa bamtocounts_strict.tsv
speed: speed.md

ref.fa: target.bed
	$(SCRIPTS)/benchmarking/bed-to-fasta.py -i target.bed -o ref.fa

speed.md: 
	hyperfine --export-markdown speed.md "$(BIN)/bamtocounts --id gene_id --type exon target.gtf seqs.bam > bamtocounts.tsv" "featureCounts -a target.gtf -o fc seqs.bam"
compare.txt: bamtocounts.tsv fc.tsv
	head  bamtocounts.tsv fc.tsv > compare.txt
	cat compare.txt

seqs.bam: seqs.bed
	$(SCRIPTS)/benchmarking/make-bam-from-bed.py -i seqs.bed -o seqs.bam

seqs.bam.bai: seqs.bam
	samtools index seqs.bam

target.gtf: target.bed
	$(SCRIPTS)/benchmarking/bed-to-gtf.py -i target.bed -o target.gtf

bamtocounts.tsv: seqs.bam target.gtf
	$(BIN)/bamtocounts --id gene_id --type exon target.gtf seqs.bam > bamtocounts.tsv

bamtocounts_strict.tsv: seqs.bam target.gtf
	$(BIN)/bamtocounts --id gene_id --type exon --strict target.gtf seqs.bam > bamtocounts_strict.tsv

bamtocounts_legacy.tsv: seqs.bam seqs.bam.bai target.gtf
	$(BIN)/bamtocounts_legacy --id gene_id --type exon target.gtf seqs.bam > bamtocounts_legacy.tsv

fc: seqs.bam target.gtf
	featureCounts -a target.gtf -o fc seqs.bam

fc.tsv: fc
	cut -f 1,7 fc | grep -v "# Program:featureCounts" > fc.tsv
#bin/bamtocov:
#	nimble build
 

clean:
	@echo "Cleaning..."
	@for i in fc fc.summary *.tsv; \
	do \
		if [ -e "$$i" ]; then rm -f $$i; echo "Removing $$i"; else echo "$$i Not found"; fi \
	done
