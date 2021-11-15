#!/usr/bin/env python3

import sys, os

def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)


class GFFFeature:

    SOFA = """
        SO:0000000 ! Sequence Ontology
    SO:0000001 ! region
    SO:0000004 ! interior coding exon
    SO:0000005 ! satellite DNA
    SO:0000006 ! PCR product
    SO:0000007 ! read pair
    SO:0000013 ! scRNA
    SO:0000038 ! match set
    SO:0000039 ! match part
    SO:0000050 ! gene part
    SO:0000057 ! operator
    SO:0000059 ! nuclease binding site
    SO:0000101 ! transposable element
    SO:0000102 ! expressed sequence match
    SO:0000103 ! clone insert end
    SO:0000104 ! polypeptide
    SO:0000109 ! sequence variant obs
    SO:0000110 ! sequence feature
    SO:0000112 ! primer
    SO:0000113 ! proviral region
    SO:0000114 ! methylated C
    SO:0000120 ! protein coding primary transcript
    SO:0000139 ! ribosome entry site
    SO:0000140 ! attenuator
    SO:0000141 ! terminator
    SO:0000143 ! assembly component
    SO:0000147 ! exon
    SO:0000148 ! supercontig
    SO:0000149 ! contig
    SO:0000150 ! read
    SO:0000151 ! clone
    SO:0000159 ! deletion
    SO:0000161 ! methylated A
    SO:0000162 ! splice site
    SO:0000163 ! five prime cis splice site
    SO:0000163 ! five prime splice site
    SO:0000164 ! three prime cis splice site
    SO:0000164 ! three prime splice site
    SO:0000165 ! enhancer
    SO:0000167 ! promoter
    SO:0000177 ! cross genome match
    SO:0000178 ! operon
    SO:0000179 ! clone insert start
    SO:0000181 ! translated nucleotide match
    SO:0000183 ! non transcribed region
    SO:0000185 ! primary transcript
    SO:0000187 ! repeat family
    SO:0000188 ! intron
    SO:0000193 ! RFLP fragment
    SO:0000195 ! coding exon
    SO:0000196 ! five prime coding exon coding region
    SO:0000196 ! five prime exon coding region
    SO:0000197 ! three prime coding exon coding region
    SO:0000197 ! three prime exon coding region
    SO:0000198 ! noncoding exon
    SO:0000200 ! five prime coding exon
    SO:0000203 ! UTR
    SO:0000204 ! five prime UTR
    SO:0000205 ! three prime UTR
    SO:0000209 ! rRNA primary transcript
    SO:0000233 ! mature transcript
    SO:0000234 ! mRNA
    SO:0000235 ! TF binding site
    SO:0000236 ! ORF
    SO:0000239 ! flanking region
    SO:0000252 ! rRNA
    SO:0000253 ! tRNA
    SO:0000274 ! snRNA
    SO:0000275 ! snoRNA
    SO:0000276 ! miRNA
    SO:0000289 ! microsatellite
    SO:0000294 ! inverted repeat
    SO:0000296 ! origin of replication
    SO:0000303 ! clip
    SO:0000305 ! modified base
    SO:0000305 ! modified base site
    SO:0000306 ! methylated base feature
    SO:0000307 ! CpG island
    SO:0000314 ! direct repeat
    SO:0000315 ! TSS
    SO:0000316 ! CDS
    SO:0000318 ! start codon
    SO:0000319 ! stop codon
    SO:0000324 ! tag
    SO:0000325 ! rRNA large subunit primary transcript
    SO:0000326 ! SAGE tag
    SO:0000330 ! conserved region
    SO:0000331 ! STS
    SO:0000332 ! coding conserved region
    SO:0000333 ! exon junction
    SO:0000334 ! nc conserved region
    SO:0000336 ! pseudogene
    SO:0000337 ! RNAi reagent
    SO:0000340 ! chromosome
    SO:0000341 ! chromosome band
    SO:0000343 ! match
    SO:0000344 ! splice enhancer
    SO:0000345 ! EST
    SO:0000347 ! nucleotide match
    SO:0000349 ! protein match
    SO:0000353 ! sequence assembly
    SO:0000360 ! codon
    SO:0000366 ! insertion site
    SO:0000368 ! transposable element insertion site
    SO:0000370 ! small regulatory ncRNA
    SO:0000372 ! enzymatic RNA
    SO:0000374 ! ribozyme
    SO:0000375 ! rRNA 5 8S
    SO:0000375 ! rRNA 5.8S
    SO:0000380 ! hammerhead ribozyme
    SO:0000385 ! RNase MRP RNA
    SO:0000386 ! RNase P RNA
    SO:0000390 ! telomerase RNA
    SO:0000391 ! U1 snRNA
    SO:0000392 ! U2 snRNA
    SO:0000393 ! U4 snRNA
    SO:0000394 ! U4atac snRNA
    SO:0000395 ! U5 snRNA
    SO:0000396 ! U6 snRNA
    SO:0000397 ! U6atac snRNA
    SO:0000398 ! U11 snRNA
    SO:0000399 ! U12 snRNA
    SO:0000403 ! U14 snoRNA
    SO:0000404 ! vault RNA
    SO:0000405 ! Y RNA
    SO:0000407 ! rRNA 18S
    SO:0000409 ! binding site
    SO:0000410 ! protein binding site
    SO:0000412 ! restriction fragment
    SO:0000413 ! sequence difference
    SO:0000418 ! signal peptide
    SO:0000419 ! mature protein region
    SO:0000436 ! ARS
    SO:0000441 ! ss oligo
    SO:0000442 ! ds oligo
    SO:0000454 ! rasiRNA
    SO:0000462 ! pseudogenic region
    SO:0000464 ! decayed exon
    SO:0000468 ! golden path fragment
    SO:0000472 ! tiling path
    SO:0000474 ! tiling path fragment
    SO:0000483 ! nc primary transcript
    SO:0000484 ! three prime coding exon noncoding region
    SO:0000486 ! five prime coding exon noncoding region
    SO:0000499 ! virtual sequence
    SO:0000502 ! transcribed region
    SO:0000551 ! polyA signal sequence
    SO:0000553 ! polyA site
    SO:0000577 ! centromere
    SO:0000581 ! cap
    SO:0000587 ! group I intron
    SO:0000588 ! autocatalytically spliced intron
    SO:0000590 ! SRP RNA
    SO:0000593 ! C D box snoRNA
    SO:0000602 ! guide RNA
    SO:0000603 ! group II intron
    SO:0000605 ! intergenic region
    SO:0000610 ! polyA sequence
    SO:0000611 ! branch site
    SO:0000612 ! polypyrimidine tract
    SO:0000616 ! transcription end site
    SO:0000624 ! telomere
    SO:0000625 ! silencer
    SO:0000627 ! insulator
    SO:0000628 ! chromosomal structural element
    SO:0000643 ! minisatellite
    SO:0000644 ! antisense RNA
    SO:0000645 ! antisense primary transcript
    SO:0000646 ! siRNA
    SO:0000649 ! stRNA
    SO:0000650 ! small subunit rRNA
    SO:0000651 ! large subunit rRNA
    SO:0000652 ! rRNA 5S
    SO:0000653 ! rRNA 28S
    SO:0000655 ! ncRNA
    SO:0000657 ! repeat region
    SO:0000658 ! dispersed repeat
    SO:0000662 ! spliceosomal intron
    SO:0000667 ! insertion
    SO:0000668 ! EST match
    SO:0000673 ! transcript
    SO:0000684 ! nuclease sensitive site
    SO:0000687 ! deletion junction
    SO:0000688 ! golden path
    SO:0000689 ! cDNA match
    SO:0000694 ! SNP
    SO:0000695 ! reagent
    SO:0000696 ! oligo
    SO:0000699 ! junction
    SO:0000700 ! remark
    SO:0000701 ! possible base call error
    SO:0000702 ! possible assembly error
    SO:0000703 ! experimental result region
    SO:0000704 ! gene
    SO:0000705 ! tandem repeat
    SO:0000706 ! trans splice acceptor site
    SO:0000714 ! nucleotide motif
    SO:0000715 ! RNA motif
    SO:0000717 ! reading frame
    SO:0000719 ! ultracontig
    """
    def __init__(self, line, identifier):
        fields = line.rstrip().split('\t')
        self._identifier = identifier
        self.chrom = fields[0]
        self.source = fields[1]
        self.type = fields[2]
        self.start = fields[3]
        self.end = fields[4]
        self.score = fields[5]
        self.strand = fields[6]
        self.frame = fields[7]
        self.attributes_line = fields[8]
        self.attributes = []
        self.validSOFA = self.parseSOFA(GFFFeature.SOFA)
        self.valid = True
        self.validate()
        self.parseAttributes()


    def parseSOFA(self,SOFA):
        """
        Parse SOFA into a dictionary
        """
        SOFA_dict = {}
        for line in SOFA.split('\n'):
            fields = line.split('!')
            if len(fields) > 1:
                SOFA_dict[fields[1].strip()] = fields[0].rstrip()

        return SOFA_dict

    def validate(self):

        if self.type not in self.validSOFA and not self.type.startswith('SO:'):
            eprint("WARNING: Unrecognized feature type: {}".format(self.type))
            self.valid = False
        else:
            self.valid = True
        # Check that start, end are integers
        try:
            self.start = int(self.start)
            self.end = int(self.end)
        except ValueError:
            eprint(  "ERROR: Start and end coordinates must be integers")
            

        if self.start > self.end:
            eprint( "ERROR: Start coordinate is greater than end coordinate: {} > {}".format(self.start, self.end))
            exit(1)
            
        if self.strand not in ["+", "-", "."]:
            eprint( "ERROR: Unrecognized strand: {}".format(self.strand))
            exit(1)
            
        if self.frame not in ["0", "1", "2", "."]:
            eprint( "ERROR: Unrecognized frame: {}".format(self.frame))
            exit(1)
   

    def parseAttributes(self):
        self.attributes = self.attributes_line.rstrip(';').split(';')
        for attribute in self.attributes:
            key, value = attribute.split('=')
            if key == self._identifier:
                self.id = value
            setattr(self, key, value)
        if not hasattr(self, self._identifier):
            eprint( "WARNING: {} attribute not found in {}".format(self._identifier, self.attributes_line))
            self.id = None
            self.valid = False
            

        

    def __str__(self):
        return '\t'.join([self.chrom, self.source, self.type, str(self.start), str(self.end), self.score, self.strand, self.frame, ";".join(self.attributes)])

def gffRecord(filename, identifier="ID"):
    """
    Iterator that returns the next GFF record from a file
    """
    header, fasta = False, False
    valid_headers = 0
    with open(filename) as f:
        for line in f:
            if fasta:
                # File is finished (sequence area)
                continue
            if line.startswith('#'):
                if line.startswith('##gff-version') or line.startswith('##sequence-region'):
                    valid_headers += 1
                if header:
                    fasta = True
                continue
            line = line.rstrip()
            if not line:
                continue
            fields = line.split('\t')
            if len(fields) == 9:
                if valid_headers >0:
                    header = True
                else:
                    eprint("ERROR: GFF Header line not found")
                    sys.exit(1)
                try:
                    r = GFFFeature(line, identifier)
                    yield r
                except Exception as e:
                    eprint("ERROR: {}".format(e))
                    sys.exit(1)
            else:
                eprint("ERROR: Invalid fields count {} in {}".format(len(fields), line))
                sys.exit(1)
            


def read_fasta(path):
	import gzip
	name = None
	with (gzip.open if path.endswith('.gz') else open)(path, 'rt') as fasta:
		for line in fasta:
			if line.startswith('>'):
				if name is not None:
					yield name, seq
				name = line[1:].rstrip()
				seq = ''
			else:
				seq += line.rstrip()
	yield name, seq

if __name__ == "__main__": 
    import argparse
    parser = argparse.ArgumentParser(description="Convert GFF to BED format.")
    parser.add_argument("GFF", help="Input GFF file.")
    parser.add_argument("--id", help="Identifier [default: %(default)s]", default="ID")
    parser.add_argument("-t", "--type", help="Feature type [default: %(default)s]", default="CDS")

    args = parser.parse_args()
    
    eprint("# Input GFF file:", args.GFF)

    tot_noid = 0
    for record in gffRecord(args.GFF, args.id):
        bed = [record.chrom, str(record.start-1), str(record.end), record.id]
        if record.id is None:
            tot_noid += 1
            eprint("WARNING: Skipping record without ID")
            continue
        if record.type != args.type:
            continue
        print("\t".join(bed))