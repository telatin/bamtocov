#!/usr/bin/env python
"""
Count reads in multiple BAM files using a target file (BED or GFF or GTF)
"""
import os, sys
import subprocess
import concurrent.futures
class Options(object):
    """
    Class to store options
    """
    def __init__(self, args):
        self.binary = args.binary
        self.feature = args.feat
        self.id = args.id
        self.verbose = args.verbose
        self.stranded = args.stranded

def eprint(*args, **kwargs):
    """
    Print to stderr
    """
    print(*args, file=sys.stderr, **kwargs)

def counts(bamfile, targetfile, options):
    """
    Count reads in a BAM file using a target file (BED or GFF or GTF)
    """
    stranded = "--stranded" if options.stranded else ""
    command = [options.binary, targetfile, stranded, bamfile]
    counts = {}
    # get command output and stderr
    p = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = p.communicate()
    # check for errors
    if p.returncode != 0:
        eprint("Error:", stderr.decode("utf-8"))
        sys.exit(1)
    # parse output
    for line in stdout.decode("utf-8").split("\n"):
        if len(line) == 0:
            continue
        if line[0] == '#':
            continue
        fields = line.split()
        if len(fields) < 1:
            continue

        # KEY: "Chromosome\tStart\tEnd\tName"
        f = fields[0]

        if options.stranded:
            try:
                fwd = int(fields[1])
                rev = int(fields[2])
                counts[f] = f"{fwd},{rev}"
            except ValueError:
                print("Error:", fields[1])

        else:
            try:
                c = int(fields[1])
                counts[f] = c
            except ValueError:
                print("Error:", fields[1])
    
    # Return FILENAME, COUNTS[feature, int]
  
    return bamfile, counts

def stranded_header(list):
    output = []
    for i in list:
        output.append(f"{i}_Fwd")
        output.append(f"{i}_Rev")
    return output

def stranded_counts(counts):
    output = []
    for i in counts:
        fwd, rev = i.split(",")
        output.append(fwd)
        output.append(rev)
    return output
    
if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Count reads in multiple BAM files using a target file (BED or GFF or GTF)")
    parser.add_argument("BAM", help="BAM file", nargs="+")
    parser.add_argument("-t", "--target", help="Target file", required=True)
    parser.add_argument("-s", "--stranded", help="Stranded", action="store_true")
    
    parser.add_argument("-o", "--output", help="Output file", required=False)
    parser.add_argument("--feat", help="Feature type [default: %(default)s]", default="exon")
    parser.add_argument("--id", help="ID attribute [default: %(default)s]", default="gene_id")
    parser.add_argument("--binary", help="Binary to bamtocounts [default: %(default)s]", default="bamtocounts")
    parser.add_argument("--verbose", help="Verbose output", action="store_true")
    opts = parser.parse_args()

    options = Options(opts)
    if opts.output is None:
        out = sys.stdout
    else:
        out = open(opts.output, 'w')
    
    matrix = {}
    first = ""

    with concurrent.futures.ThreadPoolExecutor() as executor:
        count_results = [executor.submit(counts, bam, opts.target, options) for bam in opts.BAM]
        for f in concurrent.futures.as_completed(count_results):
            sample, c = f.result()
            
            sample = "".join(os.path.basename(sample).split(".")[:-1])
            first = sample if first == "" else first
            matrix[sample] = c
            if opts.verbose:
                eprint(f"Received results: {sample}")

 
    # Print all samples names (matrix keys)
    if opts.stranded:
        # Print all samples twice (for and rev)
        print("Feature\t", "\t".join(stranded_header(matrix.keys())), file=out)
    else:
        print("Feature\t", "\t".join(matrix.keys()), file=out)
     
    # Print all features and counts
    for feature in matrix[first]:
         
        counts = [matrix[sample][feature] for sample in matrix]
        if opts.stranded:
            
            print(feature, "\t", "\t".join(map(str, stranded_counts(counts))), file=out)
        else:
            print(feature, "\t", "\t".join(map(str, counts)), file=out)
       