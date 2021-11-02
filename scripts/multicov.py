#!/usr/bin/env python
"""
Detect low coverage regions in multiple BAM files using bamtocov, megadepth or bedtools
"""

import os, sys
import argparse
import concurrent.futures
import subprocess
import tempfile

def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)


def bamtocov(bam, min, max, bin="bamtocov"):
    """
    get the STDOUT from bamtocov and return only the intervals
    between min and max 
    """
    cmd = [bin,  "--quantize", f"{min},{max+1}", bam]
    
    try:
        output = subprocess.check_output(cmd, stderr=subprocess.STDOUT)
        lines = output.decode("utf-8").split("\n") 
        # removelast line if empty
        if lines[-1] == "":
            lines.pop()
        # Make a list of intervals splitting each line on \t
        intervals = [l.split("\t") for l in lines]
        # return lines having as fourth field "{min}-{max}""
        return [l[0:3] for l in intervals if l[3] == f"{min}-{max}"]
    except subprocess.CalledProcessError as e:
        eprint("ERROR", e.output)
        raise e

def bamtocovFile(bam, min, max, bin="bamtocov"):
    tmpFile = tempfile.NamedTemporaryFile(mode="w+t", delete=False)
    cmd = [bin,  "--quantize", f"{min},{max+1}", bam]

    cmd2 = ["grep", "-w", f"{min}-{max}"]
    
    try:
        #subprocess.check_call(cmd, stdout=tmpFile)
        p1 = subprocess.Popen(cmd, stdout=subprocess.PIPE)
        p2 = subprocess.Popen(cmd2, stdin=p1.stdout, stdout=tmpFile)
        p1.stdout.close()
        p2.communicate()
        return tmpFile.name
    except subprocess.CalledProcessError as e:
        eprint("ERROR", e.output)
        raise e 


def megadepthFile(bam, min, max, bin="megadepth"):
    tmpFile = tempfile.NamedTemporaryFile(mode="w+t", delete=False)
    cmd = [bin,  "--coverage", bam]
    try:
        subprocess.check_call(cmd, stdout=tmpFile, stderr=subprocess.DEVNULL)
        return tmpFile.name
    except subprocess.CalledProcessError as e:
        eprint("ERROR", e.output)
        raise e 

def bedtoolsFile(bam, min, max, bin="bedtools"):
    tmpFile = tempfile.NamedTemporaryFile(mode="w+t", delete=False)
    cmd = [bin,  "genomecov", "-bga", "-ibam", bam]
    print(" ".join(cmd))
    try:
        subprocess.check_call(cmd, stdout=tmpFile, stderr=subprocess.DEVNULL)
        return tmpFile.name
    except subprocess.CalledProcessError as e:
        eprint("ERROR Executing bedtools:", e.output)
        raise e 
    except Exception as e:
        eprint("ERROR", e)
        raise e
        
def processBAM(bamfile, mincov, maxcov, toolpath="bamtocov"):
    """
    Call the appropriate function depending on the tool to retreive
    a set of intervals having coverage between mincov and maxcov
    """
    if toolpath.endswith("bamtocov"):
        return bamtocovFile(bamfile, mincov, maxcov, toolpath)
    elif toolpath.endswith("megadepth"):
        tmpFile = megadepthFile(bamfile, mincov, maxcov, toolpath)
        # Filter the intervals with coverage between mincov and maxcov
        intervals = [l.split("\t") for l in open(tmpFile).readlines() if int(l.split("\t")[3]) >= mincov and int(l.split("\t")[3]) <= maxcov]
        # Save to file
        newTmp = tempfile.NamedTemporaryFile(mode="w+t", delete=False)
        newTmp.write("\n".join(["\t".join(l) for l in intervals]))
        if opts.verbose:
            eprint(f"[Megadepth] {bamfile} coverage: {tmpFile} -> {newTmp.name}")
        return newTmp.name
    elif toolpath.endswith("bedtools"):
        tmpFile = bedtoolsFile(bamfile, mincov, maxcov, toolpath)
        # Filter the intervals with coverage between mincov and maxcov
        intervals = [l.split("\t") for l in open(tmpFile).readlines() if int(l.split("\t")[3]) >= mincov and int(l.split("\t")[3]) <= maxcov]
        # Save to file
        newTmp = tempfile.NamedTemporaryFile(mode="w+t", delete=False)
        newTmp.write("\n".join(["\t".join(l) for l in intervals]))
        if opts.verbose:
            eprint(f"[Bedtools] {bamfile} coverage: {tmpFile} -> {newTmp.name}")
        return newTmp.name
    else:
        raise Exception(f"Tool {toolpath} not supported")

def intersectIntervals(intervalsLists, sharedratio=1.0):
    """
    Receve a list of intervals in the format [chr, start, end]
    return a list of intervals of regions intersected across all samples
    """
    cmd = ["bedtools", "multiinter", "-i"]
    cmd.extend(intervalsLists)
    minsamples = int(len(intervalsLists) * sharedratio)
    try:
        output = subprocess.check_output(cmd, stderr=subprocess.STDOUT)
        for file in intervalsLists:
            if not opts.verbose:
                os.remove(file)
        
        """
        seq1    0       9       2       1,2     1       1
        seq1    109     200     2       1,2     1       1
        seq1    300     650     2       1,2     1       1
        """
        
        lines = output.decode("utf-8").strip().split("\n")
        #return lines having as fourth column i>=minsamples
        return [l.split("\t") for l in lines if int(l.split("\t")[3]) >= minsamples]
    except subprocess.CalledProcessError as e:
        eprint("ERROR", e.output)
        raise e

if __name__ == "__main__":
    args = argparse.ArgumentParser(description="Detect low coverage regions in multiple BAM files using bamtocov, megadepth or bedtools")
    # Coverage
    cov = args.add_argument_group("Coverage filters")
    cov.add_argument("--min", type=int, default=0, help="Minimum coverage to consider a region as low coverage [default: 10]")
    cov.add_argument("--max", type=int, default=10, help="Maximum coverage to consider a region as low coverage [default: 10]")
    cov.add_argument("--shared", type=float, default=1.0, help="Ratio of input samples sharing the interval [default: 1.0]")

    
    # Add separate group of arguments
    tools = args.add_argument_group("External tools")
    tools.add_argument("-t", "--tool", help="Tool to use for coverage detection", choices=["bamtocov", "bedtools", "megadepth"], default="bamtocov")
    tools.add_argument("--bamtocov", help="Path to bamtocov executable", default="bamtocov")
    tools.add_argument("--bedtools", help="Path to bedtools executable", default="bedtools")
    tools.add_argument("--megadepth", help="Path to megadepth executable", default="megadepth")

    
    args.add_argument("BAM", help="BAM file(s) to process", nargs="+")
    args.add_argument("--verbose", help="Print verbose output", action="store_true")
    opts = args.parse_args()

    # Coverage analysis tool
    binpath = ""
    # Temporary files output
    tmpFiles = []

    if opts.tool == "bamtocov":
        binpath = opts.bamtocov
    elif opts.tool == "bedtools":
        binpath = opts.bedtools
    elif opts.tool == "megadepth":
        binpath = opts.megadepth
    else:
        raise Exception(f"Invalid argument: Tool {opts.tool} not supported")
     
     
    with concurrent.futures.ThreadPoolExecutor() as executor:
        results = [executor.submit(processBAM, bam, int(opts.min),  int(opts.max), binpath) for bam in opts.BAM]
        for f in concurrent.futures.as_completed(results):
            tmpFiles.append(f.result())
            if opts.verbose:
                eprint(f"Received results: {f.result()}")

    # Print lines from intersectIntervals(tmpFiles) to stdout
    for line in intersectIntervals(tmpFiles):
        line[3] = f"{opts.min}X-{opts.max}X,{round(100*int(line[3])/len(opts.BAM),2)}%"
        print("\t".join(line[0:4]))