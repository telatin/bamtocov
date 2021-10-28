#!/usr/bin/env python3

import os, sys
import subprocess
import random
def getOtherPairFlag(flag):
    """
    Convert SAM flag to other pair flag:
    if read paired (0x1), return None

    If first in pair (0x40), return 0x80,
    if  read reverse strand (0x10), mate reverse strand (0x20)
    """
    # Check if bit 0x1 is set
    if not flag & 0x1:
        return flag
    else:
        # If 0x40, switch to 0x80
        if flag & 0x40 and flag & 0x10:
            return flag ^ 0x20 ^ 0x10 ^ 0x40 ^ 0x80

def intToBinary(i):
    """
    Convert integer to binary string, reversed
    """
    return "{0:b}".format(i)[::-1]
 
def removeCharsFromList(list, ch):
    """
    remove all elements from list that are equal to ch
    """
    return [x for x in list if x != ch]

def binaryToInt(string):
    """
    Convert binary string to integer
    """
    return int(string[::-1], 2)

 
if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Make flag file for bamtocov.")
    parser.add_argument("-f", "--flags", required=True, help="Input file with list of flags.")
    parser.add_argument("-o", "--output", required=True, help="Output prefix")
    parser.add_argument("-p", "--paired", action="store_true", help="Make paired files.")
    
    args = parser.parse_args()

    with open(args.flags, "r") as f:
        flags = [int(x) for x in f.read().split()]
        for flag in flags:
            samFile = args.output + "_" + str(flag) + ".tmp.sam"
            bamFile = args.output + "_" + str(flag) + ".bam"
            pos = 1 #random.randint(0, 800)
            read = ["r." + str(flag), str(flag), "seq", str(pos), "60", "1" + str(flag) + "M", "*", "0", "0", "*", "*"] 
            sam = "@SQ\tSN:seq\tLN:500000\n"
            sam += "\t".join(read) + "\n"
            with open(samFile, "w") as out:
                out.write(sam)

            cmd = ["samtools", "view", "-bS", samFile, "-o", bamFile]
            subprocess.call(cmd)

            # Remove sam file
            os.remove(samFile)


            
            
    