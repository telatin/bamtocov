#!/usr/bin/env python3
import sys
inputfile = sys.argv[1]

with open(inputfile) as f:
    totaldelta = 0
    for line in f:
        line = line.strip()
        if line:
            fields = line.split()
            id1, c1, id2, c2 = fields[0], fields[1], fields[2], fields[3]
            if id1 == id2:
                delta = int(c2) - int(c1)
                p = delta / (int(c2))  * 100 if int(c2) > 0 else 0
       
                totaldelta += delta
                print(f"{id1}\t{delta}\t{p:.2f}\t{c1},{c2}")

    print(f"Total delta: {totaldelta}")