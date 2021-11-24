#!/usr/bin/env python3

import subprocess
import sys
import argparse

def checkRscript(rscript):
    try:
        subprocess.check_call([rscript, "--version"], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    except subprocess.CalledProcessError:
        print("ERROR: Rscript not found")
        sys.exit(1)

def rCmd(commands, bin):
    """
    # universalnewlines=TRUE
    # bufsize = 1
    """
    endstr = "<END_LINE>"
    # append "string" to all comands in list
    commands = [command + f"; cat('{endstr}\n')\n" for command in commands]
     
    r = ["stdbuf",  '-oL', bin, "--vanilla", "--no-save", "--no-restore", "-"]

    # popen R
    p = subprocess.Popen(r, stdin=subprocess.PIPE, stdout=subprocess.PIPE,   universal_newlines=True, bufsize=1, close_fds=True)
    c = 0
    # Run
    for line in commands:
        # send the line to stdin
        p.stdin.write(line)
        c += 1
        s = " " * (len(str(c)) + 2)
        print(f"{s}>\t", line.replace("\n", "\\n").replace("\t", "\\t"))
        reader = True
        o = 0
        while reader:
            o += 1
            output = p.stdout.readline().strip()
            
            if endstr in output or o > 6:
                reader = False
                print("-------------")
                continue
            else:
                print(f"{c}:{o}>", output)
 

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Test subprocess. Requires R and Rscript installed.")
    parser.add_argument("commands", nargs="*", help="Commands to run in R")
    parser.add_argument("-r", "--rscript", help="Path to Rscript [default: %(default)s]", default="Rscript")
    parser.add_argument("-v", "--verbose", action="store_true", help="Print R output")
    opts   = parser.parse_args()

    if len(opts.commands) == 0:
        commands = ["a <- c(1,2,3)", "cat(a, '\n')", "a[1]", "cat('one\ntwo\nthree\nfour')", "cat(length(a), '\n')" ]
    else:
        commands = opts.commands

    if opts.verbose:
        print(len(commands), "commands")
    
    checkRscript(opts.rscript)
    rCmd(commands, bin=opts.rscript)



quit()

