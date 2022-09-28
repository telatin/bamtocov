#!/usr/bin/env python
"""
Read speed and memory usage for a set of benchmark file with the
filename format: mem.TOOLNAME.DatasetName.txt
"""
import sys, os, argparse
import pandas
from enum import Enum

class OutputFormat(Enum):
    csv = 'csv'
    markdown = 'md'
    dataframe = 'df'

    def __str__(self):
        return self.value

def print_df(df, title, format):
    if format == OutputFormat.csv:
        print("\n# %s" % title)
        print(df.to_csv(index=False))
    elif format == OutputFormat.markdown:
        print("\n## %s" % title)
        print(df.to_markdown())
    elif format == OutputFormat.dataframe:
        print("\n> %s" % title)
        print(df)
    else:
        raise Exception("Unknown output format: %s" % format)

def prepare_df(df, strip):
    df = df.sort_index(axis=1)
    # Strip all strings in "strip" from the column names
    for s in strip:
        df.columns = [x.replace(s, "") for x in df.columns]
    return df
    
def loadStats(filename):
    """
    Parse file in the format:
    elapsed time:   4.087
    peak rss:       190076
    """
    keys = {
        'elapsed time:': 'time',
        'peak rss:': 'mem'
    }
    with open(filename) as f:
        lines = f.readlines()
        stats = {
            "mem": -1,
            "time": -1
        }
        for line in lines:
            try:
                key, value = line.strip().split('\t')
                stats[keys[key]] = value
            except:
                continue
    return stats
if __name__ == "__main__":
    args = argparse.ArgumentParser(description="Read speed and memory usage for a set of benchmark file with the filename format: mem.TOOLNAME.DatasetName.txt")
    args.add_argument("-i", "--input", help="Input directory", required=True)
    args.add_argument("-o", "--output", help="Output file", required=False)
    args.add_argument("-f", "--format", type=OutputFormat, choices=list(OutputFormat), default=OutputFormat.csv)
    args.add_argument("--strip", help="Strip strings from dataset names", nargs="*", default=[])
    args = args.parse_args()

    # Get list of mem.*.txt files in input directory
    files = [f for f in os.listdir(args.input) if f.startswith("mem.")]
    
    memoryDf = pandas.DataFrame()
    speedDf  = pandas.DataFrame()
    for f in files:
        datasetName = f.split('.')[1]
        toolName    = "-".join(f.split('.')[2:-1])
        data = loadStats(os.path.join(args.input, f))
        print(data["mem"])
        # Populate memoryDf with columns (datasetName) and rows (toolName) and values (data["mem"])
        memoryDf.loc[datasetName, toolName] = data["mem"] 
        # Populate speedDf with columns (datasetName) and rows (toolName) and values (data["time"])
        speedDf.loc[datasetName, toolName] = data["time"]

 
    memoryDf = prepare_df(memoryDf, args.strip)
    speedDf = prepare_df(speedDf, args.strip)

    print_df(memoryDf, "Memory usage (Kb)", args.format)
    print_df(speedDf, "Speed (s)", args.format)
 