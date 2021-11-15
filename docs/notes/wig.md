# Wig format

The Wiggle format is used to plot quantitaitve data in a genome browser, and a detailed description
can be found [at this website (genome.ucsc.edu)](http://genome.ucsc.edu/goldenPath/help/wiggle.html).

The native output of "BamToCov" is BED, and specifically [bedGraph](http://genome.ucsc.edu/goldenPath/help/bedgraph.html), that
can be converted to bigWig via `bedGraphToBigWig`.

For convenience we report some succint examples taken from the website aforementioned.

## Examples

### Variable step

This format is used for data with irregular intervals between new data points, and is the more commonly used wiggle format.
The `span` parameter is optional.

```text
variableStep chrom=chr21 span=5
9411191	50
9411196	40
9411201	60
9411206	20
9411211	20
```

```text
track type=wiggle_0 name="variableStep" description="variableStep format" visibility=full autoScale=off viewLimits=0.0:25.0 color=50,150,255 yLineMark=11.76 yLineOnOff=on priority=10
variableStep chrom=chr19 span=150
49304701 10.0
49304901 12.5
49305401 15.0
49305601 17.5
49305901 20.0
49306081 17.5
49306301 15.0
49306691 12.5
```

```note
The variableStep format becomes very inefficient when there are only a few data points per 1,024 bases.
```

### Fixed step

```text
fixedStep chrom=chr3 start=400601 step=100 span=5
11
22
33 
```

start, step and span are all optional.