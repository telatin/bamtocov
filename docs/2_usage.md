---
sort: 2
permalink: /usage
---
# General usage

`qax` is composed by five subprogram, and the general syntax is:

```
qax [program] parameters
```

The programs are:
- **list** (it's the default action and can be omitteed)
- **extract** or **x**
- **citations** or **c**
- **provenance** or **p**
- **view** or **v**

## list



This is the default module, and can be used to list the properties of one or more artifacts.

Some features:
* Supports multiple files at once
* 100X times faster than Qiime2
* Can be used to find an artifact given the ID

Example:
```
qax_mac -b -u input/*.*
┌───────────────────────────┬────────────────┬─────────────────────────┬─────────────────────────────┐
│ ID                        │ Basename       │ Type                    │ Format                      │
├───────────────────────────┼────────────────┼─────────────────────────┼─────────────────────────────┤
│ bb1b2e93-...-2afa2110b5fb │ rep-seqs.qza   │ FeatureData[Sequence]   │ DNASequencesDirectoryFormat │
│ 313a0cf3-...-befad4ebf2f3 │ table.qza      │ FeatureTable[Frequency] │ BIOMV210DirFmt              │
│ 35c32fe7-...-85ef27545f00 │ taxonomy.qzv   │ Visualization           │ HTML                        │
└───────────────────────────┴────────────────┴─────────────────────────┴─────────────────────────────┘
```

## extract



This program extract the content of an artifact. By default, if a single file is present it will be extracted in the specified path. If multiple files are present, a directory containing them will be created instead.

_Example:_
```
# Extract representative sequences (will be called rep-seqs.fasta)
qax x -o ./ rep-seqs.qza

# Extract a visualization (a folder called "taxonomy" will be created)
qax x -o ./ taxonomy.qzv
```

## citations

Each Qiime module provides the citations for the software and resources that it uses, storing the citations in BibTeX format inside the artifacts. The cite module allows to extract all the citations from a list of artifacts, removing the duplicates, thus effectively allowing to prepare the bibliography for a complete Qiime2 analysis.

_Example:_
```
qax c files/*.qza > bibliography.bib
```

## provenance

This program allows to print the provenance of an artifact, or to produce a [publication grade graph](docs/qax-provenance.png) of the provenance.

_Example:_
```
# To view a summary
qax p taxonomy.qzv

# To save the plot
qax p -o graph.dot taxonomy.qza
```


## view

This program allows to print the content of an artifact data file to the terminal.
If the artifact contains a single file, it will be printed. Otherwise the user can specify one or multiple files to be printed, and if none
is specified, a list of files will be printed.

```
# Example: count the number of representative sequences
qax view rep-seqs.qza | grep -c '>'
```

## make

This program converts a directory containing a website (index.html) into a
_visualization_ artifact.

```
qax make -o report.qzv /path/to/webpage/
```
