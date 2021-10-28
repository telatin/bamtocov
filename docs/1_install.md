---
sort: 1
permalink: /installation
---

# Installation

## Install via Miniconda

You can install _bamtocov_ from BioConda, if you have 
[_conda_](https://docs.conda.io/en/latest/miniconda.html) installed:

```bash
conda install -c conda-forge -c bioconda bamtocov
```

## Compiling from source

1. Install [nim](https://nim-lang.org/) and nimble.
1. Install [hts-lib](https://www.htslib.org)
1. Compile with `nimble build`.