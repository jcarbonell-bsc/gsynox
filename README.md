# GSynoX (a Gene ID translator python package)

![GSynoX logo](docs/logo.png)

## Overview

GSynoX is a lightweight python package designed to simplify a common task in bioinformatics: translating gene IDs between different formats. GSynoX integrates the most popular IDs, such as HGNC, ENSEMBL, and ENTREZ IDs. To face this task efficiently, GSynoX handles the various aliases and deprecated IDs associated with each gene SYMBOL. Leveraging data from the HGNC database, GSynoX enables consistent use of gene symbols across datasets, even when different aliases are employed for the same genes. This is particularly useful for harmonizing data during integrative analyses. Unlike tools like BioMart that perform on-the-fly translations via internet connectivity, GSynoX operates entirely offline. It maintains a local copy of gene ID mappings, ensuring portability and optimal performance in environments such as high-performance computing (HPC) clusters, where internet access may be restricted for security reasons.


## Installation

Open a terminal and type:

```
> git clone https://github.com/jcarbonell-bsc/gsynox.git
> cd gsynox
> pip install .
```

## Getting started

To illustrate how GSynoX works under different scenarios, you can use [this tutorial](docs/getting_started_with_gsynox.md), also available both as a [jupyter notebook](docs/getting_started_with_gsynox.ipynb) and as a plain [python script](docs/getting_started_with_gsynox.py).

## Help and support

For help and support, please send an email to `jcarbonell-bsc@gmail.com` with the subject `GSynoX support`.
