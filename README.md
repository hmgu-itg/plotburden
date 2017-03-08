# PlotBurden : an interactive tool for visualising results of rare variant burden associations

This tool is aimed at helping researchers explore the results of their gene-based burden testing analyses. It links together results from a burden test, single point association data at the locus and online resources in a single customisable interactive plot. It connects to Ensembl to fetch relevant information about the variants (MAF, rsIDs) and the functional context (overlapping regulatory or genic features, prior evidence of association...). More precisely, it can help answer the following questions:
* Which variants drive my burden?
* Which variants are up/down-weighted?
* What is the LD structure between my variants? 
* Is the burden really driven by a single point signal?
* Do any of the variants recapitulate existing single-point associations?
* Are the driver variants intronic, exonic or regulatory?

## Example
An image speaks better than a thousand words:
![Example image](example.png)
The interactive version can be found [here](example.html).

## Prerequisites
In order to run this software you need to have working copies of the following tools installed:
* Plink 1.9 or newer ([available here](https://www.cog-genomics.org/plink2/index))
* Tabix ([available as part of bcftools/htslib](http://www.htslib.org/download/))

### Python libraries
* `pandas`
* `numpy`
* `bokeh`
* `pybedtools`
* `requests`
* `urllib` (`urllib.request` and `urllib.parse`)


## Installation
Hopefully, this tool should work quasi out-of-the-box. It **needs `tabix` and `plink` to be in your path**, since it calls these directly from inside the script. This is done like:

```bash
export PATH=/path/to/tabix:/path/to/plink:$PATH
```

If you want to make these changes permanent, do:
```bash
echo 'export PATH=/path/to/tabix:/path/to/plink:$PATH' >> ~/.bashrc
```

## Input
So far, the program works only with the output of [GEMMA](http://www.xzlab.org/software.html) for the single-point data and [MONSTER](https://www.stat.uchicago.edu/~mcpeek/software/MONSTER/index.html) for burden testing. Compatibility with more tools can be provided upon request.

## Usage

**Warning** all arguments are required and positional.

```bash
./plotburden.sh [gene_name] [input_monster] [output_monster] [sp_results] [vcf] [window] [output]
```
* **gene_name** is the name of the gene (e.g. _CETP_).
* **input_monster** is the location of the variant file that is given to MONSTER as input.
* **output_monster** : location of the output of the MONSTER run.
* **sp_results** : location of the single-point results file. Expected to be **bgzipped and tabixed**.
* **vcf** : VCF containing source variants. Used for LD calculation.
* **window** : Window, in base pairs, to extend gene boundaries by. If this is too large, the program will be very, very slow. We have experienced good results with 30kb.
* **output** : output file in html format.

## Genome build
The program is currently only compatible with build 38 (GRCh38) of the human genome. Compatibility with b37 can be provided upon request.
