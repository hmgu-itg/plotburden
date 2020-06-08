# PlotBurden : an interactive tool for visualising results of rare variant burden associations

> Warning : This tool is based on the `bokeh` Python plotting library. As the developers do not seem to care much about reverse compatibility, a lot of the features used by `plotburden` are broken as of 0.12.5 (points appearing out of the plotting region, display bugs in the LD blocks, etc.) We therefore strongly advise users to install version 0.12.4 of the library, e.g. with `pip install bokeh==0.12.4`. We will try to make `plotburden` compatible with the latest version when `bokeh` 1.0 is released.

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
The interactive version can be found [here](http://rawgit.com/wtsi-team144/plotburden/master/example.html).

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

### Usage for SMMAT
The scripts were recently upgraded to understand the SMMAT format of burden tests. Please call the python script directly like so:

```bash
plotburden_SMMAT.py [pheno_name] [gene_symbol] [condition] [variant_set_file] [results_file] [sp_results] [vcf] [window] [output_file] [linked_features] [chopping]
```

* **pheno_name** The phenotype name as present in your output file (e.g. _LDL_).
* **gene_symbol** is the name of the gene (e.g. _CETP_). It has to be recognizable and mappable by Ensembl.
* **condition** is the testing condition. The script expects to find variant sets named _ENSGXXX.condition_ in the set file, where _ENSGXXX_ maps to `gene_symbol` in Ensembl.
* **variant_set_file** describes the variant sets defining each burden. It should be of the form `set_name\tchr\tpos\tallele1\tallele2\tweight` (no header).
* **results_file** : is the SMMAT output file. It should be tab delimited with the following columns : `"protein","group","n_variants","miss_min","miss_mean","miss_max","freq_min","freq_mean","freq_max","B_score","B_var","B_pval","S_pval","O_pval","O_minp","O_minp.rho","E_pval"`. **At the moment this file must be bgzip2-ed or uncompressed, not gzipped**.
* **sp_results** : location of the single-point results file. Expected to be **bgzipped and tabixed**.
* **vcf** : VCF containing source variants. Used for LD calculation.
* **window** : Window, in base pairs, to extend gene boundaries by. If this is too large, the program will be very, very slow. We have experienced good results with 30kb.
* **output** : output file in html format. Please add the `.html` extension yourself.
* **linked_features**: linked features file as generated by [MUMMY](https://github.com/hmgu-itg/burden_testing).
* **chopping**: please set this to `False`

## Genome build
The program is currently only compatible with build 38 (GRCh38) of the human genome. Compatibility with b37 can be provided upon request.
