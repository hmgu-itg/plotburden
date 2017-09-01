#!/bin/bash
folder=$1
ph=$2
gene=$3
suffix=$4
sp=/lustre/scratch115/projects/t144_helic_15x/analysis/HA/single_point/output/missing_chunks
vcf=/lustre/scratch115/projects/t144_helic_15x/analysis/HA/release
SRCDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
$SRCDIR/plotburden.py $gene $folder/${gene}_output_variants $folder/MONSTER_out_${gene}.out $sp/$ph/MANOLIS.$ph.assoc.missingfilled.txt.bgz $vcf/chr$(wget --no-check-certificate -q -O - "https://rest.ensembl.org/lookup/symbol/homo_sapiens/${gene}?content-type=application/json;expand=0" | sed 's/.*region_name...//;s/\".*//').vcf.gz 100000 ${ph}-${gene}$suffix.html
