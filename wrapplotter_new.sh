#!/bin/bash
folder=$1 # From which the data will be extracted
ph=$2 # trait
gene=$3 # gene
suffix=$4 # Name of the condition we are extarcting the data


sp=/lustre/scratch115/projects/t144_helic_15x/analysis/HP/single_point/output
vcf=/lustre/scratch115/projects/t144_helic_15x/analysis/HP/release

chunk=$(grep -lw $gene $folder/Pheno.$ph/MONSTER.*.out | sed 's/.out//;s/.*\.//')
if [ -z "$chunk" ]; then
    continue
fi


echo Found gene $gene in chunk $chunk.

mkdir -p $ph.$gene
tar -xzf $folder/Pheno.$ph/gene_set.$chunk.tar.gz -C $ph.$gene MONSTER.out snpfile.mod.nomono.txt gene_set_output_SNPinfo_file.txt
cat <(head -1 $ph.$gene/MONSTER.out) <(fgrep -w $gene $ph.$gene/MONSTER.out) > $ph.$gene.MONSTER.out
cat <(head -n1 $ph.$gene/gene_set_output_SNPinfo_file.txt ) <(fgrep -w $gene  $ph.$gene/gene_set_output_SNPinfo_file.txt ) > ${ph}.${gene}_SNPinfo.txt
fgrep -w $gene $ph.$gene/snpfile.mod.nomono.txt > $ph.$gene.snpfile
rm -r $ph.$gene

SRCDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
$SRCDIR/plotburden.py $gene $ph.$gene.snpfile $ph.$gene.MONSTER.out $sp/$ph/Pomak.$ph.assoc.txt.gz $vcf/chr$(wget --no-check-certificate -q -O - "https://rest.ensembl.org/lookup/symbol/homo_sapiens/${gene}?content-type=application/json;expand=0" | sed 's/.*region_name...//;s/\".*//').reheaded.vcf.gz 100000 ${ph}-${gene}$suffix.html
