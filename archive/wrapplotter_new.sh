#!/bin/bash
folder=$1
ph=$2
gene=$3
suffix=$4


sp=/lustre/scratch115/projects/t144_helic_15x/analysis/HA/single_point/release
vcf=/lustre/scratch115/projects/t144_helic_15x/analysis/HA/release


mkdir -p $ph.$gene.$suffix
cd $ph.$gene.$suffix

# Extract chunk number:
chunk=$(grep -lw $gene $folder/Pheno.$ph/MONSTER.*.out | sed 's/.out//;s/.*\.//')
if [ -z "$chunk" ]; then
    echo "[Error] Chunk was not found for $gene. Exiting." >&2
    exit;    
fi
echo Found gene $gene in chunk $chunk.

# Do we want to chop the gene track?
if [[ ! -z "${5}" && "${5}" == "Yes" ]] ; then chop="Yes"; else chop="No"; fi
echo "[Info] Chopping: $chop"

mkdir $ph.$gene
tar -xzf $folder/Pheno.$ph/gene_set.$chunk.tar.gz -C $ph.$gene MONSTER.out snpfile.mod.nomono.txt
cat <(head -1 $ph.$gene/MONSTER.out) <(fgrep -w $gene $ph.$gene/MONSTER.out) > $ph.$gene.MONSTER.out
fgrep -w $gene $ph.$gene/snpfile.mod.nomono.txt > $ph.$gene.snpfile
rm -r $ph.$gene

SRCDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
$SRCDIR/plotburden.py $gene $ph.$gene.snpfile $ph.$gene.MONSTER.out $sp/$ph/MANOLIS.$ph.assoc.txt.bgz $vcf/chr$(wget --no-check-certificate -q -O - "https://rest.ensembl.org/lookup/symbol/homo_sapiens/${gene}?content-type=application/json;expand=0" | sed 's/.*region_name...//;s/\".*//').vcf.gz 100000 ${ph}-${gene}$suffix.html $chop

mv * ..
cd ..
rm -rf $ph.$gene.$suffix
