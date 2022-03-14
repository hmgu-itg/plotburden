#!/bin/bash

OIFS=$IFS;
IFS=",";

conames=$1
files=$2
conames=($conames)
files=($files)
IFS=$OIFS

for ((i=0; i<${#conames[@]}; ++i)); do
  tabix -h ${files[$i]} $3 |bgzip> ${conames[$i]}.vcf.gz
  tabix -f -p vcf ${conames[$i]}.vcf.gz
  conames[$i]=${conames[$i]}.vcf.gz
done
bcftools merge $(echo "${conames[*]}") > meta.vcf

rm ${conames[*]}


plink --vcf meta.vcf --r2 dprime --ld-window $4 --ld-window-kb $5 --ld-window-r2 0.1 1>&2 > /dev/null
cat plink.ld
rm plink.ld plink.log plink.nosex
