#!/bin/bash
tabix -h $1 $2 > tmp.vcf
plink --vcf tmp.vcf --r2 dprime --ld-window $3 --ld-window-kb $4 --ld-window-r2 0.1 2>&1 > /dev/null
cat plink.ld
rm tmp.vcf plink.ld
