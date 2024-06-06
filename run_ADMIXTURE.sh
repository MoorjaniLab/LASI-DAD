#!/bin/sh
##Script to prepare the file for and run ADMIXTURE software, can be download here : https://dalexander.github.io/admixture/download.html

##file in vcf format
file=$1

##freq for the MAF
freq=0.05

##Transform VCF to plink format + LD pruning and MAF selection
# perform linkage pruning - i.e. identify prune sites
plink2 --vcf $file --indep-pairwise 50 10 0.5 --set-missing-var-ids @:# --snps-only   --out plink_${file}
plink2 --vcf $file --maf ${freq} --extract plink_${freq}.prune.in --set-missing-var-ids @:#  --make-bed --out pruned05_maf${freq}_${file}



for K in {2..8}
do
  #run ADMIXTURE for different K
  admixture_linux-1.3.0/admixture --cv=10 pruned05_maf${freq}_${file}".bed" ${K} -j20  -C 0.1

done

