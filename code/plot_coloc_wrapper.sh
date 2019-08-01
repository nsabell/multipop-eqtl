#!/bin/bash

Rscript scripts/get_coloc_stats.R $1 $2 $3

for i in $3_*lead* ; do 

	name=`echo $i | cut -d _ -f 1-5` ;
	
	vcftools --gzvcf /users/nsabell/multipopeqtl/genotypes/CEU.chr1-22.n69.p3.MAF01.3A.sorted.vcf.gz --positions $name.txt --hap-r2-positions $i --out $name ;

	Rscript ./scripts/plot_coloc.R $name.list.hap.ld $name.txt $name ;

done