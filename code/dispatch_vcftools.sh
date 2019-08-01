#!/bin/bash

module load vcftools

vcftools --vcf /srv/gsfs0/projects/montgomery/mdegorte/data/phase3/n69/$1/$1.chr1-22.n69.p3.MAF01.sorted.vcf --bed  /srv/gsfs0/projects/montgomery/nsabell/mpra/output/$1_intervals.bed --hap-r2-positions  /srv/gsfs0/projects/montgomery/nsabell/mpra/output/$1_topvars.txt --ld-window-bp 100000 --chr $2 --min-r2 0.1 --out /srv/gsfs0/projects/montgomery/nsabell/mpra/output/$1.chr${2} #--min-r2 0.6
