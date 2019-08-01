#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

setwd("~/multipopeqtl")

library(data.table)
library(tidyverse)
library(coloc)

## Read in global eQTL statistics
eQTL_all = fread("eQTL/CEU.eQTLAnyPop.n69.reMOAT.CEUsig.txt", header = T)
eQTL_all = separate(eQTL_all, "SNP", c("Chromosome","Position"), "_")

# Read in lead eQTL variants
eQTL_leadVars = fread("eQTL/CEU.n69.reMOAT.bestPairs.FDR05.txt", header = T)

print("eQTL input complete")

# Read in GWAS
# fread("zcat colocalization/gwas/GWAS_Type-2-Diabetes_Scott_2017_nonBMIadj.txt.gz", header = T)
GWAS = fread(paste0("zcat ", args[1]), header = T)
GWAS = GWAS[,1:4]
names(GWAS) = c("rsid", "Chromosome","Position","pvalue")
GWAS$Chromosome = as.character(GWAS$Chromosome)
GWAS$Position = as.character(GWAS$Position)

print("GWAS input complete")

# Read in reults file
#results = read.table("colocalization/output/Scott_T2D_coloc_results.txt", header = T)
results = read.table(args[2], header = T)
results = results[which(results$PP.H4.abf >= 0.75),]

# Read in allele frequencies
AF = fread("AF/AF.eQTL.FDR05.merged.MAF.txt", header = T)
AF$chr = as.character(AF$chr)
AF$pos = as.character(AF$pos)

print("Allele frequencies and coloc results input complete")

# Loop through the significant co-localizations
for (gene in results$gene) {

	# Get the test probe and subset the eQTL summary statistics
	gene_name = as.character(eQTL_leadVars[which(eQTL_leadVars$ENSG == gene),2])
	lead_snp = as.character(eQTL_leadVars[which(eQTL_leadVars$ENSG == gene),1])

	lead_file = data.frame("chr" = str_split_fixed(lead_snp, "_", 2)[1], "pos" = str_split_fixed(lead_snp, "_", 2)[2])
	write.table(lead_file, paste0(args[3],"_",as.character(gene),"_lead.txt"), quote = F, sep = "\t", row.names = F, col.names = T)
	#write.table(lead_file, paste0(f,"_",as.character(gene),"_lead.txt"), quote = F, sep = "\t", row.names = F, col.names = T)

	lead_snp = str_split_fixed(lead_snp, "_", 2)[2]

	index = which(eQTL_all$gene == gene_name)
	eQTL_subset = eQTL_all[index,]

	# Merge in GWAS statistics and allele frequencies

	min = min(eQTL_subset$Position)
	max = max(eQTL_subset$Position)
	chr = unique(eQTL_subset$Chromosome)

	GWAS_subset = GWAS[which(GWAS$Position >= min & GWAS$Position <= max & GWAS$chr == chr)]

	merged = merge(eQTL_subset, GWAS, by = c("Chromosome","Position"), all = T)
	merged = merged[,c(1,2,6,9)]
	names(merged) = c("chr","pos","eQTL_p","GWAS_p")
	merged = merge(merged, AF, by = c("chr","pos"))

	#Write the table of positions and run vcftools on it
	write.table(merged, paste0(args[3],"_",as.character(gene),".txt"), quote = F, sep = "\t", row.names = F, col.names = T)

}









