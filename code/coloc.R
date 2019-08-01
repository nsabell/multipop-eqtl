setwd("~/multipopeqtl")

library(data.table)
library(tidyverse)
library(coloc)

## Read in GWAS and eQTL files and specify sample sizes

eQTL_all = fread("eQTL/CEU.eQTLAnyPop.n69.reMOAT.CEUsig.txt", header = T)
eQTL_all = separate(eQTL_all, "SNP", c("Chromosome","Position"), "_")

eQTL_leadVars = fread("eQTL/CEU.n69.reMOAT.bestPairs.FDR05.txt", header = T)

GWAS = fread("colocalization/gwas/Kunkle_etal_Stage1_results.txt", header = T)
GWAS$Chromosome = as.character(GWAS$Chromosome)
GWAS$Position = as.character(GWAS$Position)

eQTL_N = 69
GWAS_N = 63926
prop = 0.34386634546

## Read in allele frequencies

AF = fread("AF/AF.eQTL.FDR05.merged.txt", header = F)
names(AF) = c("chr","pos","KG_AF1","KG_AF2","CEU_AF1","CEU_AF2")
AF = transform(AF, KG_MAF = pmin(KG_AF1, KG_AF2), CEU_MAF = pmin(CEU_AF1, CEU_AF2))
AF$chr = as.character(AF$chr)
AF$pos = as.character(AF$pos)

## Construct the coloc inputs

output = data.frame("gene" = character(),"nsnps" = integer(),"PP.H0.abf" = double(),"PP.H1.abf"= double(),"PP.H2.abf"= double(),"PP.H3.abf"= double(),"PP.H4.abf"= double())

for (gene in eQTL_leadVars$ILMN) {
#for (gene in a$ILMN) {

	gene = as.character(gene)
	gene_name = as.character(eQTL_leadVars[which(eQTL_leadVars$ILMN == gene),5])

	index = which(eQTL_all$gene == gene)
	eQTL_subset = eQTL_all[index,]

	merged = merge(eQTL_subset, GWAS, by = c("Chromosome","Position"))
	merged = merged[,c(1,2,4,5,6,11,12,13)]
	merged[,4] = merged$beta / merged$`t-stat`
	names(merged) = c("chr","pos","eQTL_beta","eQTL_se","eQTL_p","GWAS_beta","GWAS_se","GWAS_p")
	merged = merge(merged, AF, by = c("chr","pos"))

	d1 = list(pvalues = as.numeric(merged$GWAS_p), N = GWAS_N, MAF = as.numeric(merged$KG_MAF), type = "cc", s = prop)
	d2 = list(pvalues = as.numeric(merged$eQTL_p), N = eQTL_N, MAF = as.numeric(merged$CEU_MAF), type = "quant")

	print(gene_name)
	results = coloc.abf(d1, d2)
	print("")
	toAdd = data.frame("gene" = gene_name, t(results$summary))
	output = rbind(output, toAdd)

}

write.table(output, "colocalization/Kunkle_coloc_results.txt",quote = F, sep = "\t", row.names = F, col.names = T)





