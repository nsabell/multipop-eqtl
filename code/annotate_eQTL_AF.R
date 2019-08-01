#!/usr/bin/env Rscript

library(tidyverse)

setwd("/users/nsabell/multipopeqtl")

#counts = read.table("variant.counts.FDR05.txt")
#counts = counts %>% separate(V2, c("chr","pos"), "_") %>% separate(V3, c("ENSG","gene_name"), "_", extra = "merge")
#names(counts) = c("numPop","chr","pos","ENSG","gene_name")

read_eqtl = function(eqtls) {

	data = read.table(eqtls, header = T)
	data = data[,c(1,4,5)] %>% separate(ENSG, c("ENSG","gene_name"), "_", extra = "merge") %>% separate(SNP, c("chr", "pos"), "_")

}

ceuEQTL = read_eqtl("eQTL/CEU.n69.reMOAT.bestPairs.FDR05.txt")
chbEQTL = read_eqtl("eQTL/CHB.n69.reMOAT.bestPairs.FDR05.txt")
gihEQTL = read_eqtl("eQTL/GIH.n69.reMOAT.bestPairs.FDR05.txt")
jptEQTL = read_eqtl("eQTL/JPT.n69.reMOAT.bestPairs.FDR05.txt")
lwkEQTL = read_eqtl("eQTL/LWK.n69.reMOAT.bestPairs.FDR05.txt")
yriEQTL = read_eqtl("eQTL/YRI.n69.reMOAT.bestPairs.FDR05.txt")

mergeEQTL = list(ceuEQTL, chbEQTL, gihEQTL, jptEQTL, lwkEQTL, yriEQTL) %>% 
	reduce(full_join, by = c("chr","pos","ENSG","gene_name"))
mergeEQTL = mergeEQTL[,c(1,2,4,5,3,6,7,8,9,10)]
names(mergeEQTL) = c("chr","pos","ensg","gene_name","ceu","chb","gih","jpt","lwk","yri")
mergeEQTL$numPop = rowSums(!is.na(mergeEQTL[,5:10]))

read_frequencies = function(frq){

	data = read.table(frq, col.names = paste0("V",seq_len(10)), fill = TRUE, row.names = NULL, header = T)
	data[,5:10] = apply( data[,5:10], 2, function(x) as.numeric(gsub(".*:", "", x, perl = T))  )
	data = data.frame("chr" = data[,1], "pos" = data[,2], "MAF" = apply(data[,5:10], 1, function(x) min(x, na.rm = T)))

}

ceuFreq = read_frequencies("AF/pop/variant.afCEU.FDR05.txt.frq")
chbFreq = read_frequencies("AF/pop/variant.afCHB.FDR05.txt.frq")
gihFreq = read_frequencies("AF/pop/variant.afGIH.FDR05.txt.frq")
jptFreq = read_frequencies("AF/pop/variant.afJPT.FDR05.txt.frq")
lwkFreq = read_frequencies("AF/pop/variant.afLWK.FDR05.txt.frq")
yriFreq = read_frequencies("AF/pop/variant.afYRI.FDR05.txt.frq")

mergeFreq = list(ceuFreq, chbFreq, gihFreq, jptFreq, lwkFreq, yriFreq) %>% 
	reduce(full_join, by = c("chr","pos"))
names(mergeFreq) = c("chr","pos","ceu","chb","gih","jpt","lwk","yri")

#meanAF = merged %>% group_by(numPop) %>% summarize(meanAF = mean(MAF, na.rm = T), stdAF = sd(MAF, na.rm = T))

merged = merge(mergeEQTL, mergeFreq, by = c("chr","pos"))
names(merged) = c("chr","pos","ensg","gene_name","ceu.fdr","chb.fdr","gih.fdr","jpt.fdr","lwk.fdr","yri.fdr","numPop","ceu.af","chb.af","gih.af","jpt.af","lwk.af","yri.af")



