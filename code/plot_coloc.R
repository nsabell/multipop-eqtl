#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

args = c("colocalization/output/plots/Liu_UC_coloc_ENSG00000128604_IRF5.list.hap.ld",
			"colocalization/output/plots/Liu_UC_coloc_ENSG00000128604_IRF5.txt",
			"colocalization/output/plots_DESeq2_allelic/Liu_UC_coloc_ENSG00000128604_IRF5",
			"ENSG00000128604_IRF5",
			"/users/nsabell/multipopeqtl/colocalization/1KG_deseq2_allele_coloc.txt",
			"0.05")

setwd("~/multipopeqtl")

library(data.table)
library(plyr)
library(tidyverse)
library(cowplot)

# Read in LD and colocalization statistics and format
ld = fread(args[1], header =T)
merged = fread(args[2], header = T)
merged = merge(merged, ld, by.x = c("chr","pos"), by.y = c("CHR2","POS2"), all = T)
merged[which(merged$chr == unique(ld$CHR1) & merged$pos == unique(ld$POS1)),"R^2"] = 1
merged$LD = cut(merged$R^2, breaks = c(-0.001,0.2,0.4,0.6,0.8,1))

# Read in finemapped statistics by first grabbing the lead association, then all tests on those probes
eQTL_leadVars = fread("eQTL/CEU.n69.reMOAT.bestPairs.FDR05.txt", header = T)
gene = as.character(unique(eQTL_leadVars[which(eQTL_leadVars$ENSG == args[4]),"ILMN"]))
finemapped = fread("eQTL/finemapped/split.pos.ILMN.reMOAT.CEU.ALL.fisher.txt", header = T)
finemapped = finemapped[which(finemapped$ILMN == gene),]

# Merge the GWAS and eQTL stats with the finemapped stats
merged = merge(merged, finemapped, by = c("chr","pos"), all = T)

# Read in the relevant MPRA results
mpra = read.table(args[5], sep = "\t", header = T)
names(mpra)[16] = "padj"

mpra = mpra[grep(str_split_fixed(args[4], "_", n = 2)[1], mpra$Gene),]
merged_mpra = merge(merged, mpra, by = c("chr","pos"), all= T)

## SECOND PARAMETER = what cutoff to use for coloring
merged_mpra$Hit = "Black"
merged_mpra$Hit[which(merged_mpra$padj < as.numeric(args[6]))] = "Red"

minX = min(merged_mpra$pos[which(!is.na(merged_mpra$eQTL_p))])
maxX = max(merged_mpra$pos[which(!is.na(merged_mpra$eQTL_p))])
merged_mpra = merged_mpra[which(merged_mpra$pos > minX & merged_mpra$pos < maxX),]

merged_mpra$mpraTested = "a"
merged_mpra[which(merged_mpra$padj != "NA"),"mpraTested"] = "b"

#[which(!is.na(GWAS_p)),]
plot2 = ggplot(merged_mpra) + 
		geom_point(aes(x = pos, y = -log10(GWAS_p), colour = Hit, pch = mpraTested)) + 
		geom_point(data = merged_mpra[which(merged_mpra$Hit == "Red")], aes(x = pos, y = -log10(GWAS_p), color = Hit, pch = mpraTested)) + 
		theme_bw() + xlim(c(minX, maxX)) + ylim(c(0, max(-log10(merged_mpra$GWAS_p)))) + scale_color_manual(values = c("Black" = "black","Blue" = "blue","Red" = "red"))

plot3 = ggplot(merged_mpra) + 
		geom_point(aes(x = pos, y = -log10(eQTL_p),  colour = Hit, pch = mpraTested)) +
		geom_point(data = merged_mpra[which(merged_mpra$Hit == "Red")], aes(x = pos, y = -log10(eQTL_p), color = Hit, pch = mpraTested)) +
		theme_bw() + xlim(c(minX, maxX)) + ylim(c(0, max(-log10(merged_mpra$eQTL_p)))) + scale_color_manual(values = c("Black" = "black","Blue" = "blue","Red" = "red"))

plot4 = ggplot(merged_mpra) + 
		geom_point(aes(x = pos, y = -log10(fisher),  colour = Hit, pch = mpraTested)) + 
		geom_point(data = merged_mpra[which(merged_mpra$Hit == "Red")], aes(x = pos, y = -log10(fisher), color = Hit, pch = mpraTested)) +
		theme_bw() + xlim(c(minX, maxX)) + ylim(c(0, max(-log10(merged_mpra$fisher)))) + scale_color_manual(values = c("Black" = "black","Blue" = "blue","Red" = "red"))


pdf(paste0(args[3], ".pdf"), width = 8)

p = plot_grid(plot2, plot3, plot4, align = "v", ncol = 1, axis = "l")
title = ggdraw() + draw_label(str_split_fixed(args[3], "/", n = 4)[4], fontface='bold')
plot_grid(title, p, ncol=1, rel_heights=c(0.1, 1))

dev.off()


