library(data.table)
library(tidyverse)

a = fread("eQTL/CEU.n69.reMOAT.bestPairs.FDR05.col2.txt", header = F)
b = fread("eQTL/CEU.eQTLAnyPop.n69.reMOAT.txt")

c = b[which(b$gene %in% a$V1),]
write.table(c, "eQTL/CEU.eQTLAnyPop.n69.reMOAT.CEUsig.txt", quote = F, sep = "\t", row.names = F)

d = separate(c, col = "SNP", into = c("chr","pos"), sep = "_")
write.table(d[,1:2], "CEU.eQTLAnyPop.n69.reMOAT.CEUsig.positions.txt",col.names=F, row.names=F, quote = F, sep = "\t")
