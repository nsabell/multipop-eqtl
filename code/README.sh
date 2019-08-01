###
### Multi-population manuscript scratchpad
###

# Get all positions with a FDR < 0.05 variant in any population
cat eQTL/*FDR05* | grep -v SNP | cut -f 1 | sort | uniq | tr "_" "\t" > variant.positions.FDR05.txt

# Get global 1KG allele frequencies for all those variants
for i in {12..22} ; do
	
	vcftools --positions variant.positions.FDR05.txt \
	--gzvcf /users/nsabell/lab_data_shared/1KG/ALL.chr${i}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz \
	--freq --out AF/global/chr${i}.AF.eQTL.FDR05 &

done

cat AF/global/*frq | grep -v CHROM > AF/variant.afGlobal.FDR05.txt

# Get population specific allele frequencies for all those variants

for file in genotypes/*vcf ; do

	name=`basename $file | cut -d . -f 1`;
	vcftools --positions variant.positions.FDR05.txt \
	--vcf $file --freq --out AF/pop/variant.af${name}.FDR05.txt&

done

### Re-do focusing on just CEU

# Filter CEU all-eQTLs to just those with a significant CEU test
Rscript ./scripts/get_allvars_CEU_eQTLs.R

# Get CEU and 1KG allele frequencies for this table
mkdir AF/1KG
for i in {1..22} ; do
	
	vcftools --positions eQTL/CEU.eQTLAnyPop.n69.reMOAT.CEUsig.positions.txt \
	--gzvcf /users/nsabell/lab_data_shared/1KG/ALL.chr${i}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz \
	--freq --out AF/1KG/chr${i}.AF.eQTL.FDR05 &

done

mkdir AF/CEU
vcftools --positions eQTL/CEU.eQTLAnyPop.n69.reMOAT.CEUsig.positions.txt \
	--vcf genotypes/CEU.chr1-22.n69.p3.MAF01.3A.sorted.vcf \
	--freq --out AF/CEU/CEU.AF.eQTL.FDR05 &

# Assess how many non-biallelic sites there are
cat AF/1KG/*frq | grep -v CHROM | awk '{ print NF}' | sort | uniq -c

    155 10
     39 11
      5 12
      2 13
      1 14
2621111 6
  45436 7
   3443 8
    647 9

# 6 corresponds to biallelic, so we don't lose much by filtering other sites

# Get just biallelic sites

cat AF/1KG/*frq | grep -v CHROM | awk '{ if (NF == 6) print $0 }' | tr ":" "\t" | cut -f 1,2,6,8 > AF/1KG.AF.eQTL.FDR05.txt
cat AF/CEU/*frq | grep -v CHROM | awk '{ if (NF == 6) print $0 }' | tr ":" "\t" | cut -f 1,2,6,8 > AF/CEU.AF.eQTL.FDR05.txt
sort -k 1,1 -k 2,2n 1KG.AF.eQTL.FDR05.txt > c1
sort -k 1,1 -k 2,2n CEU.AF.eQTL.FDR05.txt > c2
paste c1 c2 | cut -f 1,2,3,4,7,8 > AF.eQTL.FDR05.merged.txt

# Now, we have (i) allele frequencies in 1KG, (ii) allele frequencies in CEU,
# (iii) summary stats in CEU including beta, varbeta, pvalue, and N, (iv) GWAS

mkdir -p colocalization/gwas colocalization/output

Rscript scripts/coloc.R > colocalization/Kunkle_coloc_results.log

Rscript scripts/coloc_all.R \
	colocalization/gwas/GWAS_Coronary-Artery-Disease_Nelson_2017.txt.gz \
	148172 \
	0.07862649321 \
	colocalization/output/CAD_nelson_coloc_results.txt > colocalization/output/CAD_nelson_coloc_results.log

Rscript scripts/coloc_all.R \
	colocalization/gwas/GWAS_Inflammatory-Bowel-Disease-European_Liu_2015_CD.txt.gz \
	69268 \
	0.32590806721 \
	colocalization/output/Liu_CD_coloc_results.txt > colocalization/output/Liu_CD_coloc_results.log 2>&1 &

Rscript scripts/coloc_all.R \
	colocalization/gwas/GWAS_Inflammatory-Bowel-Disease-European_Liu_2015_IBD.txt.gz \
	96486 \
	0.44514230043 \
	colocalization/output/Liu_IBD_coloc_results.txt > colocalization/output/Liu_IBD_coloc_results.log 2>&1 &

Rscript scripts/coloc_all.R \
	colocalization/gwas/GWAS_Inflammatory-Bowel-Disease-European_Liu_2015_UC.txt.gz \
	72647 \
	0.28104395226 \
	colocalization/output/Liu_UC_coloc_results.txt > colocalization/output/Liu_UC_coloc_results.log 2>&1 &

Rscript scripts/coloc_all.R \
	colocalization/gwas/GWAS_Type-2-Diabetes_Scott_2017_nonBMIadj.txt.gz \
	159208 \
	0.16755439425 \
	colocalization/output/Scott_T2D_coloc_results.txt > colocalization/output/Scott_T2D_coloc_results.log 2>&1 &

Rscript scripts/coloc_all.R \
	colocalization/gwas/GWAS_Ulcerative-Colitis_Anderson_2011.txt.gz \
	26405 \
	0.253247491 \
	colocalization/output/Anderson_UC_coloc_results.txt > colocalization/output/Anderson_UC_coloc_results.log 2>&1 &

Rscript scripts/coloc_all.R \
	colocalization/gwas/GWAS_Asthma_Moffet_2007_hg19.txt.gz \
	2237 \
	0.44434510505 \
	colocalization/output/Moffat_asthma_coloc_results.txt > colocalization/output/Moffat_asthma_coloc_results.log

# Get summary statistics for colocalization hits

Rscript scripts/get_coloc_stats.R colocalization/gwas/GWAS_Coronary-Artery-Disease_Nelson_2017.txt.gz \
	colocalization/output/CAD_nelson_coloc_results.txt \
	colocalization/output/plots/CAD_nelson_coloc &

Rscript scripts/get_coloc_stats.R colocalization/gwas/GWAS_Inflammatory-Bowel-Disease-European_Liu_2015_CD.txt.gz \
	colocalization/output/Liu_CD_coloc_results.txt \
	colocalization/output/plots/Liu_CD_coloc &

Rscript scripts/get_coloc_stats.R colocalization/gwas/GWAS_Inflammatory-Bowel-Disease-European_Liu_2015_IBD.txt.gz  \
	colocalization/output/Liu_IBD_coloc_results.txt \
	colocalization/output/plots/Liu_IBD_coloc &

Rscript scripts/get_coloc_stats.R colocalization/gwas/GWAS_Inflammatory-Bowel-Disease-European_Liu_2015_UC.txt.gz \
	colocalization/output/Liu_UC_coloc_results.txt \
	colocalization/output/plots/Liu_UC_coloc &

Rscript scripts/get_coloc_stats.R colocalization/gwas/GWAS_Type-2-Diabetes_Scott_2017_nonBMIadj.txt.gz \
	colocalization/output/Scott_T2D_coloc_results.txt \
	colocalization/output/plots/Scott_T2D_coloc &

Rscript scripts/get_coloc_stats.R colocalization/gwas/GWAS_Ulcerative-Colitis_Anderson_2011.txt.gz \
	colocalization/output/Anderson_UC_coloc_results.txt \
	colocalization/output/plots/Anderson_UC_coloc &

Rscript scripts/get_coloc_stats.R colocalization/gwas/GWAS_Asthma_Moffet_2007_hg19.txt.gz \
	colocalization/output/Moffat_asthma_coloc_results.txt \
	colocalization/output/plots/Moffat_asthma_coloc &

# Compute LD for each colocalization with lead variant and plot
# Note - these can be separated to only compute LD once and then refine plotting

for i in colocalization/output/plots/*lead* ; do 
	name=${i/_lead.txt/};
	vcftools --gzvcf /users/nsabell/multipopeqtl/genotypes/CEU.chr1-22.n69.p3.MAF01.3A.sorted.vcf.gz --positions $name.txt --hap-r2-positions $i --out $name &
done

# Plot each set of summary statistics
for i in colocalization/output/plots/*lead* ; do 
	name=${i/_lead.txt/};
	output=${name/plots/plots_DESeq2_allelic};
	gene=`basename $name | cut -d _ -f 4-`;
	echo "Rscript ./scripts/plot_coloc.R $name.list.hap.ld $name.txt $output $gene /users/nsabell/multipopeqtl/colocalization/1KG_deseq2_allele_coloc.txt 0.05 &";
done


for i in colocalization/output/plots/*lead* ; do 
	name=${i/_lead.txt/};
	output=${name/plots/plots_DESeq2_expr};
	gene=`basename $name | cut -d _ -f 4-`;
	Rscript ./scripts/plot_coloc.R $name.list.hap.ld $name.txt $output $gene /users/nsabell/multipopeqtl/colocalization/1KG_deseq2_expr_coloc.txt 0.05 &
done


for i in colocalization/output/plots/*lead* ; do 
	name=${i/_lead.txt/};
	output=${name/plots/plots_betaReg};
	gene=`basename $name | cut -d _ -f 4-`;
	Rscript ./scripts/plot_coloc.R $name.list.hap.ld $name.txt $output $gene /users/nsabell/multipopeqtl/colocalization/1KG_betaReg_coloc.txt 0.05 &
done




