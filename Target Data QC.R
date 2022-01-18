### QC of target file###

# Standard GWAS QC
# MAF > 0.01, HWE > 1e-6, genotype rate > 0.01, missingness >0.01


# ./plink \
#   --bfile EUR \
#   --maf 0.01 \
#   --hwe 1e-6 \
#   --geno 0.01 \
#   --mind 0.01 \
#   --write-snplist \
#   --make-just-fam \
#   --out EUR.QC


# Remove highly correlated SNPs (Linkage disequilibrium)

# ./plink \
#   --bfile EUR \
#   --keep EUR.QC.fam \
#   --extract EUR.QC.snplist \
#   --indep-pairwise 200 50 0.25 \
#   --out EUR.QC

# Output: EUR.QC.prune.in, EUR.QC.prune.out
# EUR.QC.prune.in contains SNPs with LD below specified r^2


# Remove individuals with high or low heterozygosity rates

# Use plink to create .het file which contains F statistic for each individual 
# that compares expected and observed heterozygsity counts

# ./plink \
#   --bfile EUR \
#   --extract EUR.QC.prune.in \
#   --keep EUR.QC.fam \
#   --het \
#   --out EUR.QC

#Output EUR.QC.het