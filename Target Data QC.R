                         #### QC of target file ###



                             ### Section 1 ###

## Standard GWAS QC
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




                             ### Section 2 ###

## Remove highly correlated SNPs (Linkage disequilibrium)

# ./plink \
#   --bfile EUR \
#   --keep EUR.QC.fam \
#   --extract EUR.QC.snplist \
#   --indep-pairwise 200 50 0.25 \
#   --out EUR.QC

# Output: EUR.QC.prune.in, EUR.QC.prune.out
# EUR.QC.prune.in contains SNPs with LD below specified r^2




                              ### Section 3 ###

## Remove individuals with high or low heterozygosity rates

# Use plink to create .het file which contains F statistic for each individual 
# that compares expected and observed heterozygsity counts

# ./plink \
#   --bfile EUR \
#   --extract EUR.QC.prune.in \
#   --keep EUR.QC.fam \
#   --het \
#   --out EUR.QC

# Output: EUR.QC.het

# Eliminate samples >3 SD from mean
het <- read.table("EUR.QC.het", header=T)
m <- mean(het$F) #mean F
s <- sd(het$F) #std dev. F
valid <- subset(het, F <= m+3*s & F >= m-3*s) #Get samples within 3 SD
#Write valid sample IDs into file
write.table(valid[, c(1,2)], "EUR.valid.sample", quote=F, row.names=F)

#Output: IID/FID sample list of samples within 3SD of mean and thus does not
# have unusualy levels of heterozygosity




                             ### Section 4 ###

## Resolve mismatching SNPs

# Read in target bim file
bim <- read.table("EUR.bim") #location and allele of SNPs
colnames(bim) <- c("CHR", "SNP", "CM", "BP", "B.A1", "B.A2")
# Read in target QCed SNPs (SNP id only)
qc <- read.table("EUR.QC.snplist", header = F, stringsAsFactors = F)
# Read in the base data (SNP list with summary statistics)
height <-
  read.table(gzfile("Height.QC.gz"),
             header = T,
             stringsAsFactors = F, 
             sep="\t")
# Change all alleles to upper case for easy comparison
height$A1 <- toupper(height$A1)
height$A2 <- toupper(height$A2)
bim$B.A1 <- toupper(bim$B.A1)
bim$B.A2 <- toupper(bim$B.A2)

# Identify SNPs that require strand flipping
# Merge base summary statistics with target reported allele
info <- merge(bim, height, by=c("SNP", "CHR", "BP"))
# Only include QC'd SNPs
info <- info[info$SNP %in% qc$V1, ]
# Create function to find complementary allele
complement <- function(x) {
  switch (
    x,
    "A" = "T",
    "C" = "G",
    "T" = "A",
    "G" = "C",
    return(NA)
  )
}
# A -> Base data, B -> Target data
# Get SNPs that have the same alleles across base and target
info.match <- subset(info, A1 == B.A1 & A2 == B.A2)
# Identify SNPs that are complementary between base and target. First create
# what would be the complementary alleles of the target
info$C.A1 <- sapply(info$B.A1, complement)
info$C.A2 <- sapply(info$B.A2, complement)
# Extract SNPs where the base data matches the complement of the target
info.complement <- subset(info, A1 == C.A1 & A2 == C.A2)
# Correct mismatched target SNPs to correct ones
complement.snps <- bim$SNP %in% info.complement$SNP
bim[complement.snps,]$B.A1 <-
  sapply(bim[complement.snps,]$B.A1, complement)
bim[complement.snps,]$B.A2 <-
  sapply(bim[complement.snps,]$B.A2, complement)


# Identify SNPs that require recoding in the target (to ensure the coding 
# allele in the target data is the effective allele in the base summary 
# statistic

# identify SNPs that need recoding
info.recode <- subset(info, A1 == B.A2 & A2 == B.A1)
# Update the recode SNPs
recode.snps <- bim$SNP %in% info.recode$SNP
tmp <- bim[recode.snps,]$B.A1
bim[recode.snps,]$B.A1 <- bim[recode.snps,]$B.A2
bim[recode.snps,]$B.A2 <- tmp

# identify SNPs that need recoding & complement
info.crecode <- subset(info, A1 == C.A2 & A2 == C.A1)
# Update the recode + strand flip SNPs
com.snps <- bim$SNP %in% info.crecode$SNP
tmp <- bim[com.snps,]$B.A1
bim[com.snps,]$B.A1 <- as.character(sapply(bim[com.snps,]$B.A2, complement))
bim[com.snps,]$B.A2 <- as.character(sapply(tmp, complement))

# Output updated bim file
write.table(
  bim[,c("SNP", "B.A1")],
  "EUR.a1",
  quote = F,
  row.names = F,
  col.names = F,
  sep="\t"
)


# Identify SNPs that have different allele in base and target 
# (usually due to difference in genome build or Indel
mismatch <-
  bim$SNP[!(bim$SNP %in% info.match$SNP |
              bim$SNP %in% info.complement$SNP | 
              bim$SNP %in% info.recode$SNP |
              bim$SNP %in% info.crecode$SNP)]
write.table(
  mismatch,
  "EUR.mismatch",
  quote = F,
  row.names = F,
  col.names = F
)




                            ### Section 5 ###

##  Sex and gender mismatch
# The following plink command will generate an F-statistic for each individual
# for the homozygosity of the X chromosome. F < 0.2 = Female, F > 0.8 = Male

# Note that EUR.valid.sample was created in Section 3

# ./plink \
# --bfile EUR \
# --extract EUR.QC.prune.in \
# --keep EUR.valid.sample \
# --check-sex \
# --out EUR.QC

#Output: EUR.QC.sexcheck

# Exclude mismatching sex and gender

# Original: EUR.valid.sample
# Sex info.: EUR.QC.sexcehck

data <- read.table("EUR.valid.sample", header=T)
sex.info <- read.table("EUR.QC.sexcheck", header=T)
data <- subset(sex.info, STATUS=="OK" & FID %in% data$FID)
write.table(data[,c("FID", "IID")], "EUR.QC.valid", 
            row.names=F, col.names=F, sep="\t", quote=F) 

#Output: "EUR.QC.valid" which recently excluded mismatching gender and sex




                         ### Section 6 ###

# Exclude any individuals in the target group that are related

# ./plink \
# --bfile EUR \
# --extract EUR.QC.prune.in \
# --keep EUR.QC.valid \
# --rel-cutoff 0.125 \
# --out EUR.QC

#Output: "EUR.QC.rel.id" 




                          ### Section 7 ###

# Here we create the final target QC which will now make a new .bed file

#./plink \
#  --bfile EUR \
#  --make-bed \
#  --keep EUR.QC.rel.id \
#  --out EUR.QC \
#  --extract EUR.QC.snplist \
#  --exclude EUR.mismatch \
#  --a1-allele EUR.a1

## Final output: "EUR.QC.bed", "EUR.QC.bam", "EUR.QC.fam"



