                   #### Base model plink PRS ####

# Necessary files: 
#  Based data:
#    - "Height.QC.gz" 
#  Target data:
#    - "EUR.QC.bed"
#    - "EUR.QC.bim"
#    - "EUR.QC.fam"
#    - "EUR.height"  (target data height)
#    - "EUR.cov" (covariates) (will need to be created)


                   
                   
## Convert Odds Ratio to BETA scores using log transformation
data <- read.table(gzfile("Height.QC.gz"), header=T)
data$BETA <- log(data$OR)
write.table(data, "Height.QC.Transformed", quote = F, row.names = F)




## Clumping: Remove highly correlated SNPs leaving behind a list of 
## relatively independent SNPs. Low p-value SNPs (highly correlated with 
## effect) are prioritized to be preserved.

# Using plink
# ./plink \
# --bfile EUR.QC \
# --clump-p1 1 \
# --clump-r2 0.1 \
# --clump-kb 250 \
# --clump Height.QC.Transformed \
# --clump-snp-field SNP \
# --clump-field P \
# --out EUR

# Output: EUR.clumped

# Extract SNP ID and paste into new file

# Using terminal
# awk 'NR!=1{print $3}' EUR.clumped >  EUR.valid.snp

# Output: EUR.valid.snp
# 193758 SNPs




## Gather necessary PRS files
## Three files are needed:
##  1. Base data: Height.QC.Transformed
##  2. File containing ALL the SNP's IDs and their P-value
##  3. File containing various P-value thresholds


## Generate second file which will be named "SNP.pvalue"

# Using terminal
# awk '{print $3,$8}' Height.QC.Transformed > SNP.pvalue


## Generate third file which will be named "range_list"

# Using terminal
# echo "0.001 0 0.001" > range_list 
# echo "0.05 0 0.05" >> range_list
# echo "0.1 0 0.1" >> range_list
# echo "0.2 0 0.2" >> range_list
# echo "0.3 0 0.3" >> range_list
# echo "0.4 0 0.4" >> range_list
# echo "0.5 0 0.5" >> range_list

# The first column is threshold name, second is lower bound, third is upper.
# Range are inclusive.




## Generate PRS
## plink will generate a PRS using two main functions
##  --score will extract SNP ID ($3), observed genotype ($4) and 
##  SNP effect size ($12)
##  -q-score-range will make sure that --score uses SNPs from certain ranges
##  of SNP p-value since the optimum p-value range is unknown. Recall that 
##  low p-value will increase the chances of over fitting.

#Using plink

# ./plink \
# --bfile EUR.QC \
# --score Height.QC.Transformed 3 4 12 header \
# --q-score-range range_list SNP.pvalue \
# --extract EUR.valid.snp \
# --out EUR

#Output: EUR.0.5.profile to EUR.0.0.001.profile




## Accounting for population structure
## Population structure such as race, gender and age are confounding variables
## and interferes with isolating purely genetic causality.

# Prune target data by eliminating SNPs in high LD
# Using plink:
#./plink \
# --bfile EUR.QC \
# --indep-pairwise 200 50 0.25 \
# --out EUR

#Output: EUR.QC.prune.in, EUR.QC.prune.out

#Calculate first six PC
# Using plink:
# ./plink \
# --bfile EUR.QC \
# --extract EUR.prune.in \
# --pca 6 \
# --out EUR

#Output: EUR.eigenval, EUR.eigenvec




## Find the best P-value threshold and linear regression model

#p-value thresholds to test
p.threshold <- c(0.001, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5)

#Create table of actual heights and covariates for linear regression training
phenotype <- read.table("EUR.height", header=T) #Heights
pcs <-  read.table("EUR.eigenvec", header=F) #Pop. structure PCs
colnames(pcs) <- c("FID", "IID", paste0("PC",1:6)) 
covariate <- read.table("EUR.cov", header=T) #Sex or gender covariate
#Merge all together
pheno <- merge(merge(phenotype, covariate, by=c("FID", "IID")), 
               pcs, by=c("FID", "IID"))

#Create null hypothesis model for linear regression. i.e. Not including 
#PRS calculated earlier. Baseline that only considers population structure
#covariates and other covariates such as sex
null.r2 <- summary(lm(Height~., data = pheno[, !colnames(pheno) %in% 
                                              c("FID", "IID")]))$r.squared


#Include PRS into linear regression and find best threshold
prs.result <- NULL
for (i in p.threshold) {
  pheno.prs <- merge(pheno, 
                     read.table(paste0("EUR.",i,".profile"), 
                                header=T)[,c("FID","IID", "SCORE")],
                     by=c("FID", "IID"))
  model <- summary(lm(Height~., 
                      data=pheno.prs[,!colnames(pheno.prs)%in%c("FID","IID")]))
  model.r2 <- model$r.squared
  prs.r2 <- model.r2 - null.r2 #Performance of a p-value threshold
  print(paste0("EUR.",i,".profile"))
  print(model$coeff)
  prs.coef <- model$coeff["SCORE",] #Linear model 
  prs.result <- rbind(prs.result, 
                      data.frame(Threshold=i, R2=prs.r2, 
                                 P=as.numeric(prs.coef[4]), 
                                 BETA=as.numeric(prs.coef[1]),
                                 SE=as.numeric(prs.coef[2])))
}

# The best p-value threshold
print(prs.result[which.max(prs.result$R2),])


## Example predictions
pheno.best <- merge(pheno, read.table("EUR.0.3.profile", header=T),
                    by=c("FID", "IID"))
model.best <- lm(Height~., 
                 data=pheno.best[, !colnames(pheno.best)%in%c("FID",
                                                              "IID","PHENO", 
                                                              "CNT", "CNT2")])

SSR <- sum((pheno.best$Height - model.best$fitted.values)^2)
SST <- sum((pheno.best$Height - c(rep(mean(pheno.best$Height), 472)))^2)
R.squared <- 1 - (SSR / SST)
summary(model.best)









