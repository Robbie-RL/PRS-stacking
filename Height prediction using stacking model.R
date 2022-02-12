##########################################
# Predicting heights using SCT (bigsnpr) #
##########################################

## 4 comparisons
# 1. Standard PRS (default clumping) optimised for best p-value threshold
# 2. Stacked clumping and thresholding
# 3. Best clumping and thresholding combination
# 4. Lassosum




dir <- "/Users/robryan/Desktop/Masters project/Height stacking model"
setwd(dir)


## Data cleaning
base.id <- read.table("EUR.QC.fam") #Sample ID
base.id <- base.id[, c(1,2)]
colnames(base.id) <- c("FID", "IID")
tgt.data <- read.table("EUR.height", header=T) #Sample phenotype (height)
# Remove sample ID without height data
valid <- tgt.data[which(tgt.data$FID %in% base.id$FID), c(1,2)]
# Write into txt file
write.table(valid, file="valid.txt", sep=" ", 
            row.names=FALSE, col.names=FALSE, quote=FALSE)
# From EUR.QC.bed, make new bed file containing valid samples
system(paste0("./plink --bfile EUR.QC \\",
              "--keep valid.txt \\",
              "--make-bed \\",
              "--out EUR.QC.stacking"))
#Output: EUR.QC.stacking.bed + EUR.QC.stacking.bim + EUR.QC.stacking.fam




## Create bigSNP file from QC Base .bed file
snp_readBed("EUR.QC.stacking.bed") #create rds
base.bigSNP <- snp_attach("EUR.QC.stacking.rds")





## Create matching sample height
pheno <- merge(base.bigSNP$fam[c(1,2)], tgt.data, by=c(1,2))




## Useful variables
G <- base.bigSNP$genotypes #[472, 489805]
CHR <- base.bigSNP$map$chromosome #[489805]
POS <- base.bigSNP$map$physical.pos #[489805]
y <- pheno$Height #[472]
NCORES <- nb_cores()




## Load in summary statistics. Merge with base data
sumstats <- read.table(gzfile("Height.QC.gz"), header=T)
sumstats <- sumstats[, c(1, 3, 2, 4, 5, 9, 8)]
# Convert odds ratio to beta
sumstats$OR <- log(sumstats$OR)
# Update and normalise names
names(sumstats) <- c("chr", "rsid", "pos", "a0", "a1", "beta", "p")
map <- base.bigSNP$map[,-(2:3)] #base map data
names(map) <- c("chr", "pos", "a0", "a1") #rename columns to standardise
info_snp <- snp_match(sumstats, map)




# Create BETA and p-value (log) vectors
beta <- rep(NA, ncol(G))
beta[info_snp$`_NUM_ID_`] <- info_snp$beta
lpval <- rep(NA, ncol(G))
lpval[info_snp$`_NUM_ID_`] <- -log10(info_snp$p)




# Set up training/testing split
set.seed(1)
# From 559 individuals, randomly choose 330 (70/30 split)
ind.train <- sample(nrow(G), 330)
# Testing data is everyone not in training data
ind.test <- setdiff(rows_along(G), ind.train)





# Perform clumping using multiple combinations of hyper-parameters
#Note that G includes variants not found in sumstats. The exclude parameter
#will not include them
all_keep <- snp_grid_clumping(G, CHR, POS, ind.row = ind.train,
                              lpS = lpval, exclude = which(is.na(lpval)),
                              ncores = NCORES)
attr(all_keep, "grid")




# Perform thresholding and PRS calculation. Store large matrix in 'loc'
loc <- "/Users/robryan/Desktop/Masters project/Height stacking model/tmp-data-choi/base-data-scores"
multi_PRS <- snp_grid_PRS(G, all_keep, beta, lpval, ind.row = ind.train,
                          backingfile = loc, n_thr_lpS = 50, ncores = NCORES)
dim(multi_PRS) #22chr * 28clumps *50p-val thresholds = 30800 PRS




# Perform stacking linear reg.
final_mod <- snp_grid_stacking(multi_PRS, y[ind.train], ncores = NCORES, K = 4)
summary(final_mod$mod)
str(final_mod, strict.width = "cut")




# Visualise BETA score penalisation degree
new_beta <- final_mod$beta.G
ind <- which(new_beta != 0)
library(ggplot2)
library(bigstatsr)
ggplot(data.frame(y = new_beta, x = beta)[ind, ]) +
  geom_abline(slope = 1, intercept = 0, color = "red") +
  geom_abline(slope = 0, intercept = 0, color = "blue") +
  geom_point(aes(x, y), size = 0.6) +
  theme_bigstatsr() +
  labs(x = "Effect sizes from GWAS", y = "Non-zero effect sizes from SCT")




# Predict affection on testing data
pred <- final_mod$intercept + big_prodVec(G, new_beta[ind],
                                          ind.row=ind.test,
                                          ind.col=ind)
rss <- sum((y[ind.test] - pred)^2) #residual sum of squares
tss <- sum((y[ind.test] - mean(y[ind.test]))^2) #total sum of squares
rsq <- 1 - rss/tss



## Add in population coviariate and sex:
train.cov <- covar_from_df(pheno[ind.train, c(4:10)])
final_mod_cov <- snp_grid_stacking(multi_PRS,
                               y[ind.train], 
                               ncores = NCORES, 
                               K = 4,
                               covar.train=train.cov,
                               pf.covar=c(rep(0, 7))) #covar are not penalised 

summary(final_mod$mod)
str(final_mod, strict.width = "cut")


beta_cov <- final_mod_cov$beta.G #New beta scores
beta.covar <- final_mod_cov$beta.covar #Covariate weightings
ind_cov <- which(new_beta != 0) #Skip any beta scores with value of zero
#Get covariance of testing data
test.cov <- covar_from_df(pheno[ind.test, c(4:10)]) 


#combine new beta scores and covariance into single matrix and convert to FBM
x <- G[,]
x <- cbind(x, pheno[,c(4:10)])
x <- as_FBM(x)

pred.cov <- final_mod_cov$intercept + big_prodVec(x,
                                              c(beta_cov[ind], beta.covar),
                                              ind.row=ind.test,
                                              ind.col=c(ind, 489806:489812))

ssr.cov <- sum((y[ind.test] - pred.cov)^2)
sst.cov <- sum((y[ind.test] - mean(y[ind.test]))^2)
cov.r <- 1 - (ssr.cov / sst.cov)





# Compare standard stacking to introducing population and sex covariates
  #compare option to penalise these covariates

## Add in population coviariate and sex:
# https://www.rdocumentation.org/packages/bigstatsr/versions/1.5.0/topics/big_spLogReg

# Add population structure into FBM of PRS














