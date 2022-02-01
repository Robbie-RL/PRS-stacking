             ### Bigsnpr stacking model tutorial ###

library(bigsnpr)




# Create bigSNP file and object
snp_readBed("/Users/robryan/Desktop/Masters project/Height stacking model/tmp-data/public-data.bed")

# Read in bigSNP object file
obj.bigSNP <- snp_attach("/Users/robryan/Desktop/Masters project/Height stacking model/tmp-data/public-data.rds")

# Main attributes are: $genotypes, $fam, $map
# 559 individuals
str(obj.bigSNP, max.level = 2, strict.width = "cut")

# Create useful variables
G   <- obj.bigSNP$genotypes #Genotype for set of var (559x130816)
CHR <- obj.bigSNP$map$chromosome #Chromosome of set of var (130816)
POS <- obj.bigSNP$map$physical.pos #var positions (129816)
y   <- obj.bigSNP$fam$affection - 1 #0:unaffected, #1:affected (559)
NCORES <- nb_cores()




#Example counts for the first 10 variants. Number of non-reference alleles
big_counts(G, ind.col = 1:10)

#Example counts for last 10 variants. Number of non-reference alleles
big_counts(G, ind.col = 130806:130816)

#Example summary statistics columns 1 to 6 and 10 and store as object
sumstats <- bigreadr::fread2("/Users/robryan/Desktop/Masters project/Height stacking model/tmp-data/public-data-sumstats.txt", 
                       select = c(1:6, 10))




# Set up training
set.seed(1)
# From 559 individuals, randomly choose 400 
ind.train <- sample(nrow(G), 400)
# Testing data is everyone not in training data
ind.test <- setdiff(rows_along(G), ind.train)



# Rename columns of sumstats
names(sumstats) <- c("chr", "rsid", "pos", "a0", "a1", "beta", "p")

# Exclude $marker.ID and $genetic.dist
map <- obj.bigSNP$map[,-(2:3)] 
names(map) <- c("chr", "pos", "a0", "a1") #rename columns to standardise

# Merge matching genotype data (map) and summary statistics (sumstats)
info_snp <- snp_match(sumstats, map)

# Ambiguous SNPs check
filter.out <- setdiff(sumstats$rsid, info_snp$rsid)
filter.out <- data.frame(filter.out)
colnames(filter.out) <- c("rsid")
filtered.out <- merge(sumstats[, -c(1, 3)], filter.out, by=c("rsid"))
filtered.out[1:10, ]





# Create BETA and p-value (log) data frame
beta <- rep(NA, ncol(G))
beta[info_snp$`_NUM_ID_`] <- info_snp$beta
lpval <- rep(NA, ncol(G))
lpval[info_snp$`_NUM_ID_`] <- -log10(info_snp$p)





# Perform clumping using multiple combinations of hyper-parameters
all_keep <- snp_grid_clumping(G, CHR, POS, ind.row = ind.train,
                              lpS = lpval, exclude = which(is.na(lpval)),
                              ncores = NCORES)
attr(all_keep, "grid")





# Perform thresholding and PRS calculation
loc <- "/Users/robryan/Desktop/Masters project/Height stacking model/tmp-data/public-data-scores"
multi_PRS <- snp_grid_PRS(G, all_keep, beta, lpval, ind.row = ind.train,
                          backingfile = loc, n_thr_lpS = 50, ncores = NCORES)
dim(multi_PRS)
multi_PRS




# Stacked C+T
# Phenotypes (response variable) is affection
final_mod <- snp_grid_stacking(multi_PRS, y[ind.train], ncores = NCORES, K = 4)
summary(final_mod$mod)
str(final_mod, strict.width = "cut")





# Compare summary stat BETA scores weights to SCT score weights
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




