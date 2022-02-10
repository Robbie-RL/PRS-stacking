##########################
## PRS model comparison ##
##########################


library(bigsnpr)




#################################################
## Load in sample data and relavant variables  ##
#################################################


# Relevant data and variables
base.bigSNP <- snp_attach("EUR.QC.stacking.rds")
G <- base.bigSNP$genotypes #[0,1,2] for each SNP, for each person (matrix)
CHR <- base.bigSNP$map$chromosome #Chr of each SNP (vector)
POS <- base.bigSNP$map$physical.pos #Pos of each SNP (vector)
NCORES <- nb_cores() 




## Set up training/testing split datasets
set.seed(42)
# From 559 individuals, randomly choose 330 (70/30 split)
ind.train <- sample(nrow(G), 330)
# Testing data is everyone not in training data
ind.test <- setdiff(rows_along(G), ind.train)




#####################################################
## Generate principal components from training set ##
#####################################################


## Pruning and clumping of SNPs. Will output a list of SNPs that are relatively
## independent
pruned.list <- snp_clumping(G, CHR, ind.row=ind.train,
                            ncores=NCORES, exclude=snp_indLRLDR(CHR, POS))




## Perform PCA using bigsnpr SVD function and convert to PC 
## eigenvectors using prcomp().
svd <- big_randomSVD(G, snp_scaleBinom(), ncores=NCORES, ind.row=ind.train, 
                     ind.col=pruned.list)
eigenvec.train <- predict(svd, G, ind.row=ind.train, ind.col=pruned.list)




###########################################
## Perform training using stacking model ##
###########################################






