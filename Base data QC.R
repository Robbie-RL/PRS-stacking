### QC of base file###


# Base data : "Height.QC.gz"
# Target data: "EUR.zip"

setwd("/Users/robryan/Desktop/Masters project/Height stacking model")


# Read in zipped base file 
library(data.table)
dat <- fread("Height.gwas.txt.gz") #529504 individuals, 11 columns


# Create subset where INFO > 0.8 and MAF > 0.01
subset <- dat[INFO > 0.8 & MAF > 0.01] 
fwrite(subset, "Height.gz", sep="\t") #tab deliminated and zipped



# Use plink to remove duplicate and ambiguous SNPs

# Remove duplicate SNPs:

#  gunzip -c Height.gz |\
#  awk '{seen[$3]++; if(seen[$3]==1){ print}}' |\
#  gzip - > Height.nodup.gz

# Remove ambiguous SNPs:

#  gunzip -c Height.nodup.gz |\
#   awk '!( ($4=="A" && $5=="T") || \
#           ($4=="T" && $5=="A") || \
#           ($4=="G" && $5=="C") || \
#           ($4=="C" && $5=="G")) {print}' |\
#   gzip > Height.QC.gz


# Base data is now QC complete: Height.QC.gz









