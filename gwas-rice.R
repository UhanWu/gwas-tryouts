library(rrBLUP)
library(BGLR)
library(DT)
library(SNPRelate)
library(tidyverse)
library(qqman)
library(poolr)


############################ READ DATA
rm(list = ls())
Geno <- read_ped("RiceDiversity44k/sativas413.ped")
FAM <- read.table("RiceDiversity44k/sativas413.fam")
MAP <- read.table("RiceDiversity44k/sativas413.map")
rice.pheno <- read.table("http://www.ricediversity.org/data/sets/44kgwas/RiceDiversity_44K_Phenotypes_34traits_PLINK.txt", 
                         header = TRUE, stringsAsFactors = FALSE, sep = "\t")

p = Geno$p
n = Geno$n
Geno = Geno$x


Geno[Geno == 2] <- NA  # Converting missing data to NA
Geno[Geno == 0] <- 0  # Converting 0 data to 0
Geno[Geno == 1] <- 1  # Converting 1 to 1
Geno[Geno == 3] <- 2  # Converting 3 to 2
# Convert the marker data into matrix and transponse and check dimensions
Geno <- matrix(Geno, nrow = p, ncol = n, byrow = TRUE)
Geno <- t(Geno)
dim(Geno)


# See first few columns and rows of the data
rice.pheno[1:5, 1:5]


# datatable(rice.pheno, rownames = FALSE, filter='top', options =
# list(pageLength = 3, scrollX=T)) assign the row names to marker file and
# compare it with phenotypic file
rownames(Geno) <- FAM$V2
table(rownames(Geno) == rice.pheno$NSFTVID)


# Now let us extract the first trait and assign it to object y
y <- matrix(rice.pheno$Flowering.time.at.Arkansas)  # # use the first trait 
rownames(y) <- rice.pheno$NSFTVID
index <- !is.na(y)
y <- y[index, 1, drop = FALSE]  # 374
Geno <- Geno[index, ]  # 374 x 36901
table(rownames(Geno) == rownames(y))

############################ Quality Control check of Marker Data

for (j in 1:ncol(Geno)) {
  Geno[, j] <- ifelse(is.na(Geno[, j]), mean(Geno[, j], na.rm = TRUE), Geno[, 
                                                                            j])
}







