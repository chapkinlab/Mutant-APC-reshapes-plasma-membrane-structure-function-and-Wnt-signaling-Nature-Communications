library(here)
#source("http://bioconductor.org/biocLite.R")
#biocLite(c("GenomicFeatures", "AnnotationDbi", "edgeR", "DESeq", "AnnotationFuncs"))
library(edgeR) #load the edgeR package
library(MASS)
library(DESeq)
#############################################################
#Import Data
#############################################################
Counts <- read.csv(here("data","raw","cell_culture-count.T.csv"), header=T)
dim(Counts)
symbl <- Counts[,1]
nam <- colnames(Counts[,-1])
Y.Counts <- as.matrix(Counts[,-1])
rownames(Y.Counts) <- Counts[,1]
#YAMC = 1, IMCE = 2, IMCE + Beta = 3
group <- factor(c(1,1,1,2,2,2,3,3,3))
#############################################################
#  Creating the Design Matrix - Ivan wants upper quartile
#############################################################

y <- DGEList(counts=Y.Counts,group=group,remove.zeros = TRUE)
#removed 17021 rows with all zero counts
design <- model.matrix(~0+group) #GLM model statement
design
#write.csv(design,here("data","Design.csv"))
y <- calcNormFactors(y,method = "upperquartile")
f <- y$samples$lib.size * y$samples$norm.factors/mean(y$samples$lib.size*y$samples$norm.factors)
scaled.counts1 <- round(scale(Y.Counts, center=FALSE, scale=f))
scaled.counts2 <- round(t(t(Y.Counts)/f)*mean(f))
dim(scaled.counts1)
write.csv(scaled.counts2, here("data", "Normalized_data_Monica_cell_culture.csv"))




