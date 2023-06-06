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
#  Creating the Design Matrix - determine which normalization is better
#############################################################
library(preprocessCore)
library(RColorBrewer)
par(mfrow=c(2,2))
colors <- brewer.pal(9, "Set1")
d <- DGEList(counts=Y.Counts,group=group,remove.zeros = TRUE)
#removed 23957 rows with all zero counts
design <- model.matrix(~0+group) #GLM model statement
design
#write.csv(design,here("data","Design.csv"))
d <- calcNormFactors(d) #no extra normalization
d <- estimateGLMRobustDisp(d,design) #Roger said this is better, more conservative
d <- as.matrix(d)
boxplot(log10(d+1),lag=2, col=colors, main="No Extra Normalization")
#Boxplot for TMM normalization
d_tmm <- DGEList(counts=Y.Counts,group=group,remove.zeros = TRUE)
#removed 23957 rows with all zero counts
design <- model.matrix(~0+group) #GLM model statement
d_tmm <- calcNormFactors(d_tmm,method="TMM") #TMM normalization (weighted trimmed mean of M-values)
d_tmm <- estimateGLMRobustDisp(d_tmm,design) #Roger said this is better, more conservative
d_tmm <- as.matrix(d_tmm)
boxplot(log10(d_tmm+1),lag=2, col=colors, main="TMM Normalization")
#Boxplot for RLE normalization
d_rle <- DGEList(counts=Y.Counts,group=group,remove.zeros = TRUE)
#removed 23957 rows with all zero counts
design <- model.matrix(~0+group) #GLM model statement
d_rle <- calcNormFactors(d_rle,method="RLE") #RLE normalization (relative log expression)
d_rle <- estimateGLMRobustDisp(d_rle,design) #Roger said this is better, more conservative
d_rle <- as.matrix(d_rle)
boxplot(log10(d_rle+1),lag=2, col=colors, main="RLE Normalization")
#Boxplot for Upper Quartile normalization
d_uq <- DGEList(counts=Y.Counts,group=group,remove.zeros = TRUE)
#removed 23957 rows with all zero counts
design <- model.matrix(~0+group) #GLM model statement
d_uq <- calcNormFactors(d_uq,method="upperquartile") #upper quartile normalization
d_uq <- estimateGLMRobustDisp(d_uq,design) #Roger said this is better, more conservative
d_uq <- as.matrix(d_uq)
boxplot(log10(d_uq+1),lag=2, col=colors, main="Upper Quartile Normalization")



#Boxplots with no log adjustments
boxplot(d,lag=2, col=colors, main="No Extra Normalization")
boxplot(d_tmm,lag=2, col=colors, main="TMM Normalization")
boxplot(d_rle,lag=2, col=colors, main="RLE Normalization")
boxplot(d_uq,lag=2, col=colors, main="Upper Quartile Normalization")


#Am not seeing any differences so will use the default
d <- DGEList(counts=Y.Counts,group=group,remove.zeros = TRUE)
#removed 23957 rows with all zero counts
design <- model.matrix(~0+group) #GLM model statement
design
#write.csv(design,here("data","Design.csv"))
d <- calcNormFactors(d) #no extra normalization
d <- estimateGLMRobustDisp(d,design) #Roger said this is better, more conservative
fit <- glmFit(d,design)
head(fit$coefficients)
head(fit$fitted.values)
values_fit <- fit$fitted.values
write.csv(values_fit, here("data","Fitted_values_all.csv"))
############################################
##Comparisons
###########################################

#compare IMCE vs YAMC
lrt <- glmLRT(fit, contrast=c(-1,1,0))
out <- topTags(lrt,sort.by = "none", n=nrow(Counts))$table
val0 <- (cbind(out,FC=2^({lrt$table[,"logFC"]})))
write.csv(val0,here("data","IMCE_vs_YAMC.csv"))

#compare IMCE + Beta vs YAMC
lrt <- glmLRT(fit, contrast=c(-1,0,1))
out <- topTags(lrt,sort.by = "none", n=nrow(Counts))$table
val0 <- (cbind(out,FC=2^({lrt$table[,"logFC"]})))
write.csv(val0,here("data","IMCE_beta_vs_YAMC.csv"))

#compare IMCE+Beta vs IMCE
lrt <- glmLRT(fit, contrast=c(0,-1,1))
out <- topTags(lrt,sort.by = "none", n=nrow(Counts))$table
val0 <- (cbind(out,FC=2^({lrt$table[,"logFC"]})))
write.csv(val0,here("data","IMCE_beta_vs_IMCE.csv"))

