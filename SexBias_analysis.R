
####Read in a File S2####
significant <- read.table("PATH TO/FileS2_PooledCensored_SignificantGenes.txt", header = TRUE, sep = "\t")



####Fit linear model of differential expression by sex bias####
model.x <- lm(log2FoldChange ~ log2(FemExp_Allbrittonetal / MaleExp_Albrittonetal), data = subset(significant, CHR == "X"))
summary(model.x)

model.a <- lm(log2FoldChange ~ log2(FemExp_Allbrittonetal / MaleExp_Albrittonetal), data = subset(significant, CHR != "X"))
summary(model.a)



####Read in Albritton et al. 2014####
albritton <-read.table("PATH TO ALBRITTON FOG-2 DATA", header = TRUE, sep = "\t")



####Compare sex-biased expression for significant DE genes###
m <- merge(significnat, albritton, by.x = "WBGeneID", by.y = "WBID")

m$Bias[m$Bias=="LowFem"]    <- "FemBiased"
m$Bias[m$Bias=="LowMale"]   <- "MaleBiased"
m$Bias[m$Bias=="HighFem"]   <- "FemBiased"
m$Bias[m$Bias=="HighMale"]  <- "MaleBiased"
m$Bias[m$Bias=="FemSpec"]   <- "FemBiased"
m$Bias[m$Bias=="MaleSpec"]  <- "MaleBiased"

mat.data <- c(length(subset(m, padj < bonf & Bias == "FemBiased")$WBGeneID),
              length(subset(m, padj < bonf & Bias == "MaleBiased")$WBGeneID),
              length(subset(m, padj < bonf & Bias == "NoBias")$WBGeneID))
bias.matrix <- matrix(mat.data, nrow = 1, ncol = 3, byrow = TRUE)
chisq.test(bias.matrix)
