####All Libraries####
library(DESeq2)
library(ggplot2)
library(RColorBrewer)
library(scales)
library(cowplot)
library(VennDiagram)
library(goseq)
library(TxDb.Celegans.UCSC.ce11.ensGene)
library(GO.db)

####Read in female count matrix####
gene.f <- read.table("Documents/UO.academics/Phillips_Lab/experiments/experimental_evolution/2018_FullExperiment/transcriptomics/20220617_Gene_table_counts.txt",
                     sep = "\t", header = FALSE, quote = "", row.names = 1)

coldata.f <- read.table("Dropbox/Experimental_Evolution/transcriptomics/name_headers.txt", sep = "\t", header = TRUE)
coldata.f$treatment <- interaction(coldata.f$female, coldata.f$male, sep = ".")

coldata.f$female <- as.factor(coldata.f$female)
coldata.f$male <- as.factor(coldata.f$male)

z <- c();
for (i in 1:length(coldata.f$female)) {
  if (coldata.f$female[i] == coldata.f$male[i]){
    out <- coldata.f$female[i]
  }
  else {
    if (coldata.f$female[i] == "anc"){
      out <- coldata.f$male[i]
    }
    else {
      out <- coldata.f$female[i]
    }
  }
  z <- rbind(z, as.character(out))
}
coldata.f$EErep <- z
coldata.f$EErep <- as.factor(coldata.f$EErep)

f2 <- c();
for (f in 1:length(coldata.f$female)) {
  if (coldata.f$female[f] == "anc") {
    out <- "anc"
  }
  else {
    out <- "evo"
  }
  f2 <- rbind(f2, out)
}
coldata.f$female2 <- f2
coldata.f$female2 <- as.factor(coldata.f$female2)

m2 <- c();
for (m in 1:length(coldata.f$male)) {
  if (coldata.f$male[m] == "anc") {
    out <- "anc"
  }
  else {
    out <- "evo"
  }
  m2 <- rbind(m2, out)
}
coldata.f$male2 <- m2
coldata.f$male2 <- as.factor(coldata.f$male2)

coldata.f$treatment2 <- interaction(coldata.f$female2, coldata.f$male2)

names(gene.f) <- coldata.f$treatment
gene.f <- gene.f[1:(length(gene.f$anc.A2)-5), ]


#DESeq2 Analysis
#Censor: B3.anc (20_S36_L005Aligned.sortedByCoord.out.ba) and B3.B3 (21_S37_L005Aligned.sortedByCoord.out.ba)
censoredCts <- gene.f[, c(1:10, 13:39)]
censoredCol <- coldata.f[c(1:10, 13:39), ]
names(censoredCts) <- censoredCol$treatment

#Differential expression analysis
ddsF <- DESeqDataSetFromMatrix(countData = censoredCts, colData = censoredCol, design = ~ female2 + male2 + EErep)
ddsF <- DESeq(ddsF)
resultsNames(ddsF)

#results: female anc vs evo
resF1 <- results(ddsF, contrast = list("female2_evo_vs_anc"))
#results: male anc vs evo
resF2 <- results(ddsF, contrast = list("male2_evo_vs_anc"))

#GLM
ctsF <- as.data.frame(counts(ddsF, normalized = TRUE))

check <- c();
pooledFemale <- c();
for (i in 1:length(ctsF$anc.A2)) {
  x <- as.data.frame(t(ctsF[i, ]))
  colnames(x) <- "Gene"
  
  x$treatment <- row.names(x)
  x$female <- substr(x$treatment, start = 1, stop = 3)
  x$male <- substr(x$treatment, start = 4, stop = 6)
  x$EErep <- x$female
  
  x$female[x$female=="A2."] <- "evo"
  x$female[x$female=="A3."] <- "evo"
  x$female[x$female=="B2."] <- "evo"
  x$female[x$female=="B3."] <- "evo"
  
  x$male[x$male==".A2"] <- "evo"
  x$male[x$male==".A3"] <- "evo"
  x$male[x$male==".B2"] <- "evo"
  x$male[x$male==".B3"] <- "evo"
  x$male[x$male=="A2"] <- "evo"
  x$male[x$male=="A3"] <- "evo"
  x$male[x$male=="B2"] <- "evo"
  x$male[x$male=="B3"] <- "evo"
  x$male[x$male=="A2."] <- "evo"
  x$male[x$male=="A3."] <- "evo"
  x$male[x$male=="B2."] <- "evo"
  x$male[x$male=="B3."] <- "evo"
  x$male[x$male==".an"] <- "anc"
  
  x$EErep <- as.factor(x$EErep)
  
  #foldchangeF <- log2(mean(subset(x, female == "evo")$Gene) / mean(subset(x, female == "anc")$Gene))
  #foldchangeM <- log2(mean(subset(x, male == "evo")$Gene) / mean(subset(x, male == "anc")$Gene))
  
  #model <- lmer(Gene ~ female + male + 1|EErep, data = x)
  model <- glm(Gene ~ female + male + EErep, data = x)
  r <- summary(model)
  
  if(length(r$coefficients[,1]) > 6){
    a <- substr(row.names(ctsF)[i], start = 6, stop = 19)
    rbind(check, a)
  }
  else{
    out <- data.frame(WBGeneID  = substr(row.names(ctsF)[i], start = 6, stop = 19),
                      estimateF   = r$coefficients[2,1],
                      estimateM   = r$coefficients[3,1],
                      estimateA3  = r$coefficients[4,1],
                      estimateB2 = r$coefficients[5,1],
                      estimateB3  = r$coefficients[6,1],
                      pvalueF    = r$coefficients[2,4],
                      pvalueM    = r$coefficients[3,4],
                      pvalueA3   = r$coefficients[4,4],
                      pvalueB2  = r$coefficients[5,4],
                      pvalueB3   = r$coefficients[6,4])
    pooledFemale <- rbind(pooledFemale, out)
  }
}

#merge analyses
resF1 <- as.data.frame(resF1)
resF1$WBGeneID <- substr(row.names(resF1), start = 6, stop = 19)
resF1 <- resF1[, c(2, 7)]
names(resF1) <- c("log2FoldChange.resF1", "WBGeneID")

resF2 <- as.data.frame(resF2)
resF2$WBGeneID <- substr(row.names(resF2), start = 6, stop = 19)
resF2 <- resF2[, c(2, 7)]
names(resF2) <- c("log2FoldChange.resF2", "WBGeneID")

pooledFemale <- merge(pooledFemale, resF1, by.x = "WBGeneID", by.y = "WBGeneID")
pooledFemale <- merge(pooledFemale, resF2, by.x = "WBGeneID", by.y = "WBGeneID")


#volcano plot
bonf <- 0.05 / length(pooledFemale$WBGeneID)

femcolors <- c("#FA9FB5", "#F768A1", "#DD3497", "#AE017E","#41AB5D")
plot(pooledFemale$log2FoldChange.resF1, -log10(pooledFemale$pvalueF), las = 1, pch = 19,
     xlim = c(-6, 6), ylim = c(0, 25),
     xlab = "log2FoldChange\nFemale", ylab = "-log10pvalue", main = "Female Contrast",
     col = ifelse(pooledFemale$pvalueF < bonf, alpha("#5AAE61", 0.75), alpha("black", 0.4)))
abline(h = -log10(bonf))
abline(v = 2)
abline(v = -2)

plot(pooledFemale$log2FoldChange.resF2, -log10(pooledFemale$pvalueM), las = 1, pch = 19,
     xlim = c(-6, 6), ylim = c(0, 25),
     xlab = "log2FoldChange\nFemale", ylab = "-log10pvalue", main = "Male Contrast",
     col = ifelse(pooledFemale$pvalueM < bonf, alpha("#5AAE61", 0.75), alpha("black", 0.4)))
abline(h = -log10(bonf))
abline(v = 2)
abline(v = -2)

dev.copy2pdf(file="~/Desktop/2024.07.23_FemaleVolcano_Pooled.pdf", useDingbats=FALSE, family="sans", width = 6, height = 6)

#female overlap of significant genes
sigF <- subset(pooledFemale, pvalueF < bonf)

venn.plotF <- venn.diagram(
  x = list(sigA2.f$WBGeneID, sigA3.f$WBGeneID, sigB2.f$WBGeneID, sigB3.f$WBGeneID, sigF$WBGeneID),
  category.names = c("A2" , "A3" , "B2", "B3", "pooled"),
  filename = NULL,
  fill = alpha(femcolors, 0.5), col = femcolors, lwd = 2, cex = 0.75
)
plot_grid(venn.plotF)
dev.copy2pdf(file="~/Desktop/2024.07.25_Female_PooledSigOverlap.pdf", useDingbats=FALSE, family="sans", width = 3, height = 3)

venn.diagram(
  x = list(sigF$WBGeneID, sigF2$WBGeneID),
  category.names = c("lm" , "dds"),
  filename = "~/Desktop/2024.05.10_FemaleOverlap.png",
  imagetype = "png",
  height = 5.5, width = 6, units = "in",
  fill = alpha(femcolors[1:2], 0.5), col = femcolors[1:2], lwd = 2, cex = 0.75,
  output=TRUE
)





####Read in male count matrix####
gene.m <- read.table("~/Documents/UO.academics/Phillips_Lab/experiments/experimental_evolution/2018_FullExperiment/transcriptomics/20240501_Gene_table_countsM.txt",
                     sep = "\t", header = FALSE, quote = "", row.names = 1)

coldata.m <- read.table("~/Dropbox/Experimental_Evolution/transcriptomics/name_headers_male.txt", sep = "\t", header = TRUE)
coldata.m$treatment <- interaction(coldata.m$female, coldata.m$male, sep = ".")

coldata.m$female <- as.factor(coldata.m$female)
coldata.m$male   <- as.factor(coldata.m$male)

coldata.m$EErep <- coldata.m$male

m2 <- c();
for (m in 1:length(coldata.m$male)) {
  if (coldata.m$male[m] == "anc") {
    out <- "anc"
  }
  else {
    out <- "evo"
  }
  m2 <- rbind(m2, out)
}
coldata.m$male2 <- m2
coldata.m$male2 <- as.factor(coldata.m$male2)

names(gene.m) <- coldata.m$treatment
gene.m <- gene.m[1:(length(gene.m$A2.A2)-5), ]


#DESeq2 Analysis
ddsM <- DESeqDataSetFromMatrix(countData = gene.m, colData = coldata.m, design = ~ male2)
ddsM <- DESeq(ddsM)
resultsNames(ddsM)

#results: male anc vs evo
resM1 <- results(ddsM, contrast = list("male2_evo_vs_anc"))


#GLM
ctsM <- as.data.frame(counts(ddsM, normalized = TRUE))

check <- c();
pooledMale <- c();
for (i in 1:length(ctsM$A2.A2)) {
  x <- as.data.frame(t(ctsM[i, ]))
  colnames(x) <- "Gene"
  
  x$treatment <- row.names(x)
  x$male <- substr(x$treatment, start = 4, stop = 6)
  x$EErep <- substr(x$treatment, start = 1, stop = 3)
  
  x$male[x$male==".A2"] <- "evo"
  x$male[x$male==".A3"] <- "evo"
  x$male[x$male==".B2"] <- "evo"
  x$male[x$male==".B3"] <- "evo"
  x$male[x$male=="A2"] <- "evo"
  x$male[x$male=="A3"] <- "evo"
  x$male[x$male=="B2"] <- "evo"
  x$male[x$male=="B3"] <- "evo"
  x$male[x$male=="A2."] <- "evo"
  x$male[x$male=="A3."] <- "evo"
  x$male[x$male=="B2."] <- "evo"
  x$male[x$male=="B3."] <- "evo"
  x$male[x$male==".an"] <- "anc"
  
  x$EErep <- as.factor(x$EErep)
  
  #foldchangeM <- log2(mean(subset(x, male == "evo")$Gene) / mean(subset(x, male == "anc")$Gene))
  
  #model <- lmer(Gene ~ female + male + 1|EErep, data = x)
  model <- glm(Gene ~ male + EErep, data = x)
  r <- summary(model)
  
  if(length(r$coefficients[,1]) > 5){
    a <- substr(row.names(ctsM)[i], start = 6, stop = 19)
    rbind(check, a)
  }
  else{
    out <- data.frame(WBGeneID  = substr(row.names(ctsM)[i], start = 6, stop = 19),
                      estimateM  = r$coefficients[2,1],
                      estimateA3 = r$coefficients[3,1],
                      estimateB2 = r$coefficients[4,1],
                      estimateB3 = r$coefficients[5,1],
                      pvalueM  = r$coefficients[2,4],
                      pvalueA3 = r$coefficients[3,4],
                      pvalueB2 = r$coefficients[4,4],
                      pvalueB3 = r$coefficients[5,4])
    pooledMale <- rbind(pooledMale, out)
  }
}

#merge analyses
resM1 <- as.data.frame(resM1)
resM1$WBGeneID <- substr(row.names(resM1), start = 6, stop = 19)
resM1 <- resM1[, c(2, 7)]
names(resM1) <- c("log2FoldChange.resM1", "WBGeneID")

pooledMale <- merge(pooledMale, resM1, by.x = "WBGeneID", by.y = "WBGeneID")


#volcano plot
bonf2 <- 0.05 / length(pooledMale$WBGeneID)

malecolors <- c("#9ECAE1", "#6BAED6", "#4292C6", "#2171B5", "#EF3B2C")
plot(pooledMale$log2FoldChange.resM1, -log10(pooledMale$pvalueM), las = 1, pch = 19,
     xlim = c(-8, 8), ylim = c(0, 12),
     xlab = "log2FoldChange", ylab = "-log10pvalue", main = "Male Contrast",
     col = ifelse(pooledMale$pvalueM < bonf2, alpha("#5AAE61", 0.75), alpha("black", 0.4)))
abline(h = -log10(bonf2))
abline(v = 2)
abline(v = -2)

dev.copy2pdf(file="~/Desktop/2024.06.04_MaleVolcano_pooled.pdf", useDingbats=FALSE, family="sans", width = 3, height = 6)

#male overlap of significant genes
sigM <- subset(pooledMale, pvalueM < bonf2)

venn.plotM <- venn.diagram(
  x = list(sigA2.m$WBGeneID, sigA3.m$WBGeneID, sigB2.m$WBGeneID, sigB3.m$WBGeneID, sigM$WBGeneID),
  category.names = c("A2" , "A3" , "B2", "B3", "pooled"),
  filename = NULL,
  fill = alpha(malecolors, 0.5), col = malecolors, lwd = 2, cex = 0.5
)
plot_grid(venn.plotM)
dev.copy2pdf(file="~/Desktop/2024.07.25_Male_PooledSigOverlap.pdf", useDingbats=FALSE, family="sans", width = 3, height = 3)

venn.diagram(
  x = list(sigM$WBGeneID, sigM2$WBGeneID),
  category.names = c("lm" , "dds"),
  filename = "~/Desktop/2024.05.10_MaleOverlap.png",
  imagetype = "png",
  height = 5.5, width = 6, units = "in",
  fill = alpha(malecolors[1:2], 0.5), col = malecolors[1:2], lwd = 2, cex = 0.75,
  output=TRUE
)




####PCA Analysis -- Both Sexes####
#variance stabilizing transformation
vsdF <- vst(ddsF, blind = FALSE)
vsdM <- vst(ddsM, blind = FALSE)

pcaF <- plotPCA(vsdF, intgroup = c("treatment2"), returnData = TRUE)
fem <- ggplot(pcaF, aes(PC1, PC2, color = treatment2)) +
  geom_point(size = 2.5) +
  scale_color_manual(values = c("#9970AB", "#A6DBA0", "#E7D4E8", "#5AAE61")) +
  xlab(paste0("PC1: ", round(100 * attr(pcaF, "percentVar"))[1], "% variance")) +
  ylab(paste0("PC2: ", round(100 * attr(pcaF, "percentVar"))[2], "% variance")) + 
  labs(title = "Female") +
  coord_fixed() +
  theme_bw()

pcaM <- plotPCA(vsdM, intgroup = c("male2"), returnData = TRUE)
mal <- ggplot(pcaM, aes(PC1, PC2, color = male2)) +
  geom_point(size = 2.5) +
  scale_color_manual(values = c("#9970AB", "#5AAE61")) +
  xlab(paste0("PC1: ", round(100 * attr(pcaM, "percentVar"))[1], "% variance")) +
  ylab(paste0("PC2: ", round(100 * attr(pcaM, "percentVar"))[2], "% variance")) + 
  labs(title = "Male") +
  coord_fixed() +
  theme_bw()
plot_grid(mal, fem)

dev.copy2pdf(file="~/Desktop/2024.07.23_PCA_pooled.pdf", useDingbats=FALSE, family="sans", width = 8, height = 3)

####Normalized Counts####
resF1 <- results(ddsF, contrast = list("female2_evo_vs_anc"))
resM1 <- results(ddsM, contrast = list("male2_evo_vs_anc"))

par(mfrow = c(1,2))
plotMA(resF1, alpha = 1e-136, colNonSig = alpha("gray70", 0.5), colSig = femcolors[3], colLine = alpha("black", 0.5), las = 1, ylim = c(-6, 10), main = "Female")
plotMA(resM1, alpha = 1e-49,  colNonSig = alpha("gray70", 0.5), colSig = malecolors[3], colLine = alpha("black", 0.5), las = 1, ylim = c(-6, 10), main = "Male")

dev.copy2pdf(file="~/Desktop/2024.07.25_NormalizedCounts_pooled.pdf", useDingbats=FALSE, family="sans", width = 6, height = 4)


####Shared Analysis####
all <- merge(pooledFemale, pooledMale, by.x = "WBGeneID", by.y = "WBGeneID")

plot(subset(all, pvalueF > bonf & pvalueM.y > bonf2)$log2FoldChange.resF1,
     subset(all, pvalueF > bonf & pvalueM.y > bonf2)$log2FoldChange.resM1,
     xlim = c(-5, 5), ylim = c(-6, 8),
     las = 1, pch = 1, col = alpha("gray70", 0.5), cex = 1.5,
     xlab = "log2FoldChange in Females", ylab = "log2FoldChange in Males")
points(subset(all, pvalueF < bonf & pvalueM.y > bonf2)$log2FoldChange.resF1,
       subset(all, pvalueF < bonf & pvalueM.y > bonf2)$log2FoldChange.resM1,
       pch = 19, col = alpha(femcolors[2], 0.5), cex = 1.5)
points(subset(all, pvalueF > bonf & pvalueM.y < bonf2)$log2FoldChange.resF1,
       subset(all, pvalueF > bonf & pvalueM.y < bonf2)$log2FoldChange.resM1,
       pch = 19, col = alpha(malecolors[2], 0.5), cex = 1.5)
points(subset(all, pvalueF < bonf & pvalueM.y < bonf2)$log2FoldChange.resF1,
       subset(all, pvalueF < bonf & pvalueM.y < bonf2)$log2FoldChange.resM1,
       pch = 19, col = alpha("#6A3D9A", 0.5), cex = 1.5)
abline(h = 0, lty = "dashed")
abline(v = 0, lty = "dashed")

legend("topright", pch = c(1, 19, 19, 19), col = c(alpha("gray70", 0.5), alpha(femcolors[2], 0.5), alpha(malecolors[2], 0.5), alpha("#6A3D9A", 0.5)),
       c("NS", "Female DE", "Male DE", "Shared DE"),
       cex = 0.8, bty = "n")

model1 <- lm(subset(all, pvalueF < bonf & pvalueM.y > bonf2)$log2FoldChange.resM1 ~ subset(all, pvalueF < bonf & pvalueM.y > bonf2)$log2FoldChange.resF1)
x <- summary(model1)
abline(x$coefficients[1,1], x$coefficients[2,1], col = femcolors[2])

model2 <- lm(subset(all, pvalueF < bonf)$log2FoldChange.resM1 ~ subset(all, pvalueF < bonf)$log2FoldChange.resF1)
y <- summary(model2)
abline(y$coefficients[1,1], y$coefficients[2,1], col = femcolors[4])

dev.copy2pdf(file="~/Desktop/2024.07.23_SharedFoldChange.pdf", useDingbats=FALSE, family="sans", width = 5, height = 6)



####Chromosome Mapping####
genes <- read.table("~/Dropbox/Experimental_Evolution/transcriptomics/mart_export.txt", sep = "\t", header = TRUE)

all <- merge(all, genes, by.x = "WBGeneID", by.y = "WBGeneID")

chisq.test(table(subset(all, pvalueF < bonf)$chromosome))
chisq.test(table(subset(all, pvalueM.y < bonf2)$chromosome))

fch <- ggplot(subset(all, pvalueF < bonf), aes(x = chromosome, y = log2FoldChange.resF1,
                            color = ifelse(WBGeneID %in% subset(all, pvalueF < bonf & pvalueM.y < bonf2)$WBGeneID, alpha("#6A3D9A", 0.6), alpha(femcolors[1], 0.5)))) + 
  geom_point(size = 3, position = position_jitter(0.05)) +
  scale_color_identity() +
  ylim(-3, 3) +
  labs(title = "Female") +
  theme_classic()

mch <- ggplot(subset(all, pvalueM.y < bonf2), aes(x = chromosome, y = log2FoldChange.resM1,
                            color = ifelse(WBGeneID %in% subset(all, pvalueF < bonf & pvalueM.y < bonf2)$WBGeneID, alpha("#6A3D9A", 0.6), alpha(malecolors[1], 0.5)))) + 
  geom_point(size = 3, position = position_jitter(0.05)) +
  scale_color_identity() +
  ylim(-5, 3) +
  labs(title = "Male") +
  theme_classic()
plot_grid(fch, mch)
dev.copy2pdf(file="~/Desktop/2024.06.05_SigByChrom.pdf", useDingbats=FALSE, family="sans", width = 5, height = 5)


####Sex-biased Expression####
wtbias <- read.table("~/Documents/UofT/experiments/syntheticSAlocus/Albritton_Genetics2014_SexBiasedGenes.txt", header = TRUE, sep = "\t")

wtbias <- wtbias[, c(1, 4, 5)]
names(wtbias) <- c("WBGeneID", "maleExp", "femExp")

new.all <- merge(all, wtbias, by.x = "WBGeneID", by.y = "WBGeneID")
new.all$sexbias <- new.all$femExp / new.all$maleExp

plot(log2(subset(new.all, pvalueF < bonf)$sexbias), subset(new.all, pvalueF < bonf)$log2FoldChange.resF1, las = 1, pch = 19,
     col = ifelse(subset(new.all, pvalueF < bonf)$chromosome == "X", alpha("#F4A582", 0.75), alpha("#4393C3", 0.75)),
     xlim = c(-10, 10), ylim = c(-4, 2),
     xlab = "Background Female to Male Expression Bias", ylab = "Experimental Evolution Fold Change", main = "Females")
abline(h = 0, lty = "dashed")
abline(v = 0, lty = "dashed")
legend("topleft", col = c(alpha("#F4A582", 0.5), alpha("#4393C3", 0.5)), pch = 19, c("X chromosome", "Autosome"), bty = "n", cex = 0.8)

model2 <- lm(log2FoldChange.resF1 ~ log2(sexbias), data = subset(new.all, pvalueF < bonf & chromosome == "X"))
y <- summary(model2)
abline(y$coefficients[1,1], y$coefficients[2,1], col = "#F4A582")
legend("bottomright", pch = "", c("Female DE on X", round(y$adj.r.squared, digits = 2), y$coefficients[2,4]), cex = 0.8, bty = "n")

plot(log2(subset(new.all, pvalueM.y < bonf2)$sexbias), subset(new.all, pvalueM.y < bonf2)$log2FoldChange.resM1, las = 1, pch = 19,
     col = ifelse(subset(new.all, pvalueM.y < bonf2)$chromosome == "X", alpha("#F4A582", 0.75), alpha("#4393C3", 0.5)),
     xlim = c(-10, 10), ylim = c(-4, 2),
     xlab = "Background Female to Male Expression Bias", ylab = "Experimental Evolution Fold Change", main = "Males")
abline(h = 0, lty = "dashed")
abline(v = 0, lty = "dashed")

model3 <- lm(log2FoldChange.resM1 ~ log2(sexbias), data = subset(new.all, pvalueM.y < bonf2 & sexbias != "Inf"))
z <- summary(model3)
abline(z$coefficients[1,1], z$coefficients[2,1], col = "#053061")
legend("topright", pch = "", c("Male DE", round(z$adj.r.squared, digits = 2), z$coefficients[2,4]), cex = 0.8, bty = "n")

dev.copy2pdf(file="~/Desktop/2024.05.20_SexBias.pdf", useDingbats=FALSE, family="sans", width = 6, height = 4)







####Pathway Analysis####
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("org.Ce.eg.db")

#create gene length matrix
txdb <- TxDb.Celegans.UCSC.ce11.ensGene
txsByGene <- transcriptsBy(txdb, "gene")
lengthData <-median(width(txsByGene))
lengthData <- as.array(lengthData)
names(lengthData) <- substr(names(lengthData), start = 1, stop = 14)

tokeep <- all$WBGeneID[(all$WBGeneID %in% names(lengthData))]

newlengths <- lengthData[(names(lengthData) %in% tokeep)]


#female GO enrichment + KEGG pathways
gene.vectorF <- as.integer(tokeep %in% subset(all, pvalueF < bonf)$WBGeneID)
names(gene.vectorF) <- tokeep

pwf.F <- nullp(gene.vectorF, genome = "ce6", bias.data = newlengths)

GO.F <- goseq(pwf.F, genome = "ce6", "ensGene", use_genes_without_cat = TRUE, test.cats = "GO:MF")
enriched.GO.F <- GO.F$category[p.adjust(GO.F$over_represented_pvalue, method = "BH") < 0.05]

KEGG.F <- goseq(pwf.F, genome = "ce6", "ensGene", use_genes_without_cat = TRUE, test.cats = "KEGG")
enriched.KEGG.F <- KEGG.F$category[p.adjust(KEGG.F$over_represented_pvalue, method = "BH") < 0.05]


#male GO enrichment
gene.vectorM <- as.integer(tokeep %in% subset(all, pvalueM.y < bonf2)$WBGeneID)
names(gene.vectorM) <- tokeep

pwf.M <- nullp(gene.vectorM, genome = "ce6", bias.data = newlengths)

GO.M <- goseq(pwf.M, genome = "ce6", "ensGene", use_genes_without_cat = TRUE, test.cats = "GO:BP")
enriched.GO.M <- GO.M$category[p.adjust(GO.M$over_represented_pvalue, method = "BH") < 0.05]
for(go in enriched.GO.M[1:20]){
  print(GOTERM[[go]])
  cat("--------------------------------------\n")
}

KEGG.M <- goseq(pwf.M, genome = "ce6", "ensGene", use_genes_without_cat = TRUE, test.cats = "KEGG")
enriched.KEGG.M <- KEGG.M$category[p.adjust(KEGG.M$over_represented_pvalue, method = "BH") < 0.05]









####Dosage Compensation####
dc <- read.table("Documents/UO.academics/Phillips_Lab/experiments/experimental_evolution/2018_FullExperiment/transcriptomics/DosageCompensated.txt",
                 header = TRUE, sep = "\t")
notdc <- read.table("Documents/UO.academics/Phillips_Lab/experiments/experimental_evolution/2018_FullExperiment/transcriptomics/NotDC.txt",
                    header = TRUE, sep = "\t")

dc <- dc[, c(1, 3)]
notdc <- notdc[, c(1, 3, 4)]

all_dc <- merge(x = all, y = dc, by = "WBGeneID")
all_notdc <- merge(x = all, y = notdc, by = "WBGeneID")


plot(all_dc$ExpLevel_WT_XX, all_dc$log2FoldChange.resM1, las = 1, xlab = "Normalized Expression in WT XX", ylab = "Expression Change in Evolved Males", bty = "l",
     xlim = c(5, 14), ylim = c(-5, 5),
     pch = ifelse(all_dc$pvalueM.y < bonf2, 19, 1), cex = 1.2,
     col = ifelse(all_dc$pvalueF < bonf & all_dc$pvalueM.y < bonf2, "#6A3D9A", ifelse(all_dc$pvalueM.y < bonf2, alpha(malecolors[3], 0.75), alpha("gray70", 0.5))))
points(all_notdc$ExpLevel_WT_XX, all_notdc$log2FoldChange.resM1,
       pch = ifelse(all_notdc$pvalueM.y < bonf2, 17, 2), cex = 1,
       col = ifelse(all_notdc$pvalueF < bonf & all_notdc$pvalueM.y < bonf2, "#6A3D9A",
                    ifelse(all_notdc$pvalueM.y < bonf2, alpha(malecolors[3], 0.75), alpha("gray70", 0.5))))
abline(h = 0)
abline(h = -mean(notdc$FoldDecrease_XO))


plot(all_dc$ExpLevel_WT_XX, all_dc$log2FoldChange.resF1, las = 1, xlab = "Normalized Expression in WT XX", ylab = "Expression Change in Evolved Females",
     bty = "l", xlim = c(5, 14), ylim = c(-3, 3),
     pch = ifelse(all_dc$pvalueF < bonf, 19, 1), cex = 1.2,
     col = ifelse(all_dc$pvalueF < bonf & all_dc$pvalueM.y < bonf2, "#6A3D9A", ifelse(all_dc$pvalueF < bonf, alpha(femcolors[3], 0.75), alpha("gray70", 0.5))))
points(all_notdc$ExpLevel_WT_XX, all_notdc$log2FoldChange.resF1,
       pch = ifelse(all_notdc$pvalueF < bonf, 17, 2), cex = 1,
       col = ifelse(all_notdc$pvalueF < bonf & all_notdc$pvalueM.y < bonf2, "#6A3D9A",
                    ifelse(all_notdc$pvalueF < bonf, alpha(femcolors[3], 0.75), alpha("gray70", 0.5))))
abline(h = 0)
abline(h = -mean(notdc$FoldDecrease_XO))


plot(all_dc$log2FoldChange.resF1, all_dc$log2FoldChange.resM1, las = 1,
     xlim = c(-3, 3), ylim = c(-5, 5),
     pch = 19, cex = 1.2,
     col = ifelse(all_dc$pvalueF < bonf & all_dc$pvalueM.y < bonf2, "#6A3D9A",
                  ifelse(all_dc$pvalueF < bonf, alpha(femcolors[3], 0.75),
                         ifelse(all_dc$pvalueM.y < bonf2, alpha(malecolors[3]) , alpha("gray70", 0.5)))))
points(all_notdc$log2FoldChange.resF1, all_notdc$log2FoldChange.resM1,
       pch = 17, cex = 1.2,
       col = ifelse(all_notdc$pvalueF < bonf & all_notdc$pvalueM.y < bonf2, "#6A3D9A",
                    ifelse(all_notdc$pvalueF < bonf, alpha(femcolors[3], 0.75),
                           ifelse(all_notdc$pvalueM.y < bonf2, alpha(malecolors[3]) , alpha("gray70", 0.5)))))

dev.copy2pdf(file="~/Desktop/2024.06.10_Dosage.pdf", useDingbats=FALSE, family="sans", width = 6, height = 4)



####Genomic Overlap & cis-SNPs####
#genomic peak overlap
peaks <- read.table("Documents/UO.academics/Phillips_Lab/experiments/experimental_evolution/2018_FullExperiment/transcriptomics/BS_PO_peaks.txt",
                    header = TRUE, sep = "\t")

sharedF <- peaks$WBGeneID %in% subset(all, pvalueF < bonf)$WBGeneID
sharedM <- peaks$WBGeneID %in% subset(all, pvalueM.y < bonf2)$WBGeneID


position <- read.table("Documents/UO.academics/Phillips_Lab/experiments/experimental_evolution/2018_FullExperiment/transcriptomics/genePOS.txt",
                       header = TRUE, sep = "\t")
names(position) <- c("genome", "WBGeneID", "startPOS", "stopPOS", "gene.name")
position <- position[, c(2:5)]

sigFM <- subset(all, pvalueF < bonf | pvalueM.y < bonf2)
sigFM <- merge(sigFM, position, by.x = "WBGeneID", by.y = "WBGeneID")

#cis SNPs
sc <- read.table("Documents/UO.academics/Phillips_Lab/experiments/experimental_evolution/2018_FullExperiment/genomic_data/2020.08.18_GLM_SC_AllAnc_table.txt", header = TRUE, sep = "\t")
sc.genomic <- subset(sc, CHR != "MtDNA")

cisSNPs <- c();
for (i in 1:2) {
  if (dim(subset(sc.genomic, POS < sigFM$startPOS[i] & POS > sigFM$startPOS[i]-2000 & CHR == sigFM$chromosome[i]))[1] != 0){
    out <- subset(sc.genomic, POS < sigFM$startPOS[i] & POS > sigFM$startPOS[i]-2000 & CHR == sigFM$chromosome[i])
    out$ID <- sigFM$WBGeneID[i]
    
    cisSNPs <- rbind(cisSNPs, out)
  }
  else {
      z <- subset(sc.genomic, POS < genes$startPos[i] & POS > genes$startPos[i]-2000 & CHR == genes$CHR[i])
    }
}

sc.bonf <- 0.05 / length(sc.genomic$CHR)

subset(sigFM, WBGeneID == subset(cisSNPs, pvalue < sc.bonf)$ID)
