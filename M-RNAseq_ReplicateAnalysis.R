library(DESeq2)
library(ggplot2)
library(RColorBrewer)
library(scales)
library(cowplot)
library(VennDiagram)

####Read in a count matrix####
gene.m <- read.table("~/Documents/UO.academics/Phillips_Lab/experiments/experimental_evolution/2018_FullExperiment/transcriptomics/20240501_Gene_table_countsM.txt",
                     sep = "\t", header = FALSE, quote = "", row.names = 1)

coldata.m <- read.table("~/Dropbox/Experimental_Evolution/transcriptomics/name_headers_male.txt", sep = "\t", header = TRUE)
coldata.m$treatment <- interaction(coldata.m$female, coldata.m$male, sep = ".")

coldata.m$female <- as.factor(coldata.m$female)
coldata.m$male   <- as.factor(coldata.m$male)

coldata.m$EErep <- coldata.m$male

names(gene.m) <- coldata.m$treatment
gene.m <- gene.m[1:(length(gene.m$A2.A2)-5), ]


####Analyze A2 data####
A2coldata <- subset(coldata.m, EErep == "A2" | EErep == "anc")
toMatch   <- c("A2.A2", "anc.anc")
A2cts     <- gene.m[, grep(paste(toMatch,collapse="|"), names(gene.m))]
names(A2cts) <- A2coldata$treatment

#Differential expression analysis
ddsA2 <- DESeqDataSetFromMatrix(countData = A2cts, colData = A2coldata, design = ~ male)
ddsA2 <- DESeq(ddsA2)
ddsA2$male <- relevel(ddsA2$male, ref = "anc")
ddsA2 <-DESeq(ddsA2)
resultsNames(ddsA2)

#result: male only (anc vs A2)
resA2 <- results(ddsA2, contrast = list("male_A2_vs_anc"))

#variance stabilizing transformation
vsdA2 <- vst(ddsA2, blind=FALSE)


####Analyze A3 data####
A3coldata <- subset(coldata.m, EErep == "A3" | EErep == "anc")
toMatch   <- c("A3.A3", "anc.anc")
A3cts     <- gene.m[, grep(paste(toMatch,collapse="|"), names(gene.m))]
names(A3cts) <- A3coldata$treatment

#Differential expression analysis
ddsA3 <- DESeqDataSetFromMatrix(countData = A3cts, colData = A3coldata, design = ~ male)
ddsA3 <- DESeq(ddsA3)
ddsA3$male <- relevel(ddsA3$male, ref = "anc")
ddsA3 <-DESeq(ddsA3)
resultsNames(ddsA3)

#result: male only (anc vs A3)
resA3 <- results(ddsA3, contrast = list("male_A3_vs_anc"))

#variance stabilizing transformation
vsdA3 <- vst(ddsA3, blind=FALSE)


####Analyze B2 data####
B2coldata <- subset(coldata.m, EErep == "B2" | EErep == "anc")
toMatch   <- c("B2.B2", "anc.anc")
B2cts     <- gene.m[, grep(paste(toMatch,collapse="|"), names(gene.m))]
names(B2cts) <- B2coldata$treatment

#Differential expression analysis
ddsB2 <- DESeqDataSetFromMatrix(countData = B2cts, colData = B2coldata, design = ~ male)
ddsB2 <- DESeq(ddsB2)
resultsNames(ddsB2)

#result: male only (anc vs B2)
resB2 <- results(ddsB2, contrast = list("male_B2_vs_anc"))

#variance stabilizing transformation
vsdB2 <- vst(ddsB2, blind=FALSE)


####Analyze B3 data####
B3coldata <- subset(coldata.m, EErep == "B3" | EErep == "anc")
toMatch   <- c("B3.B3", "anc.anc")
B3cts     <- gene.m[, grep(paste(toMatch,collapse="|"), names(gene.m))]
names(B3cts) <- B3coldata$treatment

#Differential expression analysis
ddsB3 <- DESeqDataSetFromMatrix(countData = B3cts, colData = B3coldata, design = ~ male)
ddsB3 <- DESeq(ddsB3)
resultsNames(ddsB3)

#result: male only (anc vs B3)
resB3 <- results(ddsB3, contrast = list("male_B3_vs_anc"))

#variance stabilizing transformation
vsdB3 <- vst(ddsB3, blind=FALSE)



####PCA####
pcaA2 <- plotPCA(vsdA2, intgroup = c("treatment"), returnData = TRUE)
a2 <- ggplot(pcaA2, aes(PC1, PC2, color = treatment)) +
  geom_point(size = 2) +
  scale_color_manual(values = c("#5AAE61", "#9970AB")) +
  xlab(paste0("PC1: ", round(100 * attr(pcaA2, "percentVar"))[1], "% variance")) +
  ylab(paste0("PC2: ", round(100 * attr(pcaA2, "percentVar"))[2], "% variance")) + 
  labs(title = "A2") +
  coord_fixed() +
  theme_bw()

pcaA3 <- plotPCA(vsdA3, intgroup = c("treatment"), returnData = TRUE)
a3 <- ggplot(pcaA3, aes(PC1, PC2, color = treatment)) +
  geom_point(size = 2) +
  scale_color_manual(values = c("#5AAE61", "#9970AB")) +
  xlab(paste0("PC1: ", round(100 * attr(pcaA3, "percentVar"))[1], "% variance")) +
  ylab(paste0("PC2: ", round(100 * attr(pcaA3, "percentVar"))[2], "% variance")) + 
  labs(title = "A3") +
  coord_fixed() +
  theme_bw()

pcaB2 <- plotPCA(vsdB2, intgroup = c("treatment"), returnData = TRUE)
b2 <- ggplot(pcaB2, aes(PC1, PC2, color = treatment)) +
  geom_point(size = 2) +
  scale_color_manual(values = c("#9970AB", "#5AAE61")) +
  xlab(paste0("PC1: ", round(100 * attr(pcaB2, "percentVar"))[1], "% variance")) +
  ylab(paste0("PC2: ", round(100 * attr(pcaB2, "percentVar"))[2], "% variance")) + 
  labs(title = "B2") +
  coord_fixed() +
  theme_bw()

pcaB3 <- plotPCA(vsdB3, intgroup = c("treatment"), returnData = TRUE)
b3 <- ggplot(pcaB3, aes(PC1, PC2, color = treatment)) +
  scale_color_manual(values = c("#9970AB", "#5AAE61")) +
  geom_point(size = 2) +
  xlab(paste0("PC1: ", round(100 * attr(pcaB3, "percentVar"))[1], "% variance")) +
  ylab(paste0("PC2: ", round(100 * attr(pcaB3, "percentVar"))[2], "% variance")) + 
  labs(title = "B3") +
  coord_fixed() +
  theme_bw()

plot_grid(a2, a3, b2, b3)
dev.copy2pdf(file="~/Desktop/2024.07.23_Male_PCAbyRep.pdf", useDingbats=FALSE, family="sans", width = 8, height = 6)


####Significant Differential Expression + Overlap####
par(mfrow = c(4,1))
bonf <- 0.05 / length(ddsA2)

plot(resA2$log2FoldChange, -log10(resA2$padj), las = 1, pch = 19, xlim = c(-10, 10), ylim = c(0, 80),
     xlab = "log2FoldChange\nA2", ylab = "-log10Padj", main = "Male Contrast",
     col = ifelse(-log10(resA2$padj) > -log10(bonf), alpha("#5AAE61", 0.75), alpha("black", 0.4)))
abline(h = -log10(bonf))
abline(v = 2)
abline(v = -2)

plot(resA3$log2FoldChange, -log10(resA3$padj), las = 1, pch = 19, xlim = c(-10, 10), ylim = c(0, 80),
     xlab = "log2FoldChange\nA3", ylab = "-log10Padj",
     col = ifelse(-log10(resA3$padj) > -log10(bonf), alpha("#5AAE61", 0.75), alpha("black", 0.4)))
abline(h = -log10(bonf))
abline(v = 2)
abline(v = -2)

plot(resB2$log2FoldChange, -log10(resB2$padj), las = 1, pch = 19, xlim = c(-10, 10), ylim = c(0, 80),
     xlab = "log2FoldChange\nB2", ylab = "-log10Padj",
     col = ifelse(-log10(resB2$padj) > -log10(bonf), alpha("#5AAE61", 0.75), alpha("black", 0.4)))
abline(h = -log10(bonf))
abline(v = 2)
abline(v = -2)

plot(resB3$log2FoldChange, -log10(resB3$padj), las = 1, pch = 19, xlim = c(-10, 10), ylim = c(0, 80),
     xlab = "log2FoldChange\nB3", ylab = "-log10Padj",
     col = ifelse(-log10(resB3$padj) > -log10(bonf), alpha("#5AAE61", 0.75), alpha("black", 0.4)))
abline(h = -log10(bonf))
abline(v = 2)
abline(v = -2)
dev.copy2pdf(file="~/Desktop/2024.07.23_Male_VolcanoPlotByRep.pdf", useDingbats=FALSE, family="sans", width = 3, height = 8)


#overlap of significant genes
sigA2.m <- as.data.frame(subset(resA2, padj < bonf))
sigA2.m$WBGeneID <- substr(row.names(sigA2.m), start = 6, stop = 19)

sigA3.m <- as.data.frame(subset(resA3, padj < bonf))
sigA3.m$WBGeneID <- substr(row.names(sigA3.m), start = 6, stop = 19)

sigB2.m <- as.data.frame(subset(resB2, padj < bonf))
sigB2.m$WBGeneID <- substr(row.names(sigB2.m), start = 6, stop = 19)

sigB3.m <- as.data.frame(subset(resB3, padj < bonf))
sigB3.m$WBGeneID <- substr(row.names(sigB3.m), start = 6, stop = 19)

malecolors <- c("#9ECAE1", "#6BAED6", "#4292C6", "#2171B5")

venn.plot <- venn.diagram(
  x = list(sigA2.m$WBGeneID, sigA3.m$WBGeneID, sigB2.m$WBGeneID, sigB3.m$WBGeneID),
  category.names = c("A2" , "A3" , "B2", "B3"),
  filename = NULL,
  fill = alpha(malecolors, 0.5), col = malecolors, lwd = 2
)
plot_grid(venn.plot)
dev.copy2pdf(file="~/Desktop/2024.07.25_Male_SigOverlapByRep.pdf", useDingbats=FALSE, family="sans", width = 3, height = 3)




