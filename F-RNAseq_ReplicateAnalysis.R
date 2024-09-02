library(DESeq2)
library(ggplot2)
library(RColorBrewer)
library(scales)
library(cowplot)
library(VennDiagram)

####Read in a count matrix####
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

names(gene.f) <- coldata.f$treatment
gene.f <- gene.f[1:(length(gene.f$anc.A2)-5), ]


####Analyze A2 data####
A2coldata <- subset(coldata.f, EErep == "A2" | EErep == "anc")
toMatch   <- c("anc.anc", "anc.A2", "A2.anc", "A2.A2")
A2cts     <- gene.f[, grep(paste(toMatch,collapse="|"), names(gene.f))]
names(A2cts) <- A2coldata$treatment

#Differential expression analysis
ddsA2 <- DESeqDataSetFromMatrix(countData = A2cts, colData = A2coldata, design = ~ female + male + female:male)
ddsA2 <- DESeq(ddsA2)
ddsA2$female <- relevel(ddsA2$female, ref = "anc")
ddsA2$male <- relevel(ddsA2$male, ref = "anc")
ddsA2 <-DESeq(ddsA2)
resultsNames(ddsA2)

#results: female only (anc vs A2)
resA2.1 <- results(ddsA2, contrast = list("female_A2_vs_anc"))
#results: male only (anc vs A2)
resA2.2 <- results(ddsA2, contrast = list("male_A2_vs_anc"))
#results: interaction
resA2.3 <- results(ddsA2, contrast = list("femaleA2.maleA2"))

#variance stabilizing transformation
vsdA2 <- vst(ddsA2, blind=FALSE)


####Analyze A3 data####
A3coldata <- subset(coldata.f, EErep == "A3" | EErep == "anc")
toMatch   <- c("anc.anc", "anc.A3", "A3.anc", "A3.A3")
A3cts     <- gene.f[, grep(paste(toMatch,collapse="|"), names(gene.f))]
names(A3cts) <- A3coldata$treatment

#Differential expression analysis
ddsA3 <- DESeqDataSetFromMatrix(countData = A3cts, colData = A3coldata, design = ~ female + male + female:male)
ddsA3 <- DESeq(ddsA3)
ddsA3$female <- relevel(ddsA3$female, ref = "anc")
ddsA3$male <- relevel(ddsA3$male, ref = "anc")
ddsA3 <-DESeq(ddsA3)
resultsNames(ddsA3)

#results: female only (anc vs A3)
resA3.1 <- results(ddsA3, contrast = list("female_A3_vs_anc"))
#results: male only (anc vs A3)
resA3.2 <- results(ddsA3, contrast = list("male_A3_vs_anc"))
#results: interaction
resA3.3 <- results(ddsA3, contrast = list("femaleA3.maleA3"))

#variance stabilizing transformation
vsdA3 <- vst(ddsA3, blind=FALSE)


####Analyze B2 data####
B2coldata <- subset(coldata.f, EErep == "B2" | EErep == "anc")
toMatch   <- c("anc.anc", "anc.B2", "B2.anc", "B2.B2")
B2cts     <- gene.f[, grep(paste(toMatch,collapse="|"), names(gene.f))]
names(B2cts) <- B2coldata$treatment

#Differential expression analysis
ddsB2 <- DESeqDataSetFromMatrix(countData = B2cts, colData = B2coldata, design = ~ female + male + female:male)
ddsB2 <- DESeq(ddsB2)
resultsNames(ddsB2)

#results: female only (anc vs B2)
resB2.1 <- results(ddsB2, contrast = list("female_B2_vs_anc"))
#results: male only (anc vs A3)
resB2.2 <- results(ddsB2, contrast = list("male_B2_vs_anc"))
#results: interaction
resB2.3 <- results(ddsB2, contrast = list("femaleB2.maleB2"))

#variance stabilizing transformation
vsdB2 <- vst(ddsB2, blind=FALSE)


####Analyze B3 data####
B3coldata <- subset(coldata.f, EErep == "B3" | EErep == "anc")
toMatch   <- c("anc.anc", "anc.B3", "B3.anc", "B3.B3")
B3cts     <- gene.f[, grep(paste(toMatch,collapse="|"), names(gene.f))]

#Differential expression analysis
ddsB3 <- DESeqDataSetFromMatrix(countData = B3cts, colData = B3coldata, design = ~ female + male + female:male)
ddsB3 <- DESeq(ddsB3)
resultsNames(ddsB3)

#variance stabilizing transformation
vsdB3 <- vst(ddsB3, blind=FALSE)

#Censor: B3.anc (20_S36_L005Aligned.sortedByCoord.out.ba) and B3.B3 (21_S37_L005Aligned.sortedByCoord.out.ba)
B3ctsCensor        <- B3cts[, c(1:3, 6:12)]
B3coldataCensor    <- B3coldata[c(1:3, 6:12),]
names(B3ctsCensor) <- B3coldataCensor$treatment

#Differential expression analysis
ddsB3C <- DESeqDataSetFromMatrix(countData = B3ctsCensor, colData = B3coldataCensor, design = ~ female + male + female:male)
ddsB3C <- DESeq(ddsB3C)
resultsNames(ddsB3C)

#results: female only (anc vs B3)
resB3.1C <- results(ddsB3C, contrast = list("female_B3_vs_anc"))
#results: male only (anc vs B3)
resB3.2C <- results(ddsB3C, contrast = list("male_B3_vs_anc"))
#results: interaction
resB3.3C <- results(ddsB3C, contrast = list("femaleB3.maleB3"))

#variance stabilizing transformation
vsdB3C <- vst(ddsB3C, blind=FALSE)



####PCA####
pcaB3full <- plotPCA(vsdB3, intgroup = c("treatment"), returnData = TRUE)
b3full <- ggplot(pcaB3full, aes(PC1, PC2, color = treatment)) +
  geom_point(size = 2) +
  scale_color_manual(values = c("#9970AB", "#A6DBA0", "#E7D4E8", "#5AAE61")) +
  xlab(paste0("PC1: ", round(100 * attr(pcaB3full, "percentVar"))[1], "% variance")) +
  ylab(paste0("PC2: ", round(100 * attr(pcaB3full, "percentVar"))[2], "% variance")) + 
  labs(title = "B3 full") +
  coord_fixed() +
  theme_bw()

pcaB3cens <- plotPCA(vsdB3C, intgroup = c("treatment"), returnData = TRUE)
b3cens <- ggplot(pcaB3cens, aes(PC1, PC2, color = treatment)) +
  geom_point(size = 2) +
  scale_color_manual(values = c("#9970AB", "#A6DBA0", "#E7D4E8", "#5AAE61")) +
  xlab(paste0("PC1: ", round(100 * attr(pcaB3cens, "percentVar"))[1], "% variance")) +
  ylab(paste0("PC2: ", round(100 * attr(pcaB3cens, "percentVar"))[2], "% variance")) + 
  labs(title = "B3 censor") +
  coord_fixed() +
  theme_bw()
plot_grid(b3full, b3cens)
dev.copy2pdf(file="~/Desktop/2024.07.23_Female_B3_CensorComparison.pdf", useDingbats=FALSE, family="sans", width = 6, height = 3)


pcaA2 <- plotPCA(vsdA2, intgroup = c("treatment"), returnData = TRUE)
a2 <- ggplot(pcaA2, aes(PC1, PC2, color = treatment)) +
  geom_point(size = 2) +
  scale_color_manual(values = c("#9970AB", "#A6DBA0", "#E7D4E8", "#5AAE61")) +
  xlab(paste0("PC1: ", round(100 * attr(pcaA2, "percentVar"))[1], "% variance")) +
  ylab(paste0("PC2: ", round(100 * attr(pcaA2, "percentVar"))[2], "% variance")) + 
  labs(title = "A2") +
  coord_fixed() +
  theme_bw()

pcaA3 <- plotPCA(vsdA3, intgroup = c("treatment"), returnData = TRUE)
a3 <- ggplot(pcaA3, aes(PC1, PC2, color = treatment)) +
  geom_point(size = 2) +
  scale_color_manual(values = c("#9970AB", "#A6DBA0", "#E7D4E8", "#5AAE61")) +
  xlab(paste0("PC1: ", round(100 * attr(pcaA3, "percentVar"))[1], "% variance")) +
  ylab(paste0("PC2: ", round(100 * attr(pcaA3, "percentVar"))[2], "% variance")) + 
  labs(title = "A3") +
  coord_fixed() +
  theme_bw()
  
pcaB2 <- plotPCA(vsdB2, intgroup = c("treatment"), returnData = TRUE)
b2 <- ggplot(pcaB2, aes(PC1, PC2, color = treatment)) +
  geom_point(size = 2) +
  scale_color_manual(values = c("#9970AB", "#A6DBA0", "#E7D4E8", "#5AAE61")) +
  xlab(paste0("PC1: ", round(100 * attr(pcaB2, "percentVar"))[1], "% variance")) +
  ylab(paste0("PC2: ", round(100 * attr(pcaB2, "percentVar"))[2], "% variance")) + 
  labs(title = "B2") +
  coord_fixed() +
  theme_bw()

plot_grid(a2, a3, b2, b3cens)
dev.copy2pdf(file="~/Desktop/2024.07.23_Female_PCAbyRep.pdf", useDingbats=FALSE, family="sans", width = 8, height = 6)


####Significant Differential Expression + Overlap####
par(mfrow = c(4,2))
bonf <- 0.05 / length(ddsA2)
r <- c("Female Contrast", "Male Contrast")
for (i in 1:2) {
  plot(get(paste0("resA2.", i))$log2FoldChange, -log10(get(paste0("resA2.", i))$padj), las = 1, pch = 19, xlim = c(-6, 8), ylim = c(0, 50),
       xlab = "log2FoldChange\nA2", ylab = "-log10Padj", main = r[i],
       col = ifelse(-log10(get(paste0("resA2.", i))$padj) > -log10(bonf),alpha("#5AAE61", 0.75), alpha("black", 0.4)))
  abline(h = -log10(bonf))
  abline(v = 2)
  abline(v = -2)
}

for (i in 1:2) {
  plot(get(paste0("resA3.", i))$log2FoldChange, -log10(get(paste0("resA3.", i))$padj), las = 1, pch = 19, xlim = c(-6, 8), ylim = c(0, 50),
       xlab = "log2FoldChange\nA3", ylab = "-log10Padj",
       col = ifelse(-log10(get(paste0("resA3.", i))$padj) > -log10(bonf),alpha("#5AAE61", 0.75), alpha("black", 0.4)))
  abline(h = -log10(bonf))
  abline(v = 2)
  abline(v = -2)
}

for (i in 1:2) {
  plot(get(paste0("resB2.", i))$log2FoldChange, -log10(get(paste0("resB2.", i))$padj), las = 1, pch = 19, xlim = c(-6, 8), ylim = c(0, 50),
       xlab = "log2FoldChange\nB2", ylab = "-log10Padj",
       col = ifelse(-log10(get(paste0("resB2.", i))$padj) > -log10(bonf),alpha("#5AAE61", 0.75), alpha("black", 0.4)))
  abline(h = -log10(bonf))
  abline(v = 2)
  abline(v = -2)
}

for (i in 1:2) {
  plot(get(paste0("resB3.", i, "C"))$log2FoldChange, -log10(get(paste0("resB3.", i, "C"))$padj), las = 1, pch = 19, xlim = c(-6, 8), ylim = c(0, 50),
       xlab = "log2FoldChange\nB3", ylab = "-log10Padj",
       col = ifelse(-log10(get(paste0("resB3.", i, "C"))$padj) > -log10(bonf),alpha("#5AAE61", 0.75), alpha("black", 0.4)))
  abline(h = -log10(bonf))
  abline(v = 2)
  abline(v = -2)
}
dev.copy2pdf(file="~/Desktop/2024.07.23_Female_VolcanoPlotsByRep.pdf", useDingbats=FALSE, family="sans", width = 6, height = 8)


#overlap of significant genes
sigA2.f <- as.data.frame(subset(resA2.1, padj < bonf))
sigA2.f$WBGeneID <- substr(row.names(sigA2.f), start = 6, stop = 19)

sigA3.f <- as.data.frame(subset(resA3.1, padj < bonf))
sigA3.f$WBGeneID <- substr(row.names(sigA3.f), start = 6, stop = 19)

sigB2.f <- as.data.frame(subset(resB2.1, padj < bonf))
sigB2.f$WBGeneID <- substr(row.names(sigB2.f), start = 6, stop = 19)

sigB3.f <- as.data.frame(subset(resB3.1C, padj < bonf))
sigB3.f$WBGeneID <- substr(row.names(sigB3.f), start = 6, stop = 19)

femcolors <- c("#FA9FB5", "#F768A1", "#DD3497", "#AE017E")

venn.plot <- venn.diagram(
  x = list(sigA2.f$WBGeneID, sigA3.f$WBGeneID, sigB2.f$WBGeneID, sigB3.f$WBGeneID),
  category.names = c("A2" , "A3" , "B2", "B3"),
  filename = NULL,
  fill = alpha(femcolors, 0.5), col = femcolors, lwd = 2
)
plot_grid(venn.plot)
dev.copy2pdf(file="~/Desktop/2024.07.25_Female_SigOverlapByRep.pdf", useDingbats=FALSE, family="sans", width = 3, height = 3)


