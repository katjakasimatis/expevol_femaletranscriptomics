library(DESeq2)
library(pheatmap)
library(ggplot2)
library(scales)


####Read in a count matrix####
gene <- read.table("PATH TO/20220617_Gene_table_counts.txt", sep = "\t", header = FALSE, quote = "", row.names = 1)

coldata <- read.table("PATH TO/name_headers.txt", sep = "\t", header = TRUE)
coldata$treatment <- interaction(coldata$female, coldata$male, sep = ".")

coldata$female <- as.factor(coldata$female)
coldata$male <- as.factor(coldata$male)

z <- c();
for (i in 1:length(coldata$female)) {
  if (coldata$female[i] == coldata$male[i]){
    out <- coldata$female[i]
  }
  else {
    if (coldata$female[i] == "anc"){
      out <- coldata$male[i]
    }
    else {
      out <- coldata$female[i]
    }
  }
  z <- rbind(z, as.character(out))
}
coldata$EErep <- z
coldata$EErep <- as.factor(coldata$EErep)

f2 <- c();
for (f in 1:length(coldata$female)) {
  if (coldata$female[f] == "anc") {
    out <- "anc"
  }
  else {
    out <- "evo"
  }
  f2 <- rbind(f2, out)
}
coldata$female2 <- f2
coldata$female2 <- as.factor(coldata$female2)

m2 <- c();
for (m in 1:length(coldata$male)) {
  if (coldata$male[m] == "anc") {
    out <- "anc"
  }
  else {
    out <- "evo"
  }
  m2 <- rbind(m2, out)
}
coldata$male2 <- m2
coldata$male2 <- as.factor(coldata$male2)

names(gene) <- coldata$treatment
gene <- gene[1:(length(gene$anc.A2)-5), ]



####Analyze A2 data####
A2coldata <- subset(coldata, EErep == "A2" | EErep == "anc")
toMatch   <- c("anc.anc", "anc.A2", "A2.anc", "A2.A2")
A2cts     <- gene[, grep(paste(toMatch,collapse="|"), names(gene))]
names(A2cts) <- A2coldata$treatment

#Differential expression analysis
ddsA2 <- DESeqDataSetFromMatrix(countData = A2cts, colData = A2coldata, design = ~ female + male + female:male)
ddsA2 <- DESeq(ddsA2)
resultsNames(ddsA2)

#results: female only (anc vs A2)
resA2.1 <- results(ddsA2, contrast = list("female_anc_vs_A2"))
#results: male only (anc vs A2)
resA2.2 <- results(ddsA2, contrast = list("male_anc_vs_A2"))
#results: interaction of female anc-evo and male anc-evo
resA2.3 <- results(ddsA2, contrast = list("femaleanc.maleanc"))

#variance stabilizing transformation
vsdA2 <- vst(ddsA2, blind=FALSE)

#PCA of full model
plotPCA(vsdA2, intgroup = c("treatment")) + theme_bw() + guides(color = guide_legend("Female x Male"))



####Analyze A3 data####
A3coldata <- subset(coldata, EErep == "A3" | EErep == "anc")
toMatch   <- c("anc.anc", "anc.A3", "A3.anc", "A3.A3")
A3cts     <- gene[, grep(paste(toMatch,collapse="|"), names(gene))]
names(A3cts) <- A3coldata$treatment

#Differential expression analysis
ddsA3 <- DESeqDataSetFromMatrix(countData = A3cts, colData = A3coldata, design = ~ female + male + female:male)
ddsA3 <- DESeq(ddsA3)
resultsNames(ddsA3)

#results: female only (anc vs A3)
resA3.1 <- results(ddsA3, contrast = list("female_anc_vs_A3"))
#results: male only (anc vs A3)
resA3.2 <- results(ddsA3, contrast = list("male_anc_vs_A3"))
#results: interaction of female anc-evo and male anc-evo
resA3.3 <- results(ddsA3, contrast = list("femaleanc.maleanc"))

#variance stabilizing transformation
vsdA3 <- vst(ddsA3, blind=FALSE)

#PCA of full model
plotPCA(vsdA3, intgroup = c("treatment")) + theme_bw() + guides(color = guide_legend("Female x Male"))



####Analyze B2 data####
B2coldata <- subset(coldata, EErep == "B2" | EErep == "anc")
toMatch   <- c("anc.anc", "anc.B2", "B2.anc", "B2.B2")
B2cts     <- gene[, grep(paste(toMatch,collapse="|"), names(gene))]
names(B2cts) <- B2coldata$treatment

#Differential expression analysis
ddsB2 <- DESeqDataSetFromMatrix(countData = B2cts, colData = B2coldata, design = ~ female + male + female:male)
ddsB2 <- DESeq(ddsB2)
resultsNames(ddsB2)

#results: female only (anc vs B2)
resB2.1 <- results(ddsB2, contrast = list("female_B2_vs_anc"))
#results: male only (anc vs A3)
resB2.2 <- results(ddsB2, contrast = list("male_B2_vs_anc"))
#results: interaction of female anc-evo and male anc-evo
resB2.3 <- results(ddsB2, contrast = list("femaleB2.maleB2"))

#variance stabilizing transformation
vsdB2 <- vst(ddsB2, blind=FALSE)

#PCA of full model
plotPCA(vsdB2, intgroup = c("treatment")) + theme_bw() + guides(color = guide_legend("Female x Male"))



####Analyze B3 data####
B3coldata <- subset(coldata, EErep == "B3" | EErep == "anc")
toMatch   <- c("anc.anc", "anc.B3", "B3.anc", "B3.B3")
B3cts     <- gene[, grep(paste(toMatch,collapse="|"), names(gene))]

#Censor: B3.anc (20_S36_L005Aligned.sortedByCoord.out.ba) and B3.B3 (21_S37_L005Aligned.sortedByCoord.out.ba)
B3cts        <- B3cts[, c(1:3, 6:12)]
B3coldata    <- B3coldata[c(1:3, 6:12),]
names(B3cts) <- B3coldata$treatment

#Differential expression analysis
ddsB3 <- DESeqDataSetFromMatrix(countData = B3cts, colData = B3coldata, design = ~ female + male + female:male)
ddsB3 <- DESeq(ddsB3)
resultsNames(ddsB3)

#results: female only (anc vs B3)
resB3.1 <- results(ddsB3, contrast = list("female_B3_vs_anc"))
#results: male only (anc vs B3)
resB3.2 <- results(ddsB3, contrast = list("male_B3_vs_anc"))
#results: interaction of female anc-evo and male anc-evo
resB3.4 <- results(ddsB3, contrast = list("femaleB3.maleB3"))

#variance stabilizing transformation
vsdB3 <- vst(ddsB3, blind=FALSE)

#PCA of full model
plotPCA(vsdB3, intgroup = c("treatment")) + theme_bw() + guides(color = guide_legend("Female x Male"))



####Analyze full data as reduced ancestor vs evolved####
#Censor: B3.anc (20_S36_L005Aligned.sortedByCoord.out.ba) and B3.B3 (21_S37_L005Aligned.sortedByCoord.out.ba)
censoredCts <- gene[, c(1:10, 13:39)]
censoredCol <- coldata[c(1:10, 13:39), ]
names(censoredCts) <- censoredCol$treatment
#Differential expression analysis
dds <- DESeqDataSetFromMatrix(countData = censoredCts, colData = censoredCol, design = ~ female2 + male2 + female2:male2)
dds <- DESeq(dds)

#results: female only (anc vs evo)
res1 <- results(dds, contrast = list("female2_evo_vs_anc"))
#results: male only (anc vs evo)
res2 <- results(dds, contrast = list("male2_evo_vs_anc"))
#results: interaction of female anc-evo and male anc-evo
res3 <- results(dds, contrast = list("female2evo.male2evo"))

#variance stabilizing transformation
vsd <- vst(dds, blind=FALSE)

#PCA of full model
plotPCA(vsd, intgroup = c("female2", "male2")) + theme_bw() + geom_point(size = 4) + guides(color = guide_legend("Female x Male"))


#volcano plot
par(mfrow = c(1, 3))
bonf <- 0.05 / length(dds)
r <- c("Female Contrast", "Male Contrast", "Interaction")
for (i in 1:3) {
  plot(get(paste0("res", i))$log2FoldChange, -log10(get(paste0("res", i))$padj), las = 1, pch = 19, ylim = c(0, 50),
       xlab = "log2FoldChange", ylab = "-log10Padj", main = r[i],
       col = ifelse(-log10(get(paste0("res", i))$padj) > -log10(bonf), alpha("#DE77AE", 0.75), alpha("black", 0.4)))
  abline(h = -log10(bonf))
  abline(v = -2)
  abline(v = 2)
}

#visualize differential expression across samples with a heatmap
sig <- subset(res1, padj < bonf)
df <- as.data.frame(colData(dds)[ ,c("male2", "female2")])
mycols <- colorRampPalette(brewer.pal(9, "YlGnBu"))(20)
pheatmap(assay(vsd)[row.names(sig), ], cluster_rows = FALSE, cluster_cols = TRUE, show_colnames = TRUE, col = mycols, annotation_col = df)


# Mean and variance relationship
norm.counts <- counts(dds, normalized = TRUE)
mean.counts <- rowMeans(norm.counts)
variance.counts <- apply(norm.counts, 1, var)

mean.var.col <- densCols(x = log2(mean.counts), y = log2(variance.counts))
plot(x = log2(mean.counts), y = log2(variance.counts), pch = 20, cex = 0.8, col = mean.var.col, main = "Mean-variance relationship",
     xlab = "Mean log2(normalized counts) per gene", ylab = "Variance of log2(normalized counts)",
     panel.first = grid(), las = 1)
abline(a = 0, b = 1, col = "darkred")

dds.disp <- estimateDispersions(dds)
plotDispEsts(dds.disp)
