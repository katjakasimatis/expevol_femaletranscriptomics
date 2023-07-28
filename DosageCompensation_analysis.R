library(DESeq2)


####Read in a count matrix####
gene <- read.table("Documents/UO.academics/Phillips_Lab/experiments/experimental_evolution/2018_FullExperiment/transcriptomics/20220617_Gene_table_counts.txt",
                   sep = "\t", header = FALSE, quote = "", row.names = 1)

coldata <- read.table("Dropbox/Experimental_Evolution/transcriptomics/name_headers.txt", sep = "\t", header = TRUE)
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


####Differential expression analysis####
#Censor: B3.anc (20_S36_L005Aligned.sortedByCoord.out.ba) and B3.B3 (21_S37_L005Aligned.sortedByCoord.out.ba)
censoredCts <- gene[, c(1:10, 13:39)]
censoredCol <- coldata[c(1:10, 13:39), ]
names(censoredCts) <- censoredCol$treatment

dds <- DESeqDataSetFromMatrix(countData = censoredCts, colData = censoredCol, design = ~ female2 + male2 + female2:male2)
dds <- DESeq(dds)

res1 <- results(dds, contrast = list("female2_evo_vs_anc"))

femcontrast <- as.data.frame(res1)
femcontrast$WBGeneID <- rownames(femcontrast)
femcontrast$WBGeneID <- substr(femcontrast$WBGeneID, start = 6, stop = 19)

bonf <- 0.05 / length(dds)



###Read in dosage compensation data from Jans et al. 2009####
dc <- read.table("PATH TO DOSAGE COMPENSATED GENE LIST", header = TRUE, sep = "\t")
notdc <- read.table("PATH TO NOT COMPENSATED GENE LIST", header = TRUE, sep = "\t")



####
allDC <- merge(x = femcontrast, y = dc, by = "WBGeneID")
allDC$compensated <- rep("yes", length(allDC$WBGeneID))

allNDC <- merge(x = femcontrast, y = notdc, by = "WBGeneID")
allNDC$compensated <- rep("no", length(allNDC$WBGeneID))

x <- rbind(allDC[, c(1:9, 12)], allNDC[, c(1:9, 11)])

plot(x$WT_XX, x$log2FoldChange, las = 1, xlab = "Normalized Expression in WT XX", ylab = "Expression Change in Evolved", bty = "l",
     pch = ifelse(x$compensated=="yes", 19, 17), cex = 1.4,
     col = ifelse(-log10(x$padj) > -log10(bonf), alpha("#DE77AE", 0.95), alpha("black", 0.75)))
legend("topleft", c("DC", "not DC"), pch = c(19, 17), bty = "n")
abline(h = -1.9)

mat.data <- c(length(subset(x, padj < bonf & compensated == "yes")$WBGeneID),
              length(subset(x, padj < bonf & compensated == "no")$WBGeneID),
              length(subset(x, padj > bonf & compensated == "yes")$WBGeneID),
              length(subset(x, padj > bonf & compensated == "no")$WBGeneID))
mat <- matrix(mat.data, nrow = 2, ncol = 2, byrow = TRUE)
chisq.test(mat)

