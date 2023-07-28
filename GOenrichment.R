library(DESeq2)
library(goseq)
library(TxDb.Celegans.UCSC.ce11.ensGene)
library(GO.db)



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



####GO enrichment analysis####
assayed.genes <- femcontrast$WBGeneID

txdb <- TxDb.Celegans.UCSC.ce11.ensGene
txsByGene <- transcriptsBy(txdb, "gene")
lengthData <-median(width(txsByGene))
lengthData <- as.array(lengthData)
names(lengthData) <- substr(names(lengthData), start = 1, stop = 14)

tokeep <- assayed.genes[(assayed.genes %in% names(lengthData))]

newlengths <- lengthData[(names(lengthData) %in% tokeep)]

de.genes <- subset(femcontrast, padj < bonf)$WBGeneID
gene.vector <- as.integer(tokeep %in% de.genes)
names(gene.vector) <- tokeep

pwf <- nullp(gene.vector, genome = "ce11", bias.data = newlengths)

GO.wall <- goseq(pwf, genome = "ce6", "ensGene")
head(GO.wall)

enriched.GO <- GO.wall$category[p.adjust(GO.wall$over_represented_pvalue, method = "BH") < 0.05]

for(go in enriched.GO[1:10]){
  print(GOTERM[[go]])
  cat("--------------------------------------\n")
}
