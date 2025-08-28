---
title: "DADA2 Analysis Procedures"
author: "Anushka KC"
output: html_document
---

# load library
```{r}
library(dada2); packageVersion("dada2")
```
# set path
```{r}
path <- "fastq"
list.files(path)
```
#forward and reverse files
```{r}
fnFs <- sort(list.files(path, pattern="_1.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_2.fastq", full.names = TRUE))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
```


```{r}
plotQualityProfile(fnFs[1:2])
```

```{r}
plotQualityProfile(fnRs[1:2])
```

```{r}
# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names
```

```{r}
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(230, 120), trimLeft=c(0,20),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=FALSE) # On Windows set multithread=FALSE
```

```{r}
head(out)
```

```{r}
errF <- learnErrors(filtFs, multithread=TRUE)
```

```{r}
errR <- learnErrors(filtRs, multithread=TRUE)
```

```{r}
plotErrors(errF, nominalQ=TRUE)
```

```{r}
dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
```

```{r}
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)
```

```{r}
dadaFs[[1]]
```

```{r}
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
```

```{r}
# Inspect the merger data.frame from the first sample
head(mergers[[1]])
```

```{r}
seqtab <- makeSequenceTable(mergers)
dim(seqtab)
```

```{r}
# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))
```

```{r}
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
```

```{r}
dim(seqtab.nochim)
```

```{r}
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)
```

```{r}
taxa <- assignTaxonomy(seqtab.nochim, "../silva_nr99_v138.1_train_set.fa.gz", multithread=TRUE)
```

```{r}
taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)
```

```{r}
library(DECIPHER); packageVersion("DECIPHER")
```

```{r}
library(phangorn); packageVersion("phangorn")
```

```{r}
taxa.species <- addSpecies(taxa, "../silva_species_assignment_v138.1.fa.gz")
```

```{r}
ASVs.nochim = DNAStringSet(colnames(seqtab.nochim))
names(ASVs.nochim) = paste0("ASV", 1:ncol(seqtab.nochim))
```

```{r}
alignment = AlignSeqs(ASVs.nochim, anchor=NA, processors=30)
```

```{r}
phang.align <- phyDat(as(alignment, "matrix"), type="DNA")
dm <- dist.ml(phang.align)
treeNJ <- NJ(dm)

fit = pml(treeNJ, data=phang.align)
fitGTR <- update(fit, k=4, inv=0.2)
#fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE, rearrangement = "stochastic", control = pml.control(trace = 0))
save(fitGTR, file=file.path("fitGTR.RData"))
```

### Loading Metadata
```{r}
mdata <- read.csv("metadata.csv")
```

```{r}
mdata2 <- mdata[match(rownames(seqtab.nochim), mdata$Run), ]
rownames(mdata2) <- mdata2$Run

ASVs.nochim <- DNAStringSet(colnames(seqtab.nochim))
names(ASVs.nochim) <- paste0("ASV", 1:ncol(seqtab.nochim))

tmp.seqtab <- seqtab.nochim
colnames(tmp.seqtab) <- names(ASVs.nochim)
tmp.taxa <- taxa.species
rownames(tmp.taxa) <- names(ASVs.nochim)
```

### Creating Phyloseq Object
```{r}
library(phyloseq); packageVersion("phyloseq")
```

```{r}
ps <- phyloseq(
        otu_table(tmp.seqtab, taxa_are_rows=FALSE),
        sample_data(mdata2),
        tax_table(tmp.taxa),
        refseq(ASVs.nochim),
        phy_tree(fitGTR$tree))

ps
```

```{r}
save(ps, file=file.path("phyloseq_nochim_silva.RData"))
```

#### Saving ASV/OTU table
```{r}
OTU <- as(otu_table(ps), "matrix")
if(taxa_are_rows(ps)){OTU <- t(OTU)}
# Coerce to data.frame
OTUdf <- as.data.frame(OTU)
```

```{r}
write.table(OTUdf, "ASV_OTU_Table.tsv", row.names = TRUE)
```

#### Taxonomy Table
```{r}
taxonomy <- as(tax_table(ps), "matrix")
if(taxa_are_rows(ps)){taxonomy <- t(taxonomy)}
# Coerce to data.frame
taxonomydf <- as.data.frame(taxonomy)
```

```{r}
write.table(taxonomy, "Taxonomy.tsv", row.names = TRUE)
```

```{r}
write.csv(seqtab.nochim, "seqtab_nochim.csv", row.names = TRUE)
```

### Save SeqTab as RDS
```{r}
saveRDS(seqtab.nochim, "seqtab.nochim.rds")
```