
####import data####
library(tximport)
quant_files <- file.path("Project_borealis_sexual_antagonism/data/supertranscriptome/counts_salmon/quants", 
                         list.files("Project_borealis_sexual_antagonism/data/supertranscriptome/counts_salmon/quants"), "quant.sf")
samples <- c("BJE3896_dad_quant",
             "BJE3897_mom_quant",
             "BJE3929_boy_quant",
             "BJE4009_girl_quant",
             "BJE4017_boy_quant", 
             "BJE4039_boy_quant",
             "BJE4072_girl_quant",
             "BJE4082_girl_quant")

names(quant_files) <- samples

txi <- tximport(quant_files, type="salmon", txOut=TRUE,
                countsFromAbundance="scaledTPM")
cts <- txi$counts
cts <- cts[rowSums(cts) > 0,]


#load in the condition information
samps <- read.csv(file.path("Project_borealis_sexual_antagonism/data/supertranscriptome/DTU_analysis", "sample_condition.csv"))
samps$condition <- factor(samps$condition)
table(samps$condition)

#transcript to gene mapping
library(GenomicFeatures)
gtf <- "Project_borealis_sexual_antagonism/data/supertranscriptome/supertranscriptome_fasta/borealis_superTrans.gtf"
txdb.filename <- "borealis_superTrans_gtf.sqlite"
txdb <- makeTxDbFromGFF(gtf)
saveDb(txdb, txdb.filename)
txdb <- loadDb(txdb.filename)
txdf <- select(txdb, keys(txdb, "GENEID"), "TXNAME", "GENEID")
tab <- table(txdf$GENEID)
txdf$ntx <- tab[match(txdf$GENEID, names(tab))]

#build a DRIMSeq object
cts[1:3,1:3]
range(colSums(cts)/1e6)
all(rownames(cts) %in% txdf$TXNAME)
txdf <- txdf[match(rownames(cts),txdf$TXNAME),]
all(rownames(cts) == txdf$TXNAME)

counts <- data.frame(gene_id=txdf$GENEID,
                     feature_id=txdf$TXNAME,
                     cts)
library(DRIMSeq)
d <- dmDSdata(counts=counts, samples=samps)
d

#build a DEXSeq 
library(DEXSeq)
sample.data <- DRIMSeq::samples(d)
count.data <- round(as.matrix(counts(d)[,-c(1:2)]))
dxd <- DEXSeqDataSet(countData=count.data,
                     sampleData=sample.data,
                     design=~sample + exon + condition:exon,
                     featureID=counts(d)$feature_id,
                     groupID=counts(d)$gene_id)

system.time({
  dxd <- estimateSizeFactors(dxd)
  dxd <- estimateDispersions(dxd, quiet=TRUE)
  dxd <- testForDEU(dxd, reducedModel=~sample + exon)
})

dxr <- DEXSeqResults(dxd, independentFiltering=FALSE)
qval <- perGeneQValue(dxr)
dxr.g <- data.frame(gene=names(qval),qval)


columns <- c("featureID","groupID","pvalue")
dxr <- as.data.frame(dxr[,columns])
head(dxr)


#######StageR
ata(dex_tables)

library(stageR)
strp <- function(x) substr(x,1,15)
pConfirmation <- matrix(dxr$pvalue,ncol=1)
dimnames(pConfirmation) <- list(strp(dxr$featureID),"transcript")
pScreen <- qval
names(pScreen) <- strp(names(pScreen))
tx2gene <- as.data.frame(dxr[,c("featureID", "groupID")])
for (i in 1:2) tx2gene[,i] <- strp(tx2gene[,i])

stageRObj <- stageRTx(pScreen=pScreen, pConfirmation=pConfirmation,
                      pScreenAdjusted=TRUE, tx2gene=tx2gene)
stageRObj <- stageWiseAdjustment(stageRObj, method="dtu", alpha=0.05)
suppressWarnings({
  dex.padj <- getAdjustedPValues(stageRObj, order=FALSE,
                                 onlySignificantGenes=TRUE)
})
head(dex.padj)