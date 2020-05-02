# docker run -it --rm --mount type=bind,source=/k11e/pvdisk/,target=/k11e/pvdisk/ wangjianshou/dada2:R_3.5.0 bash
# docker run -d --rm --mount type=bind,source=/k11e/pvdisk/,target=/k11e/pvdisk/ wangjianshou/dada2:R_3.5.0 Rscript /workdir/paired_dada2.R --runid runID_file --outdir path_of_outdir --sampleNameIndex 3
# docker run -d --rm --mount type=bind,source=/k11e/pvdisk/,target=/k11e/pvdisk/ wangjianshou/dada2:R_3.5.0 Rscript /workdir/paired_dada2.R --sample read1_path_file read2_path_file --outdir path_of_outdir  --sampleNameIndex 3
# --runid 可以通过runid file找到数据，一行一个runid，不过这个没用，因为这是16S数据，一开始写这个东西的时候不知道
# --sample read1_path_file read2_path_file，直接给出原始数据的路径，一行一个，--runid和--sample二者只能选一
# --outdir path_of_outdir 输出路径
# --sampleNameIndex 3 以.为分割，第几个字段作为样本名称

library(dada2)
library(DECIPHER)
library(phyloseq)
library(ggplot2)

#dada2Database <- "/k11e/pvdisk/fastbase/Users/wangjianshou/git/16s-vs-shallow-shotgou-sequencing/result/dada2/"
dada2Database <- "/workdir/"
ar <- commandArgs(trailingOnly=TRUE)
generateFqPath <- function(runid)
{
  dataDir <- "/k11e/pvdisk/bigbase/kbdata/sampledata/"
  samplePath <- substring(runid, 1, seq(nchar(runid)))
  R1 <- paste0(samplePath[length(samplePath)], ".R1.fastq.gz")
  R2 <- paste0(samplePath[length(samplePath)], ".R2.fastq.gz")
  R1path <- file.path(dataDir, paste0(samplePath, collapse='/'), R1)
  R2path <- file.path(dataDir, paste0(samplePath, collapse='/'), R2)
  return(c(R1path, R2path))
}

if(ar[1]=='--runid')
{
  runIDs <- readLines(ar[2])
  fqPath <- lapply(runIDs, generateFqPath)
  R1path <- unlist(lapply(fqPath, `[`, 1))
  R2path <- unlist(lapply(fqPath, `[`, 2))
} else if(ar[1]=='--sample')
{
  R1path <- readLines(ar[2])
  R2path <- readLines(ar[3])
}

outdir <- if('--outdir' %in% ar) ar[which(ar=='--outdir')+1] else getwd()
sampleNameIndex <- if('--sampleNameIndex' %in% ar) as.integer(ar[which(ar=='--sampleNameIndex')+1]) else 1
sampleNames <- sapply(strsplit(basename(R1path), split='[.]'), `[`, sampleNameIndex)
if(!dir.exists(file.path(outdir, 'figure'))) dir.create(file.path(outdir, 'figure'))
if(!dir.exists(file.path(outdir, 'FilterFq'))) dir.create(file.path(outdir, 'FilterFq'))

pdf(file.path(outdir, 'figure', 'R1_Quality.pdf'))
plotQualityProfile(R1path)
dev.off()

pdf(file.path(outdir, 'figure', 'R2_Quality.pdf'))
plotQualityProfile(R2path)
dev.off()

filterR1 <- file.path(outdir, 'FilterFq', paste0(sampleNames, ".R1.Filtered.fastq.gz"))
filterR2 <- file.path(outdir, 'FilterFq', paste0(sampleNames, ".R2.Filtered.fastq.gz"))
out <- filterAndTrim(R1path, filterR1, R2path, filterR2, trimRight=c(1, 40), minLen=100, maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE, compress=TRUE, multithread=TRUE)

pdf(file.path(outdir, 'figure', 'clean_R1_Quality.pdf'))
plotQualityProfile(filterR1)
dev.off()
pdf(file.path(outdir, 'figure', 'clean_R2_Quality.pdf'))
plotQualityProfile(filterR2)
dev.off()

if('--save-filter-image' %in% ar)
{
  save.image('filtered.RDATA')
  stopifnot(FALSE)
}
errR1 <- learnErrors(filterR1, multithread=TRUE)
errR2 <- learnErrors(filterR2, multithread=TRUE)
pdf(file.path(outdir, "figure", "R1_errModel.pdf"))
plotErrors(errR1, nominalQ=TRUE)
dev.off()
pdf(file.path(outdir, "figure", "R2_errModel.pdf"))
plotErrors(errR2, nominalQ=TRUE)
dev.off()

derepR1 <- derepFastq(filterR1, verbose=TRUE)
derepR2 <- derepFastq(filterR2, verbose=TRUE)
names(derepR1) <- sampleNames
names(derepR2) <- sampleNames

dadaR1 <- dada(derepR1, err=errR1, multithread=TRUE)
dadaR2 <- dada(derepR2, err=errR2, multithread=TRUE)

mergers <- mergePairs(dadaR1, derepR1, dadaR2, derepR2, verbose=TRUE)  # justConcatenate=TRUE used for non-overlapping reads
seqtab <- makeSequenceTable(mergers)
# dim(seqtab)
# table(nchar(getSequences(seqtab)))
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
write.table(seqtab.nochim, file.path(outdir, "seqtab.nochim"), sep='\t', quote=FALSE, row.names=T, col.names=T)

getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaR1, getN), sapply(dadaR2, getN), sapply(mergers, getN), rowSums(seqtab.nochim))

colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sampleNames
write.table(track, file.path(outdir, "track.info"), sep='\t', quote=FALSE, row.names=T, col.names=T)

taxa <- assignTaxonomy(seqtab.nochim, file.path(dada2Database, "silva_nr_v132_train_set.fa.gz"), multithread=TRUE, tryRC=TRUE)
taxa <- addSpecies(taxa, file.path(dada2Database, "silva_species_assignment_v132.fa.gz"), tryRC=TRUE)
write.table(taxa, file.path(outdir, "taxa_result_by_naive_Bayesian"), sep='\t', quote=FALSE, row.names=T, col.names=T)
#save.image("/dada2/taxa_1.RDATA")

dna <- DNAStringSet(getSequences(seqtab.nochim))
load(file.path(dada2Database, "SILVA_SSU_r132_March2018.RData"))
ids <- IdTaxa(dna, trainingSet, strand="both", processors=20, verbose=FALSE)
ranks <- c("domain", "phylum", "class", "order", "family", "genus", "species")

fun <- function(x)
{
  m <- match(ranks, x$rank)
  taxa <- x$taxon[m]
  taxa[startsWith(taxa, "unclassified_")] <- NA
  taxa
}

taxid <- t(sapply(ids, fun))
colnames(taxid) <- ranks
rownames(taxid) <- getSequences(seqtab.nochim)
write.table(taxid, file.path(outdir, "taxid.result"), sep='\t', quote=FALSE, row.names=T, col.names=T)
# head(taxa.print)
save.image(file.path(outdir, "result.RDATA"))

