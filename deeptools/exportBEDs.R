library(data.table)
library(GenomicFeatures)
library(dplyr)
library(SummarizedExperiment)
library(rtracklayer)
setwd("/data/SBCS-ademendoza/02-lukesarre/amoebidium_paper_trimmed")

GENE_GTF_PATH <- "./inputs/annotations/Amoebidium.chr.v3.no_TEs.gtf"
REPEAT_GTF_PATH <- "./inputs/annotations/Apar_genome.fasta.RM.bed.gtf.cleaned.final.kimura"
VIRAL_GTF_PATH <- "./inputs/annotations/GiantRepeats.gtf"
EXPRESSION_DATA_PATH = "./inputs/expression/A_par_static.cntTable.cleaned"
FASTA_INDEX_PATH = "./inputs/genomes/Apar_genome_pUC19_lambda.fasta.fai"
BLACKLISTED_REGION = "./inputs/annotations/Blacklisted.regions.0based.bed"
load_expression_data <- function(EXPRESSION_DATA_PATH){
  EXPRESSION_DATA <- fread(EXPRESSION_DATA_PATH[[1]])
  counts <- c(EXPRESSION_DATA[[2]])
  if (length(SAMPLES) > 1){
    for (SAMPLE_NUMBER in c(2:length(SAMPLES))) {
      EXPRESSION_DATA_TEMP <- fread(EXPRESSION_DATA_PATH[[SAMPLE_NUMBER]])
      counts <- c(counts, EXPRESSION_DATA_TEMP[[2]])
    }}
  MATRIX <- matrix(counts, nrow(EXPRESSION_DATA))
  OUTPUT <- as.data.frame(MATRIX)
  OUTPUT$gene_id <- EXPRESSION_DATA[[1]]
  return(OUTPUT)
}

viral_annotation_granges <- function(VIRAL_GTF_PATH, SLENGTHS){
  VIRAL_GTF <- fread(VIRAL_GTF_PATH)
  gene_id <- sapply(strsplit(as.character(VIRAL_GTF$V9),";"), "[", 1)
  gene_id <- sapply(strsplit(as.character(gene_id),"\""), "[", 2)
  VIRAL_GTF$gene_id <- gene_id
  return(GRanges(seqnames = VIRAL_GTF$V1, 
                ranges = IRanges(start = as.numeric(VIRAL_GTF$V4), end = as.numeric(VIRAL_GTF$V5)), 
                strand = VIRAL_GTF$V6,
                seqlengths = SLENGTHS,
                gene_id = VIRAL_GTF$gene_id,
                family=VIRAL_GTF$V3))
}

repeat_annotation_granges <- function(REPEAT_GTF_PATH, SLENGTHS){
  #Read in the GTF
  REPEAT_GTF <- fread(REPEAT_GTF_PATH)
  #Parse the 9th column, to get repeat information
  gene_id <- sapply(strsplit(as.character(sapply(strsplit(as.character(REPEAT_GTF$V9),";"), "[", 2)),"\""), "[", 2)
  group_id <- sapply(strsplit(sapply(strsplit(as.character(REPEAT_GTF$V9),"gene_id \""), "[", 2),"\";"), "[", 1)
  class <- sapply(strsplit(sapply(strsplit(as.character(REPEAT_GTF$V9),"class_id \""), "[", 2),"\";"), "[", 1)
  REPEAT_GTF$gene_id <- gene_id
  REPEAT_GTF$group <- group_id
  REPEAT_GTF$class <- class
  
  REPEAT_GTF <- dplyr::select(REPEAT_GTF, V1, V3, V4, V5, V7, V8, gene_id, group, class)
  REPEAT_GTF <- dplyr::rename(REPEAT_GTF, "chr" = "V1", "family" = "V3", "start" = "V4", "end" = "V5", "strand" = "V7", "Kimura" = "V8")
  
  #Remove duplicated repeats
  REPEAT_GTF <- REPEAT_GTF[order(REPEAT_GTF[,c('chr','start','end','strand','family','Kimura')]),]
  REPEAT_GTF <- REPEAT_GTF[!duplicated(REPEAT_GTF[,c('chr','start','end','family')]),]
  
  #Create a GRanges object describing the repeats. This is one part of the SE object.
  return(GR <- GRanges(seqnames = REPEAT_GTF$chr, 
                ranges = IRanges(start = as.numeric(REPEAT_GTF$start), end = as.numeric(REPEAT_GTF$end)), 
                strand = REPEAT_GTF$strand,
                seqlengths = SLENGTHS,
                gene_id = REPEAT_GTF$gene_id,
                family = REPEAT_GTF$family,
                class = REPEAT_GTF$class,
                group = REPEAT_GTF$group,
                kimura = REPEAT_GTF$Kimura))
}

gene_annotation_summarizedexperiment <- function(GENE_GTF_PATH, EXPRESSION_DATA, SLENGTHS){
  GENE_GTF <- fread(GENE_GTF_PATH)
  GENE_GTF <- dplyr::select(GENE_GTF, V1, V3, V4, V5, V7, V9)
  gene_id <- sapply(strsplit(as.character(sapply(strsplit(as.character(GENE_GTF$V9),";"), "[", 2)),"\""), "[", 2)
  transcript_id <- sapply(strsplit(as.character(sapply(strsplit(as.character(GENE_GTF$V9),";"), "[", 1)),"\""), "[", 2)
  GENE_GTF$gene_id <- gene_id
  GENE_GTF$transcript_id <- transcript_id
  GENE_GTF <- GENE_GTF[,-c("V9")]
  GENE_GTF <- dplyr::rename(GENE_GTF, chr = V1, family = V3, start = V4, end = V5, strand = V7)
  txdb <- makeTxDbFromGFF(file = GENE_GTF_PATH, format = "gtf", taxonomyId=4881)
  GENES <- GenomicFeatures::genes(txdb)
  GENE_DT <- data.table(chr = as.character(seqnames(GENES)), family = "gene", start = start(GENES), end = end(GENES), strand = as.character(strand(GENES)), gene_id = GENES$gene_id)
  longest_transcripts <- GENE_GTF[GENE_GTF$family=="transcript"] %>% group_by(gene_id) %>% summarise(max=max(end-start),
                                                                                                       transcript_id = transcript_id[which.max(end-start)]) %>% dplyr::select(transcript_id, gene_id)
  GENE_DT <- merge(GENE_DT, longest_transcripts, by="gene_id")
  GENE_GTF <- rbindlist(list(GENE_GTF, GENE_DT), use.names = TRUE)
  GENE_GTF <- GENE_GTF[GENE_GTF$transcript_id %in% longest_transcripts$transcript_id]
    
  GENE_GTF <- (merge(GENE_GTF, EXPRESSION_DATA, by = 'gene_id')) 
  counts <- as.matrix(GENE_GTF[,(ncol(GENE_GTF)-length(SAMPLES)+1):ncol(GENE_GTF)])
  colnames(counts) <- NULL
  GR <- GRanges(seqnames = GENE_GTF$chr, 
                ranges = IRanges(start = as.numeric(GENE_GTF$start), end = as.numeric(GENE_GTF$end)), 
                strand = GENE_GTF$strand,
                seqlengths = SLENGTHS,
                gene_id = GENE_GTF$gene_id,
                family = GENE_GTF$family)
  SE <- SummarizedExperiment(assays=list(counts=counts),
                             rowRanges=GR)
  return(SE)
}

index2slengths <- function(FASTA_INDEX_PATH){
  FASTA_INDEX <- fread(FASTA_INDEX_PATH)
  SLENGTHS <- FASTA_INDEX$V2
  names(SLENGTHS) <- FASTA_INDEX$V1
  return(SLENGTHS)
}

SAMPLES="NA"
SLENGTHS <- index2slengths(FASTA_INDEX_PATH)
EXPRESSION_DATA <- load_expression_data(EXPRESSION_DATA_PATH)
BLACKLISTED_REGION_GRANGES <- import.bed(BLACKLISTED_REGION)
#Load GRanges objects to retrieve BED files from:
CURATED_REPEATS <- repeat_annotation_granges(REPEAT_GTF_PATH, SLENGTHS)
gene_SE <- gene_annotation_summarizedexperiment(GENE_GTF_PATH, EXPRESSION_DATA, SLENGTHS)
CURATED_GENES <- gene_SE[rowRanges(gene_SE)$family == "gene"]
viral_GR <- viral_annotation_granges(VIRAL_GTF_PATH, SLENGTHS)

#First we export a bed file of repeats that are not in viruses or in blacklisted regions
# CURATED_REPEATS <- subsetByOverlaps(CURATED_REPEATS, viral_GR, invert=TRUE, ignore.strand=TRUE)
CURATED_REPEATS <- subsetByOverlaps(CURATED_REPEATS, BLACKLISTED_REGION_GRANGES, invert=TRUE, ignore.strand=TRUE)
CURATED_REPEATS <- CURATED_REPEATS[width(CURATED_REPEATS)>500]
CURATED_REPEATS <- CURATED_REPEATS[CURATED_REPEATS$family %in% c("LTR", "DNA", "LINE")]

for (repeat_family in unique(CURATED_REPEATS$family)){
  export.bed(CURATED_REPEATS[CURATED_REPEATS$family==repeat_family], 
             file.path(getwd(), "deeptools/BEDs", paste(repeat_family, ".bed", sep="")))
}

#Export gene bed:
# CURATED_GENES <- subsetByOverlaps(CURATED_GENES, viral_GR, invert=TRUE)
CURATED_GENES <- subsetByOverlaps(gene_SE, BLACKLISTED_REGION_GRANGES, invert=TRUE)
export.bed(granges(CURATED_GENES), file.path(getwd(), "deeptools/BEDs", "genes.bed"))

#Split gene SE object by decile and export as bed file for each decile
unexpressed_genes <- gene_SE[assays(gene_SE)$counts == 0]
export.bed(granges(unexpressed_genes), 
           file.path(getwd(), "deeptools/BEDs", "unexpressed_genes.bed"))

expressed_genes <- gene_SE[assays(gene_SE)$counts > 0]
rowRanges(expressed_genes)$decile <- ntile(assays(expressed_genes)$counts, 10)

for (decile in unique(rowRanges(expressed_genes)$decile)) {
  export.bed(granges(expressed_genes[rowRanges(expressed_genes)$decile==decile]), 
             file.path(getwd(), "deeptools/BEDs", paste("gene_expression_decile_", decile, ".bed", sep="")))
}
