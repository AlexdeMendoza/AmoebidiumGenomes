library(ggplot2)
library(dplyr)
library(data.table) #installPackages("data.table")Allows use of fread, and allows efficient manipulation of very large tables.
library(tidyverse)     ## data wrangling + ggplot2 install.packages("tidyverse")

CG_BS <- readRDS("CG_BS.rds")
GENE_GTF_PATH <- "./inputs/annotations/Amoebidium.chr.v3.no_TEs.gtf"
REPEAT_GTF_PATH <- "./inputs/annotations/Apar_genome.fasta.RM.bed.gtf.cleaned.final.kimura"
VIRAL_GTF_PATH <- "./inputs/annotations/GiantRepeats.gtf"
FASTA_INDEX_PATH = "./inputs/genomes/Apar_genome_pUC19_lambda.fasta.fai"
BLACKLISTED_REGION <- "./inputs/annotations/Blacklisted.regions.0based.bed"
BLACKLISTED_REGION_GRANGES <- import.bed(BLACKLISTED_REGION)
SAMPLES = "NA"
TREATMENTS = "NA"
get_mc_stats_from_granges <- function(BS, GR){
  strand(BS) <- "*"
  CG_coverage_positions <- getCoverage(BS[rowRanges(BS)$dinucleotide == "CG"], regions = GR, type = "Cov", what = "perBase")
  CG_position_per_gene <- matrix(unlist(sapply(CG_coverage_positions, function(x){
    return(ifelse(is.null(nrow(x)), 0, nrow(x)))
  })), length(GR))
  CG_MetRead <- getCoverage(BS[rowRanges(BS)$dinucleotide == "CG"], regions = GR, type = "M", what = "perRegionTotal")
  CG_TotRead <- getCoverage(BS[rowRanges(BS)$dinucleotide == "CG"], regions = GR, type = "Cov", what = "perRegionTotal")
  mCG_fraction <- matrix((CG_MetRead / CG_TotRead), length(GR))
  return(list(CG_position_per_gene, CG_MetRead, CG_TotRead, mCG_fraction))
}

index2slengths <- function(FASTA_INDEX_PATH){
  FASTA_INDEX <- fread(FASTA_INDEX_PATH)
  SLENGTHS <- FASTA_INDEX$V2
  names(SLENGTHS) <- FASTA_INDEX$V1
  return(SLENGTHS)
}


SLENGTHS <- index2slengths(FASTA_INDEX_PATH)
OUTPUT_DIR <- "./broad_features"

viral_annotation_se <- function(VIRAL_GTF_PATH, BS, SLENGTHS){
  VIRAL_GTF <- fread(VIRAL_GTF_PATH)
  gene_id <- sapply(strsplit(as.character(VIRAL_GTF$V9),";"), "[", 1)
  gene_id <- sapply(strsplit(as.character(gene_id),"\""), "[", 2)
  VIRAL_GTF$gene_id <- gene_id
  GR <- GRanges(seqnames = VIRAL_GTF$V1, 
                ranges = IRanges(start = as.numeric(VIRAL_GTF$V4), end = as.numeric(VIRAL_GTF$V5)), 
                strand = VIRAL_GTF$V6,
                seqlengths = SLENGTHS,
                gene_id = VIRAL_GTF$gene_id,
                family=VIRAL_GTF$V3)
  GR_STATS <- get_mc_stats_from_granges(BS, GR)
  SE <- SummarizedExperiment(assays=list(mCG_met = GR_STATS[[2]],
                                         mCG_cov = GR_STATS[[3]],
                                         mCG_fraction=GR_STATS[[4]]),
                             rowRanges=GR)
  rowRanges(SE)$CG_position_per_gene <- GR_STATS[[1]]
  return(SE)
}

repeat_annotation_summarizedexperiment <- function(BS, REPEAT_GTF_PATH, SLENGTHS){
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
  GR <- GRanges(seqnames = REPEAT_GTF$chr, 
                ranges = IRanges(start = as.numeric(REPEAT_GTF$start), end = as.numeric(REPEAT_GTF$end)), 
                strand = REPEAT_GTF$strand,
                seqlengths = SLENGTHS,
                gene_id = REPEAT_GTF$gene_id,
                family = REPEAT_GTF$family,
                class = REPEAT_GTF$class,
                group = REPEAT_GTF$group,
                kimura = REPEAT_GTF$Kimura)
  
  GR_STATS <- get_mc_stats_from_granges(BS, GR)
  SE <- SummarizedExperiment(assays=list(mCG_met = GR_STATS[[2]],
                                         mCG_cov = GR_STATS[[3]],
                                         mCG_fraction=GR_STATS[[4]]),
                             rowRanges=GR)
  rowRanges(SE)$CG_position_per_gene <- GR_STATS[[1]]
  return(SE)
}

gene_annotation_summarizedexperiment <- function(BS, GENE_GTF_PATH, SLENGTHS, UPSTREAM, DOWNSTREAM){
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
  PR <- promoters(txdb, upstream=UPSTREAM, downstream=DOWNSTREAM, columns="gene_id")
  PR_DT <- data.table(chr = as.character(seqnames(PR)), family = "promoter", start = start(PR), end = end(PR), strand = as.character(strand(PR)), gene_id = unlist(PR$gene_id), transcript_id = names(PR))
  GENE_DT <- data.table(chr = as.character(seqnames(GENES)), family = "gene", start = start(GENES), end = end(GENES), strand = as.character(strand(GENES)), gene_id = GENES$gene_id)
  longest_transcripts <- GENE_GTF[GENE_GTF$family=="transcript"] %>% group_by(gene_id) %>% summarise(max=max(end-start),                                                                                                   transcript_id = transcript_id[which.max(end-start)]) %>% select(transcript_id, gene_id)
  GENE_DT <- merge(GENE_DT, longest_transcripts, by="gene_id")
  GENE_GTF <- rbindlist(list(GENE_GTF, PR_DT, GENE_DT), use.names = TRUE)
  GENE_GTF <- GENE_GTF[GENE_GTF$transcript_id %in% longest_transcripts$transcript_id]
  GR <- GRanges(seqnames = GENE_GTF$chr, 
                ranges = IRanges(start = as.numeric(GENE_GTF$start), end = as.numeric(GENE_GTF$end)), 
                strand = GENE_GTF$strand,
                seqlengths = SLENGTHS,
                gene_id = GENE_GTF$gene_id,
                family = GENE_GTF$family)
  
  GR_STATS <- get_mc_stats_from_granges(BS, GR)
  SE <- SummarizedExperiment(assays=list(mCG_met = GR_STATS[[2]],
                                         mCG_cov = GR_STATS[[3]],
                                         mCG_fraction=GR_STATS[[4]]),
                             rowRanges=GR)
  rowRanges(SE)$CG_position_per_gene <- GR_STATS[[1]]
  return(SE)
}

middle_base<-3
CG_BS_CG <- CG_BS[rowRanges(CG_BS)$dinucleotide=="CG"]
BS_cgc_CG <- CG_BS_CG[substr(rowRanges(CG_BS_CG)$context, middle_base, middle_base+2) == "cgc" |
                        substr(rowRanges(CG_BS_CG)$context, middle_base-1, middle_base+1) == "gcg"]
BS_non_cgc_CG <- CG_BS_CG[!(substr(rowRanges(CG_BS_CG)$context, middle_base, middle_base+2) == "cgc" |
                              substr(rowRanges(CG_BS_CG)$context, middle_base-1, middle_base+1) == "gcg")]

broad_features <- function(BS, use_ch = FALSE, title){
  #Load in Summarized Experiment for viral regions:
  function_viral_se <- viral_annotation_se(VIRAL_GTF_PATH, BS, SLENGTHS)
  
  #Load in Summarized Experiment for repeats regions:
  function_repeat_se <- repeat_annotation_summarizedexperiment(BS, REPEAT_GTF_PATH, SLENGTHS)
  
  #Subset for those not in blacklisted regions:
  function_repeat_se <- subsetByOverlaps(function_repeat_se, BLACKLISTED_REGION_GRANGES, invert = TRUE, ignore.strand=TRUE)
  #Subset for those not in viral regions:
  function_repeat_se <- subsetByOverlaps(function_repeat_se, function_viral_se, invert = TRUE, ignore.strand=TRUE)
  # Exclude >500bp
  function_repeat_se <- function_repeat_se[width(function_repeat_se)>500]
  
  #Load in Summarized Experiment for genes:
  function_gene_se <- gene_annotation_summarizedexperiment(BS, GENE_GTF_PATH, SLENGTHS, UPSTREAM = 2000, DOWNSTREAM = 400)
  #Subset for those not in blacklisted regions:
  function_gene_se <- subsetByOverlaps(function_gene_se, BLACKLISTED_REGION_GRANGES, invert = TRUE, ignore.strand=TRUE)
  #Subset for those not in viral regions:
  function_gene_se <- subsetByOverlaps(function_gene_se, function_viral_se, invert = TRUE, ignore.strand=TRUE)
  
  #Seperate into gene SE, and promoter SE
  function_only_genes <- function_gene_se[rowData(function_gene_se)$family=="gene"]
  function_only_promoters <- function_gene_se[rowData(function_gene_se)$family=="promoter"]
  
  #Make list of SEs to be looped over
  SE_LIST <- list(function_only_genes, function_only_promoters, function_repeat_se, function_viral_se)
  SE_names <- c("Genes", "Promoters", "Repeats", "Viral Insertions")
  BROAD_MC_DATA <- lapply(seq_along(SE_LIST), function(x){
    dt <- data.table(gene_id = rowData(SE_LIST[[x]])$gene_id,
                     mcg = as.vector(assays(SE_LIST[[x]])$mCG_fraction),
                     sample = rep(SAMPLES, each = nrow(SE_LIST[[x]])),
                     condition = rep(TREATMENTS, each = nrow(SE_LIST[[x]])),
                     Feature = SE_names[x])
    return(dt)
  })
  BROAD_MC_DATA <- rbindlist(BROAD_MC_DATA)
  
  Global_mCG <- sum(assays(BS[rowRanges(BS)$source=="Organism"])$M)/sum(assays(BS[rowRanges(BS)$source=="Organism"])$Cov)
  output <- ggplot(data = BROAD_MC_DATA, mapping = aes (x = Feature, y = (mcg*100))) +
    geom_boxplot(outlier.shape = NA) + ylab("m (%)") + xlab("Genomic feature") +
    ggtitle(title) +
    geom_hline(yintercept=Global_mCG*100, linetype="dashed", color = "red") + 
    theme_bw()
  ggsave(paste("MainFeatures_m", title, ".pdf", sep=""), plot = output, path = OUTPUT_DIR, width = 10, height = 10, units = "cm")
}

broad_features(BS_cgc_CG, title = "CGC")
broad_features(BS_non_cgc_CG, title = "(non)CGC")