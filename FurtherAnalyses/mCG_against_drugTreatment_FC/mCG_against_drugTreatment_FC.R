library(ggplot2)
library(dplyr)
library(data.table) #installPackages("data.table")Allows use of fread, and allows efficient manipulation of very large tables.
library(DESeq2)

GENE_GTF_PATH <- "Amoebidium.chr.v3.no_TEs.gff3" # This is Amoebidium.chr.v3.no_TEs.gff3, excluding TEs present in DifferentialExpression/TE.ORFs.ids

#These expression files contain the gene names and counts from the relevant columns in AmoebidiumTreatment.TElocal.cntTable.gz
EXPRESSION_DATA_PATH = c("DMSO1.bam.cntTable.cleaned",
                         "DMSO2.bam.cntTable.cleaned",
                         "DMSO3.bam.cntTable.cleaned",
                         "Aza1.bam.cntTable.cleaned",
                         "Aza2.bam.cntTable.cleaned",
                         "Aza3.bam.cntTable.cleaned")
FASTA_INDEX_PATH = "Aapp_genome_pUC19_lambda.fasta.fai"

SAMPLES <- c("DMSO1", "DMSO2", "DMSO3", "Aza1", "Aza2", "Aza3")
TREATMENTS <- c(rep("DMSO", 3), rep("Aza", 3))

#### Define functions ####

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

load_gene_expression_data <- function(GENE_GTF_PATH, EXPRESSION_DATA_PATH){
  EXPRESSION_DATA <- load_expression_data(EXPRESSION_DATA_PATH)
  GENE_GTF <- fread(GENE_GTF_PATH)
  GENE_GTF <- dplyr::select(GENE_GTF, V1, V3, V4, V5, V7, V9)
  gene_id <- sapply(strsplit(as.character(sapply(strsplit(as.character(GENE_GTF$V9),";"), "[", 2)),"\""), "[", 2)
  transcript_id <- sapply(strsplit(as.character(sapply(strsplit(as.character(GENE_GTF$V9),";"), "[", 1)),"\""), "[", 2)
  GENE_GTF$gene_id <- gene_id
  GENE_GTF$transcript_id <- transcript_id
  GENE_GTF <- GENE_GTF[,-c("V9")]
  return(EXPRESSION_DATA[EXPRESSION_DATA$gene_id %in% GENE_GTF$gene_id,])
}

OUTPUT_DIR <- "mC_against_drugTrt_FC"
index2slengths <- function(FASTA_INDEX_PATH){
  FASTA_INDEX <- fread(FASTA_INDEX_PATH)
  SLENGTHS <- FASTA_INDEX$V2
  names(SLENGTHS) <- FASTA_INDEX$V1
  return(SLENGTHS)
}

SLENGTHS <- index2slengths(FASTA_INDEX_PATH)

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

#### Load in drug-treatment expression data, and run DEseq2 ####

geneExpressionData <- load_gene_expression_data(GENE_GTF_PATH, EXPRESSION_DATA_PATH)

EXPRESSION_MATRIX <- as.matrix(geneExpressionData[,1:6])
dimnames(EXPRESSION_MATRIX) <- list(geneExpressionData$gene_id,
                                    SAMPLES)

coldata <- DataFrame(Treatment=TREATMENTS,
                     row.names=SAMPLES)
dds <- DESeqDataSetFromMatrix(countData = EXPRESSION_MATRIX,
                              colData = coldata,
                              design = ~ Treatment)
TREATMENT_VECTOR <- c("DMSO", "Aza")
dds$Treatment <- factor(dds$Treatment, levels = TREATMENT_VECTOR)
dds <- estimateSizeFactors(dds)
rowData(dds)$control_expr <- rowMeans(counts(dds, normalize = TRUE)[,dds$Treatment == TREATMENT_VECTOR[[1]]])
rowData(dds)$condition_expr <- rowMeans(counts(dds, normalize = TRUE)[,dds$Treatment == TREATMENT_VECTOR[[2]]])
dds <- DESeq(dds)
res <- results(dds)
fdrlevel.de <- 0.05
res$control_expr <- rowData(dds)$control_expr
res$condition_expr <- rowData(dds)$condition_expr
dds_good <- dds[!is.na(res$padj)]
res_good <- res %>% na.omit()
deseq_res <- res
deseq_df <- data.frame(baseMean = deseq_res$baseMean,
                       log2FoldChange = deseq_res$log2FoldChange,
                       lfcSE = deseq_res$lfcSE,
                       stat = deseq_res$stat,
                       pvalue = deseq_res$pvalue,
                       padj = deseq_res$padj,
                       control_exp = deseq_res$control_expr,
                       condition_expr = deseq_res$condition_expr)
deseq_df$gene_id <- row.names(deseq_df)

#### Compare FC from expression data and compare to methylation in untreated condition ####

compare_methylation_from_emseq_vs_expression_change_in_drug_trt <- function(BS, figure_title){
  GENE_SE <- gene_annotation_summarizedexperiment(BS, GENE_GTF_PATH, SLENGTHS, UPSTREAM = 2000, DOWNSTREAM = 400)
  subject_name <- "promoter"
  subject <- GENE_SE[rowData(GENE_SE)$family==subject_name]
  subject <- subject[!is.na(assays(subject)$mCG_fraction) &
                       !is.nan(assays(subject)$mCG_fraction)]
  mcols(subject)$decile <- ntile(assays(subject)$mCG_fraction, 10)
  
  untreated_deciles <- data.table(gene_id = mcols(subject)$gene_id, decile = mcols(subject)$decile, fraction_in_untreated = as.vector(assays(subject)$mCG_fraction))
  
  dt <- merge(x = deseq_df, y = untreated_deciles, by = "gene_id", all.x = TRUE)
  dt <- dt[!is.na(dt$decile) & !is.na(dt$log2FoldChange) & !is.na(dt$padj),]
  
  results <- data.table()
  fdrlevel.de <- 0.05
  for (x in 1:10) {
    dt_temp <- dt[dt$decile == x,]
    dt_temp <- dt_temp[!is.na(dt_temp$log2FoldChange),]
    dt_temp <- dt_temp[!is.na(dt_temp$padj),]
    number_upregulated <- nrow(dt_temp[dt_temp$padj <= fdrlevel.de & dt_temp$log2FoldChange > 0,])
    number_downregulated <- nrow(dt_temp[dt_temp$padj <= fdrlevel.de & dt_temp$log2FoldChange < 0,])
    number_neutral <- nrow(dt_temp[dt_temp$padj > fdrlevel.de | dt_temp$log2FoldChange == 0,])
    temp_results <- data.table(decile = rep(x, 3),
                               count = c(number_upregulated, number_downregulated, number_neutral),
                               direction = c("up", "down", "neutral"))
    results <- rbind(results, temp_results)
  }
  results <- results[results$direction!="neutral"]
  output <- ggplot(results, aes(fill=direction, y=count, x=decile)) + 
    geom_bar(position="dodge", stat="identity") +
    xlab(paste(subject_name, "m (%) decile")) + 
    ylab("No. DE GENES\nDMSO vs 5-azacytidine") +
    ggtitle(figure_title) +
    theme_bw()
  ggsave(file.path(OUTPUT_DIR, paste(subject_name, figure_title, "methylation_decile_vs_change_with_aza.pdf", sep="_")),
         plot = output,
         width = 8,
         height = 5,
         units = "in")
}

BS <- readRDS("CG_BS.rds")
middle_base=3
CG_BS_CG <- BS[rowRanges(BS)$dinucleotide=="CG"]

compare_methylation_from_emseq_vs_expression_change_in_drug_trt(CG_BS_CG, figure_title = "CG")