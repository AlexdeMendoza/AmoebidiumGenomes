library(remotes)
library(devtools)
library(data.table)
library(stringr)
library(ggplot2)
library(dplyr)
library(tidyr)
library(bsseq)
library(R.utils)
library(RIdeogram) #install.packages("RIdeogram")

OUTPUT_DIR <- "./ideogram"
FASTA_INDEX_PATH = "./inputs/genomes/Apar_genome_pUC19_lambda.fasta.fai"
SLENGTHS <- index2slengths(FASTA_INDEX_PATH)

#Prepare subset BS objects to CGC / [ACT]CG[AGT] contexts
CG_BS <- readRDS("CG_BS.rds")
middle_base=3
CG_BS_CG <- CG_BS[rowRanges(CG_BS)$dinucleotide=="CG"]
BS_cgc_CG <- CG_BS_CG[substr(rowRanges(CG_BS_CG)$context, middle_base, middle_base+2) == "cgc" |
                     substr(rowRanges(CG_BS_CG)$context, middle_base-1, middle_base+1) == "gcg"]
BS_non_cgc_CG <- CG_BS_CG[!(substr(rowRanges(CG_BS_CG)$context, middle_base, middle_base+2) == "cgc" |
                           substr(rowRanges(CG_BS_CG)$context, middle_base-1, middle_base+1) == "gcg")]

#Set up functions
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
sliding_window <- function(SLENGTHS,
                           BS,
                           WINDOW_SIZE){
  ORG_SLENGTHS <- SLENGTHS[!names(SLENGTHS) %in% c("NC_001416.1", "M77789.2")]
  #Create GRanges object describing each organism scaffold
  contig_gr <- GRanges(seqnames = names(ORG_SLENGTHS), ranges = IRanges(start = 1, end = ORG_SLENGTHS), seqlengths = ORG_SLENGTHS)
  
  #Create GRanges object describing each scaffold, split into windows of a set width
  sliding_windows_split <- slidingWindows(contig_gr, width = WINDOW_SIZE, step = WINDOW_SIZE/5) %>% unlist()
  
  #Get methylation statistics based on this sliding windows GRanges object
  Sliding_windows_mC <- get_mc_stats_from_granges(BS, sliding_windows_split)
  
  colData <- DataFrame(Treatment=TREATMENTS,
                       row.names=SAMPLES)
  SE <- SummarizedExperiment(assays=list(mCG_met = Sliding_windows_mC[[2]],
                                         mCG_cov = Sliding_windows_mC[[3]],
                                         mCG_fraction=Sliding_windows_mC[[4]]),
                             rowRanges=sliding_windows_split, colData=colData)
  rowRanges(SE)$CG_position_per_gene <- Sliding_windows_mC[[1]]
  
  
  return(SE)
}

virus_figure_choices <- data.table(Type = c("Adintovirus", "Pandoravirus", "Plavaka"),
                                   Shape = c("box", "circle", "triangle"),
                                   color = c("D81B60", "1E88E5", "FFC107"))

VIRAL_GTF_PATH <- "./inputs/annotations/GiantRepeats.gtf"

#Load in GRanges that describes viral loci
viral_GR <- viral_annotation_granges(VIRAL_GTF_PATH, SLENGTHS)

#Get data.table that describes viral loci, and merge with info on visualisation choices
virus <- data.table(Chr=as.vector(seqnames(viral_GR)),
                    Start=start(viral_GR),
                    End=end(viral_GR),
                    Type = viral_GR$family)

virus <- merge(virus, virus_figure_choices, by='Type')

#Subset to those that are in the 18 chromosomes:
virus <- virus[virus$Chr %in% names(SLENGTHS)[1:18],]
virus <- virus[,c(1,5,2,3,4,6)]
karyotype<-data.table(Chr=names(SLENGTHS),Start=0,End=SLENGTHS)
karyotype <- karyotype[1:18]

ideogram_figures <- function(BS, title, use_ch=FALSE){
  chromosomal_bs <- BS[as.vector(seqnames(BS)) %in% names(SLENGTHS)[1:18]]
  methylation_windows <- sliding_window(SLENGTHS, chromosomal_bs, 10000)
  methylation <- data.table(Chr=as.vector(seqnames(methylation_windows)),
                              Start=start(methylation_windows),
                              End=end(methylation_windows),
                              Value=as.vector(assays(methylation_windows)$mCG_fraction),
                              Color="000000")
  methylation$Value <- methylation$Value*100
  methylation <- methylation[!is.na(methylation$Value),]
  ideogram(karyotype, overlaid=methylation, label = virus, label_type = "marker", output = file.path(OUTPUT_DIR, paste("ideogram_virus_m", title, ".svg", sep="")), colorset1 = c("#ffffff", "#ff7f7f", "#ff0000"))
}
ideogram_figures(BS_cgc_CG, title = "CGC")
ideogram_figures(BS_non_cgc_CG, title = "(non)CGC")