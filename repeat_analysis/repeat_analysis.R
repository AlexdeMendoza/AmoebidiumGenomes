OUTPUT_DIR <- "./repeat_analysis"

CG_BS <- readRDS("CG_BS.rds")
#This is a simple function to retrieve methylation statistics for each element of a GRanges
#This uses a BSseq object which describes coverage and methylation levels at each cytosine base.
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

REPEAT_GTF_PATH <- "./inputs/annotations/Apar_genome.fasta.RM.bed.gtf.cleaned.final.kimura"
FASTA_INDEX_PATH = "./inputs/genomes/Apar_genome_pUC19_lambda.fasta.fai"

index2slengths <- function(FASTA_INDEX_PATH){
  FASTA_INDEX <- fread(FASTA_INDEX_PATH)
  SLENGTHS <- FASTA_INDEX$V2
  names(SLENGTHS) <- FASTA_INDEX$V1
  return(SLENGTHS)
}

##### The first step is to generate a Summarized Experiment (SE) object of repeats #####
#This will include kimura distance, methylated calls, and total calls for each element.

SLENGTHS <- index2slengths(FASTA_INDEX_PATH)

#Create the repeat SE:
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
CURATED_REPEATS <- repeat_annotation_summarizedexperiment(CG_BS, REPEAT_GTF_PATH, SLENGTHS)

#Blacklist (assembly artifacts) granges:
BLACKLISTED_REGION_GRANGES <- import.bed(BLACKLISTED_REGION)

middle_base=3
CG_BS_CG <- CG_BS[rowRanges(CG_BS)$dinucleotide=="CG"]
BS_cgc_CG <- CG_BS_CG[substr(rowRanges(CG_BS_CG)$context, middle_base, middle_base+2) == "cgc" |
                     substr(rowRanges(CG_BS_CG)$context, middle_base-1, middle_base+1) == "gcg"]
BS_non_cgc_CG <- CG_BS_CG[!(substr(rowRanges(CG_BS_CG)$context, middle_base, middle_base+2) == "cgc" |
                           substr(rowRanges(CG_BS_CG)$context, middle_base-1, middle_base+1) == "gcg")]

##### First make the kimura distribution plots (no methylation data needed): ##### 

#Filter repeats to those which are >500bp
CURATED_REPEATS <- CURATED_REPEATS[width(CURATED_REPEATS)>500]

#Only look at those which we could identify as belonging to defined TE family:
CURATED_REPEATS <- CURATED_REPEATS[rowData(CURATED_REPEATS)$family %in% c("DNA", "LINE", "LTR")]

#Remove those which lie within either endogenised viral regions, or assembly artifacts (blacklist)
CURATED_REPEATS <- subsetByOverlaps(CURATED_REPEATS, BLACKLISTED_REGION_GRANGES, invert = TRUE)

#Extract data for tables:
kimura_data <- data.table(Family = rowRanges(CURATED_REPEATS)$family,
                          kimura = rowRanges(CURATED_REPEATS)$kimura)
kimura_histogram_split_by_family <- ggplot(data=kimura_data, mapping = aes(x = kimura)) +
  geom_histogram(binwidth=1) +
  facet_wrap(Family ~ .) + 
  theme_bw() + 
  xlab("Kimura score") + 
  ylab("Count")
kimura_histogram <- ggplot(data=kimura_data, mapping = aes(x = kimura)) +
  geom_histogram(binwidth=1) + 
  theme_bw() + 
  xlab("Kimura score") + 
  ylab("Count")
ggsave(filename = "kimura_histogram.pdf", path = OUTPUT_DIR, plot = kimura_histogram, width = 75, height = 25, units = "cm")
ggsave(filename = "kimura_histogram_split_by_family.pdf", path = OUTPUT_DIR, plot = kimura_histogram_split_by_family, width = 75, height = 25, units = "cm")

cut_kimura_into_ranges <- function(kimura_vector){
  output <- sapply(kimura_vector, function(element){
    if(element <= 5){
      return("0-5")
    }else if(element <= 10){
      return("5-10")
    }else if(element <= 15){
      return("10-15")
    }else if(element <= 20){
      return("15-20")
    }else if(element <= 25){
      return("20-25")
    }else if(element <= 30){
      return("25-30")
    }else if(element <= 35){
      return("30-35")
    }else if(element <= 40){
      return("35-40")
    }else if(element > 40){
      return(">40")
    }else{
      return(NA)
    }
  })
  output <- factor(output, levels = c("0-5",
                                      "5-10",
                                      "10-15",
                                      "15-20",
                                      "20-25",
                                      "25-30",
                                      "30-35",
                                      "35-40",
                                      ">40"))
  return(output)
}

##### Now make figures relating kimura range to methylation #####

make_repeat_figures <- function(BS, figure_title){
  CURATED_REPEATS <- repeat_annotation_summarizedexperiment(BS,
                                                    REPEAT_GTF_PATH = REPEAT_GTF_PATH, 
                                                    SLENGTHS = SLENGTHS)
  
  #Exclude repeats that are <500bp, within viruses, and within the blacklisted regions (likely mitochondrial):
  CURATED_REPEATS <- CURATED_REPEATS[width(CURATED_REPEATS)>500]
  CURATED_REPEATS <- CURATED_REPEATS[rowData(CURATED_REPEATS)$family %in% c("DNA", "LINE", "LTR")]

  CURATED_REPEATS <- subsetByOverlaps(CURATED_REPEATS, BLACKLISTED_REGION_GRANGES, invert = TRUE, ignore.strand = TRUE)
  
  REPEAT_MC_DATA <- data.table(gene_id = rowData(CURATED_REPEATS)$gene_id,
                               class = rowData(CURATED_REPEATS)$class,
                               Family = rowData(CURATED_REPEATS)$family,
                               tpm = as.vector(assays(CURATED_REPEATS)$counts),
                               width=width(CURATED_REPEATS),
                               chr=as.vector(seqnames(CURATED_REPEATS)),
                               kimura=rowData(CURATED_REPEATS)$kimura)
  
  
  Global_methylation <- sum(assays(BS)$M)/sum(assays(BS)$Cov)
  REPEAT_MC_METHYLATION_DATA <- data.table(fraction = as.vector(assays(CURATED_REPEATS)$mCG_fraction*100),
                                           reads = as.vector(assays(CURATED_REPEATS)$mCG_met),
                                           cov = as.vector(assays(CURATED_REPEATS)$mCG_cov),
                                           position_per_gene = as.vector(rowData(CURATED_REPEATS)$CG_position_per_gene))
  REPEAT_MC_DATA <- cbind(REPEAT_MC_DATA, REPEAT_MC_METHYLATION_DATA)
  
  REPEAT_MC_DATA <- REPEAT_MC_DATA[!is.na(REPEAT_MC_DATA$fraction),]
  REPEAT_MC_DATA <- REPEAT_MC_DATA[!is.na(REPEAT_MC_DATA$kimura),]
  
  REPEAT_MC_DATA$kimura_range <- cut_kimura_into_ranges(REPEAT_MC_DATA$kimura)
  
  kimura_range_vs_methylation <- ggplot(data=REPEAT_MC_DATA, mapping=aes(x=kimura_range, y=fraction)) +
    geom_boxplot(outlier.shape = NA, coef=0) +
    geom_violin(alpha=0) +
    theme_bw() + 
    xlab("Kimura score") + 
    ylab("Methylation (%)") + 
    ggtitle(paste("Context:", figure_title))
  
  kimura_range_vs_methylation_split_by_family <- ggplot(data=REPEAT_MC_DATA, mapping=aes(x=kimura_range, y=fraction)) +
    geom_boxplot(outlier.shape = NA, coef=0, fill="light grey") +
    geom_violin(alpha=0) +
    facet_wrap(. ~ Family) + 
    theme_bw() + 
    xlab("Kimura score") + 
    ylab("Methylation (%)") + 
    ggtitle(paste("Context:", figure_title)) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
  
  ggsave(filename = paste("m", figure_title, "_", "kimura_range_vs_methylation.pdf", sep=""), path = OUTPUT_DIR, plot = kimura_range_vs_methylation, width = 25, height = 25, units = "cm")
  ggsave(filename = paste("m", figure_title, "_", "kimura_range_vs_methylation_split_by_family.pdf", sep=""), path = OUTPUT_DIR, plot = kimura_range_vs_methylation_split_by_family, width = 75, height = 25, units = "cm")
}

make_repeat_figures(CG_BS_CG, figure_title = "CG")
make_repeat_figures(BS_cgc_CG, figure_title = "CGC")
make_repeat_figures(BS_non_cgc_CG, figure_title = "(non)CGC")