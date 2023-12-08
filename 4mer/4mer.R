library(ggplot2)
library(ggseqlogo)

OUTPUT_DIR <- "./4mer"

FASTA_INDEX_PATH = "./inputs/genomes/Apar_genome_pUC19_lambda.fasta.fai"
BLACKLISTED_REGION <- "./inputs/annotations/Blacklisted.regions.0based.bed"
BLACKLISTED_REGION_GRANGES <- import.bed(BLACKLISTED_REGION)

##########################

dinucleotides <- c("a", "c", "t", "g")

# This is a recursive function to generate a kmer "tree": a list of lists describing nucleotide contexts.
# It is nested to allow for efficient subsetting of the BS object.
generate_kmer_tree <- function(kmer="", 
                               max_kmer_length = 4, 
                               position_of_methylated_nucleotide_in_kmer = 2,
                               methylated_nucleotide = "c"){
  position_of_methylated_nucleotide_in_kmer <- 2
  if (nchar(kmer) >= max_kmer_length){
    names(kmer) <- kmer
    return(kmer)
  }else{
    if (nchar(kmer) == position_of_methylated_nucleotide_in_kmer-1) {
      daughter_kmer_list <- list(paste0(kmer, methylated_nucleotide))
    } else {
      daughter_kmer_list <- lapply(1:4, function(ith_dinucleotide){
        daughter_kmer <- paste0(kmer,
                                dinucleotides[ith_dinucleotide])
        return(daughter_kmer)
      })
    }
    daughter_kmer_recursive_list <- lapply(daughter_kmer_list, function(daughter_kmer){
      generate_kmer_tree(daughter_kmer, 
                         max_kmer_length = max_kmer_length, 
                         position_of_methylated_nucleotide_in_kmer = position_of_methylated_nucleotide_in_kmer, 
                         methylated_nucleotide = methylated_nucleotide)
    })
    names(daughter_kmer_recursive_list) <- daughter_kmer_list
    return(daughter_kmer_recursive_list)
  }
}


tetramer_tree <- generate_kmer_tree()

#This is the function to carry out on each nucleotide context, to extract relevant statistics for that context.
#mC is the total number of methylated calls / the total number of calls as a fraction
subsetted_bs_statistics <- function(BS, context){
  C <- sum(assays(BS)$M)
  CT <- sum(assays(BS)$Cov)
  count <- nrow(BS)
  min = min(assays(BS)$Fraction)
  lowerquant = quantile(assays(BS)$Fraction, 0.25, names = FALSE)
  median = median(assays(BS)$Fraction)
  upperquant = quantile(assays(BS)$Fraction, 0.75, names = FALSE)
  max = max(assays(BS)$Fraction)
  return(data.frame(context = context,
                    mC = C/CT,
                    n = count))
}

#This function recursively searches through the kmer tree generated above. 
#At each level, if that level is a list then it goes a level deeper. 
#If it is not a list, then it extracts the statistics using the function above
#width_of_context here refers to the length of the context in the BS object (in this instance, 5)
extract_kmer_stats <- function(kmer_tree, BS, elementC, width_of_context){
  middle_base <- (width_of_context+1)/2
  parent_kmer <- substr(names(kmer_tree)[1], 0, (nchar(names(kmer_tree)[1])-1))
  kmer_start_position <- middle_base - elementC + 1
  kmer_end_position <- kmer_start_position + nchar(parent_kmer)
  if (parent_kmer != "") {
    BS <- BS[substr(rowRanges(BS)$context, kmer_start_position, kmer_end_position-1) == parent_kmer]
  }
  if(!is.list(kmer_tree)){
    BS <- BS[substr(rowRanges(BS)$context, kmer_start_position, kmer_end_position) == kmer_tree]
    print(names(kmer_tree))
    return(subsetted_bs_statistics(BS, kmer_tree))
  }else{
    return(rbindlist(lapply(kmer_tree, function(kmer_subtree){
      return(extract_kmer_stats(kmer_subtree, BS, elementC=2, width_of_context=width_of_context))
    })))
  }
}

CG_BS <- readRDS("CG_BS.rds")
CG_BS <- subsetByOverlaps(CG_BS, BLACKLISTED_REGION_GRANGES, invert = TRUE, ignore.strand=TRUE)
CG_BS <- CG_BS[rowRanges(CG_BS)$source=="Organism"]

whole_genome_results <- extract_kmer_stats(kmer_tree=tetramer_tree, BS=CG_BS, elementC = 2, width_of_context = 5)

whole_genome_results <- whole_genome_results[order(whole_genome_results$mC, decreasing = TRUE)]

whole_genome_results_CG <- whole_genome_results[substr(whole_genome_results$context, 3, 3)=="g"]

bar <- ggplot(whole_genome_results_CG, mapping = aes(x = toupper(context), y=mC*100)) +
  geom_col() + ylim(0,100) + xlab("Nucleotide context") + ylab("Global mCG (%)") + 
  theme_bw() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

ggsave(file.path(OUTPUT_DIR, "tetramer_mC_barchart.pdf"), bar)