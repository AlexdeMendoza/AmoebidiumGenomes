library(seqinr)
library(bsseq)

setwd("/data/SBCS-ademendoza/02-lukesarre/amoebidium_paper_trimmed/")

#Define functions:
build_cg_loci_from_fasta <- function(ORGANISM_FASTA_PATH, width_of_context=5){
  MET <- read.fasta(MET_FASTA_PATH)
  UNMET <- read.fasta(UNMET_FASTA_PATH)
  FASTA <- read.fasta(ORGANISM_FASTA_PATH)
  FASTA <-c(MET, UNMET, FASTA)
  distance_from_center <- (width_of_context-1)/2
  positions_to_find <- seq(-distance_from_center, distance_from_center)
  CG_maps <- lapply(FASTA, function(chr){
    chr_vector <- c(rep("n",(distance_from_center)), chr, rep("n",(distance_from_center)))
    print(attr(chr, "name"))
    C_coords <- seq_along(chr_vector)[!is.na(chr_vector) & chr_vector=="c"]
    C_list_offset_bases <- lapply(positions_to_find, function(position){
      positions <- C_coords + position
      return(chr_vector[positions])
    })
    C_upstream_downstream <- do.call(paste0, C_list_offset_bases)
    G_coords <- seq_along(chr_vector)[!is.na(chr_vector) & chr_vector=="g"]
    G_list_offset_bases <- lapply(positions_to_find, function(position){
      positions <- G_coords - position
      return(comp(chr_vector[positions], ambiguous = TRUE))
    })
    G_upstream_downstream <- do.call(paste0, G_list_offset_bases)
    if (length(C_coords) != 0) {
      C <- data.frame(chr=attr(chr, "name"),
                      position=C_coords-(distance_from_center),
                      base="C",
                      dinucleotide=toupper(substr(C_upstream_downstream, (distance_from_center)+1, (distance_from_center)+2)),
                      trinucleotide=toupper(substr(C_upstream_downstream, (distance_from_center)+1, (distance_from_center)+3)),
                      context=C_upstream_downstream)
    }else{C <- NULL}
    if (length(G_coords) != 0) {
      G <- data.frame(chr=attr(chr, "name"),
                      position=G_coords-(distance_from_center),
                      base="G",
                      dinucleotide=toupper(substr(G_upstream_downstream, (distance_from_center)+1, (distance_from_center)+2)),
                      trinucleotide=toupper(substr(G_upstream_downstream, (distance_from_center)+1, (distance_from_center)+3)),
                      context=G_upstream_downstream)
    }else{G <- NULL}
    
    CG_map <- rbind(C, G)
    return(CG_map)
  })
  CG_maps <- rbindlist(CG_maps)
  return(CG_maps)
}

loadcgmap <- function(CG_LOCI=CG_LOCI, CGMAP_PATH = CGMAP_PATH, TREATMENTS = "NA", SAMPLES = "NA"){
  #Rename the columns of CG_map
  #Fix this being a for loop in future
  for (x in seq_along(SAMPLES)) {
    CGMAP_TEMP <- fread(CGMAP_PATH[x], select = c("V1", "V3", "V6", "V7", "V8"), key = c("V1", "V3"))
    #Rename the columns of CG_map
    CGMAP_TEMP <- dplyr::rename(CGMAP_TEMP, 
                                !!paste("Fraction_Met", x, sep = "") := V6, 
                                !!paste("C_Reads", x, sep = "") := V7, 
                                !!paste("CT_Reads", x, sep = "") := V8)
    CG_LOCI <- dplyr::left_join(CG_LOCI, CGMAP_TEMP, by = c("chr" = "V1", "position" = "V3"))
  }
  #Add new column (strand), saying which strand the C is on
  
  CG_LOCI$strand <- ifelse(CG_LOCI$base == "G", "-", "+")
  
  #### Sort the CGmap by organism
  CG_LOCI$source <- ifelse(CG_LOCI$chr == "M77789.2", "MetCntrl", ifelse(CG_LOCI$chr == "NC_001416.1", "UnmetCntrl", "Organism"))
  nsamples=length(SAMPLES)
  # Take these vectors as input to the array.
  Fraction_Met <- unlist(lapply(c(1:nsamples), FUN=function(x){
    return(CG_LOCI[[paste("Fraction_Met", x, sep = "")]])
  }))
  C_Reads <- unlist(lapply(c(1:nsamples), FUN=function(x){
    return(CG_LOCI[[paste("C_Reads", x, sep = "")]])
  }))
  CT_Reads <- unlist(lapply(c(1:nsamples), FUN=function(x){
    return(CG_LOCI[[paste("CT_Reads", x, sep = "")]])
  }))
  array <- array(c(Fraction_Met, C_Reads, CT_Reads), 
                 dim = c(nrow(CG_LOCI),
                         nsamples,
                         3), 
                 dimnames = list(NULL,
                                 SAMPLES,
                                 c("Fraction_Met", "C_Reads", "CT_Reads")))
  CG_LOCI <- dplyr::select(CG_LOCI, 
                           chr, 
                           position, 
                           dinucleotide, 
                           trinucleotide, 
                           context,
                           strand, 
                           source)
  rm(CT_Reads, C_Reads, Fraction_Met)
  array[is.na(array)] <- 0
  return(list(CG_LOCI, array))
}

build_bsseq <- function(SLENGTHS, CGMAP, METH_ARRAY){
  Gr_Obj <- GRanges(seqnames = CGMAP$chr, 
                    ranges = IRanges(start = as.numeric(CGMAP$position), end = as.numeric(CGMAP$position)), 
                    strand = CGMAP$strand,
                    dinucleotide = CGMAP$dinucleotide,
                    trinucleotide = CGMAP$trinucleotide,
                    context = CGMAP$context,
                    seqlengths = SLENGTHS,
                    source = CGMAP$source)
  if(length(SAMPLES) == 1){
    BSObj <- BSseq(gr = Gr_Obj,
                   sampleNames = SAMPLES,
                   M = as.matrix(METH_ARRAY[,,2]),
                   Cov = as.matrix(METH_ARRAY[,,3]),
                   rmZeroCov = FALSE)
    assays(BSObj, withDimnames = FALSE)$Fraction <- as.matrix(METH_ARRAY[,,1])
    assays(BSObj)
  }else{
    BSObj <- BSseq(gr = Gr_Obj,
                   sampleNames = SAMPLES,
                   M = METH_ARRAY[,,2],
                   Cov = METH_ARRAY[,,3],
                   rmZeroCov = FALSE)
    assays(BSObj)$Fraction <- METH_ARRAY[,,1]
  }
  rowData(BSObj)$source <- CGMAP$source
  return(BSObj)
}

index2slengths <- function(FASTA_INDEX_PATH){
  FASTA_INDEX <- fread(FASTA_INDEX_PATH)
  SLENGTHS <- FASTA_INDEX$V2
  names(SLENGTHS) <- FASTA_INDEX$V1
  return(SLENGTHS)
}

#List of inputs:
ORGANISM_FASTA_PATH <- "./inputs/genomes/amoebidium/Apar_genome.fasta"
MET_FASTA_PATH <- "./inputs/genomes/pUC19.fasta"
UNMET_FASTA_PATH <- "./inputs/genomes/lambda.fasta"
TREATMENTS <- "NA"
SAMPLES <- "NA"
CGMAP_PATH <- c("./inputs/CGmaps/Apar_feb_assembly_merged_initial_emseq.CGmap.gz")
FASTA_INDEX_PATH <- "./inputs/genomes/Apar_genome_pUC19_lambda.fasta.fai"

#This describes the length of each chromosome / contig:
SLENGTHS <- index2slengths(FASTA_INDEX_PATH)

#This generates a large dataframe that describes the location of each C nucleotide, with strand and nucleotide context
CG_LOCI <- build_cg_loci_from_fasta(ORGANISM_FASTA_PATH, width_of_context=5)

#Here we append methylation data to these loci (structure of function allows for adding multiple CGmaps, if needed)
CGMAP <- loadcgmap(CG_LOCI, CGMAP_PATH, TREATMENTS, SAMPLES)
CG_METH_ARRAY <- CGMAP[[2]]
CGMAP_SUMMARY <- CGMAP[[1]]

#Build a BS object
CG_BS <- build_bsseq(SLENGTHS, CGMAP_SUMMARY, CG_METH_ARRAY)
saveRDS(CG_BS, "CG_BS.rds")
