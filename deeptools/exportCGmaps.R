export_cg_map <- function(BS, output_path){
  CGmap_df <- data.frame(chr = as.vector(seqnames(BS)),
                         nucleotide = ifelse(strand(BS)=="+", "C", "G"),
                         base = start(BS),
                         context = ifelse(rowRanges(BS)$dinucleotide=="CG",
                                          "CG",
                                          "CH"),
                         dinucleotide = rowRanges(BS)$dinucleotide,
                         fraction = assays(BS)$Fraction,
                         meth = assays(BS)$M,
                         cov = assays(BS)$Cov)
  write.table(CGmap_df, file = output_path, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
}

#Load BS object
BS <- readRDS("CG_BS.rds")

#Export CGmaps for different subsets of nucleotide context, for each sample
width_of_context <- nchar(rowRanges(BS[1])$context)
middle_base <- (width_of_context+1)/2

BS_cgc <- BS[substr(rowRanges(BS)$context, middle_base, middle_base+2) == "cgc" |
               substr(rowRanges(BS)$context, middle_base-1, middle_base+1) == "gcg" |
               substr(rowRanges(BS)$context, middle_base-2, middle_base) == "cgc"]
BS_non_cgc <- BS[!(substr(rowRanges(BS)$context, middle_base, middle_base+2) == "cgc" |
                     substr(rowRanges(BS)$context, middle_base-1, middle_base+1) == "gcg" |
                     substr(rowRanges(BS)$context, middle_base-2, middle_base) == "cgc")]

export_cg_map(BS_cgc, file.path(getwd(), "deeptools/CGmaps", "BS_cgc.CGmap"))
export_cg_map(BS_non_cgc, file.path(getwd(), "deeptools/CGmaps", "BS_non_cgc.CGmap"))