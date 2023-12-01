library(data.table)
library(stringr)
library(ggplot2)
library(dplyr)
library(tidyr)
library(bsseq)
library(R.utils)
library(reshape2)
library(ggseqlogo)
library(cowplot)
library(dmrseq)
library("topGO")


#####################################################
### Read CGmap files into bsseq
#####################################################

#this function reads a CGmap file, selects the CpGs, collapses the CpG information (each C in the CpG is merged into a single value)
# makes bsseq objects (type of file useful for methylation analysis) for the nuclear genome
# saves bsseq objects as "RDS" files, these are preprocessed binarised files that you can load later on very quickly
# it also saves a bsseq file for the lambda genome in rds format

read_CGmap_into_CG_lambda_files <- function(CGmap, name){
  
  # Read the file
  dat <- fread(input = CGmap, sep = "\t", select = c(1,2,3,4,5,7,8),
               col.names = c("chr", "base", "position", "context","dinucleotide",
                             "C_reads", "CT_reads"))
  # Subset to CG context only from the genome, not spike ins (pUC19 or phage lambda)
  message(paste0("processing CG..."))
  datCG <- dat[dat$context == "CG" & !(dat$chr %in% c("chrL","pUC19","M77789.2","NC_001416.1")), ]
  
  datCG$strand <- ifelse(test = datCG$base == "G",
                         yes = "-",
                         no = "+")
  
  bs_obj <- BSseq(gr = GRanges(seqnames = datCG$chr,
                               ranges = IRanges(start = as.numeric(datCG$position),
                                                end = as.numeric(datCG$position)), strand = datCG$strand ), 
                  sampleNames = name, M = as.matrix(datCG$C_reads), Cov = as.matrix(datCG$CT_reads),
                  rmZeroCov = FALSE)
  
  bs_obj_collapsed <- strandCollapse(bs_obj)
  rm(datCG)
  
  saveRDS(object = bs_obj_collapsed, file = paste0(name,".CG_bsseq.rds"))
  
  ###########
  message(paste0("processing Lambda..."))
  # Subset chrL
  dat_Lambda <- dat[dat$chr == "NC_001416.1", ]
  
  #add_sp_name
  dat_Lambda$species <- name
  
  #save lambda genome
  saveRDS(object = dat_Lambda, file = paste0(name,".lambda.rds"))
  
  ###########
  message(paste0("processing Mitochondria..."))
  # Subset chrL
  dat_mito <- dat[dat$chr == "chrM", ]
  
  #add_sp_name
  dat_mito$species <- name
  
  #save lambda genome
  saveRDS(object = dat_mito, file = paste0(name,".mito.rds"))
  
  ###########
  message(paste0("processing pUC19..."))
  # Subset chrL
  dat_puc <- dat[dat$chr == "M77789.2", ]
  
  #add_sp_name
  dat_puc$species <- name
  
  #save lambda genome
  saveRDS(object = dat_puc, file = paste0(name,".pUC19.rds"))
  
  gc()
  
  # return the aggregated data object
  message(paste0("Objects created for ",name))
}


#######################################################################
### DMR calling for Amoebidium appalachense developmental EM-seq samples
#######################################################################

# The CGmap files are available from GEO submission: GSEXXXX.

read_CGmap_into_CG_lambda_files(CGmap = "EMseq_Aapp_Dev_5h_rep1.CGmap.gz", 
                                name = "Apar_5h_rep1")
read_CGmap_into_CG_lambda_files(CGmap = "EMseq_Aapp_Dev_5h_rep2.CGmap.gz", 
                                name = "Apar_5h_rep2")
read_CGmap_into_CG_lambda_files(CGmap = "EMseq_Aapp_Dev_14h_rep1.CGmap.gz", 
                                name = "Apar_14h_rep1")
read_CGmap_into_CG_lambda_files(CGmap = "EMseq_Aapp_Dev_14h_rep2.CGmap.gz", 
                                name = "Apar_14h_rep2")
read_CGmap_into_CG_lambda_files(CGmap = "EMseq_Aapp_Dev_20h_rep1.CGmap.gz", 
                                name = "Apar_20h_rep1")
read_CGmap_into_CG_lambda_files(CGmap = "EMseq_Aapp_Dev_20h_rep2.CGmap.gz", 
                                name = "Apar_20h_rep2")
read_CGmap_into_CG_lambda_files(CGmap = "EMseq_Aapp_Dev_33h_rep1.CGmap.gz", 
                                name = "Apar_33h_rep1")
read_CGmap_into_CG_lambda_files(CGmap = "EMseq_Aapp_Dev_33h_rep2.CGmap.gz", 
                                name = "Apar_33h_rep2")


obj_files <- c("Apar_5h_rep1.CG_bsseq.rds",
               "Apar_5h_rep2.CG_bsseq.rds",
               "Apar_14h_rep1.CG_bsseq.rds",
               "Apar_14h_rep2.CG_bsseq.rds",
               "Apar_20h_rep1.CG_bsseq.rds",
               "Apar_20h_rep2.CG_bsseq.rds",
               "Apar_33h_rep1.CG_bsseq.rds",
               "Apar_33h_rep2.CG_bsseq.rds")


# Read a Bs_seq boject from .Rds file

obj_list <- lapply(X = obj_files, readRDS)
bs_all <- bsseq::combineList(x = obj_list)


pData(bs_all)$condition <- factor(c(rep("5h", times=2), rep("14h", times=2),rep("20h", times=2), rep("33h", times=2)))
pData(bs_all)$Replicate <- c(1:2,1:2,1:2,1:2)

saveRDS(object = bs_all, file = "all_methylomes_bsseq.rds")

bs_all <- readRDS("all_methylomes_bsseq.rds")

T <- getCoverage(bs_all, type="Cov", what = "perBase")
C <- getCoverage(bs_all, type="M", what = "perBase")

# obtain global levels for each sample an plot

globalLevels <- data.frame(sample = colnames(C),
           mCG = colSums(C)/colSums(T),
          replicate = rep(x = c("rep1","rep2"), times = 4), stage = c("5h","5h","14h","14h","20h","20h","33h","33h"))

globalLevels$stage <- factor(globalLevels$stage , levels = c("5h","14h","20h","33h"))

gg_dev_mC <- ggplot(globalLevels, aes(x = stage, y = 100*mCG, color = replicate)) + geom_point() +
  theme_bw() + ylim(c(0, 50)) + ylab("mCG %")

ggsave(gg_dev_mC, filename = "Developmental_mCG.pdf", width = 3, height = 2.5)


# This function adds builds a "loci" from a GRanges object

add_loci <- function(gr){
  
  loci <- str_c(seqnames(gr), start(gr), sep=":") %>%
    str_c(end(gr), sep = "-")
  
  gr$loci <- loci
  
  return(gr)
  
}

# Call DMRs with dmrseq


call_DMRs_from_object_list <- function(obj_files, condition1, condition2){
  # Load the data
  obj_list <- lapply(X = obj_files, readRDS)
  obj_list <- bsseq::combineList(x = obj_list)
  
  pData(obj_list)$condition <- factor(c(rep(condition1, times=2), rep(condition2, times=2)))
  pData(obj_list)$Replicate <- c(1:2,1:2)
  
  # Remove CpG with no coverage
  loci.idx <- which(DelayedMatrixStats::rowSums2(getCoverage(obj_list, type="Cov")==0) == 0)
  obj_list <- obj_list[loci.idx, ]

  regions <- dmrseq(obj_list, testCovariate = "condition", bpSpan = 500, maxGap = 500, maxPerms = 20, chrsPerChunk = 1)
  
  T <- getCoverage(obj_list, regions = regions, type="Cov", what = "perRegionTotal")
  C <- getCoverage(obj_list, regions = regions, type="M", what = "perRegionTotal")
  regions$delta <- rowMeans(C[,c(1,2)]) - rowMeans(C[,c(3,4)])
  gc()
  return(regions)
}


obj_files <- c("Apar_5h_rep1.CG_bsseq.rds",
               "Apar_5h_rep2.CG_bsseq.rds",
               "Apar_20h_rep1.CG_bsseq.rds",
               "Apar_20h_rep2.CG_bsseq.rds")

dmrs_5_to_20 <- call_DMRs_from_object_list(obj_files = obj_files, condition1 = "5h", condition2 = "20h")

write.table( x = data.frame(dmrs_5_to_20), file = "~/Dropbox/Amoebidium_genome/Chromosome_scale/DMRseq_5h_vs_20h.tsv",
             row.names = F, quote = F, sep = "\t")


fread("~/Dropbox/Amoebidium_genome/Chromosome_scale/DMRseq_5h_vs_20h.tsv") %>% 
  filter(pval < 0.05, beta > 5)

obj_files <- c("Apar_5h_rep1.CG_bsseq.rds",
               "Apar_5h_rep2.CG_bsseq.rds",
               "Apar_33h_rep1.CG_bsseq.rds",
               "Apar_33h_rep2.CG_bsseq.rds")

dmrs_5_to_33 <- call_DMRs_from_object_list(obj_files = obj_files, condition1 = "5h", condition2 = "33h")
write.table(x = data.frame(dmrs_5_to_33), file = "~/Dropbox/Amoebidium_genome/Chromosome_scale/DMRseq_5h_vs_33h.tsv",
row.names = F, quote = F, sep = "\t")

obj_files <- c("Apar_5h_rep1.CG_bsseq.rds",
               "Apar_5h_rep2.CG_bsseq.rds",
               "Apar_14h_rep1.CG_bsseq.rds",
               "Apar_14h_rep2.CG_bsseq.rds")

dmrs_5_to_14 <- call_DMRs_from_object_list(obj_files = obj_files, condition1 = "5h", condition2 = "14h")
write.table(x = data.frame(dmrs_5_to_14), file = "~/Dropbox/Amoebidium_genome/Chromosome_scale/DMRseq_5h_vs_14h.tsv",
            row.names = F, quote = F, sep = "\t")

obj_files <- c("Apar_14h_rep1.CG_bsseq.rds",
               "Apar_14h_rep2.CG_bsseq.rds",
               "Apar_20h_rep1.CG_bsseq.rds",
               "Apar_20h_rep2.CG_bsseq.rds")

dmrs_14_to_20 <- call_DMRs_from_object_list(obj_files = obj_files, condition1 = "14h", condition2 = "20h")
write.table(x = data.frame(dmrs_14_to_20), file = "~/Dropbox/Amoebidium_genome/Chromosome_scale/DMRseq_14h_vs_20h.tsv",
            row.names = F, quote = F, sep = "\t")


#######################################################################
### Processing the RNA-seq to visualize Giant Repeat reactivation upon 5-Azacytidine treatment
#######################################################################

# Read TPM matrix
tpm <- fread("Amoebidium.chr.v3.mRNA.tpm.tsv")
tpm <- tpm[grep(tpm$target_id, pattern = "\\.1"),]
tpm$gene_id <- tpm$target_id %>% str_remove(., pattern = "\\.1")

tpm_mtx <- tpm %>% dplyr::select(-gene_id, -target_id)
rownames(tpm_mtx) <- tpm$gene_id


# read the gene IDs from various repetitive types
adinto_ids <- scan(file = "DifferentialExpression/Adintovirus.genes.ids",what = "character")
pandora_ids <- scan(file = "DifferentialExpression/Pandoravirus.genes.ids",what = "character")
plavaka_ids <- scan(file = "DifferentialExpression/Plavaka.genes.ids",what = "character")
te_genes_ids <- scan(file = "DifferentialExpression/TE.ORFs.ids",what = "character")

# classify genes according to they expression in Aza samples (tpm > 1), or in any other sample.
tpm$expressed_in <- rowSums(tpm_mtx > 1)
tpm$expressed_in_Aza <- rowSums(tpm_mtx[,c("AZA1", "AZA2", "AZA3")] > 1)
active_aza <- tpm %>% data.frame() %>% filter(gene_id %in% c(adinto_ids, pandora_ids,plavaka_ids, te_genes_ids)) %>%
  mutate(expression_pattern = ifelse(expressed_in == 0, "Not_expressed", 
                                     ifelse(expressed_in == expressed_in_Aza, "Aza_only",
                                            ifelse(expressed_in > (expressed_in_Aza + 2), "Broad","Aza_most"))),
         type = ifelse(gene_id %in% adinto_ids, "Adintoviral", ifelse(gene_id %in% pandora_ids,"Pandoravirus",
                                                                      ifelse(gene_id %in% "Plavaka","Plavaka","TE ORFs") )))
  
  

GiantRepeatsSummaryExpr <- table(active_aza$expression_pattern,active_aza$type) %>% data.frame()


# obtain Pandoravirus genes that are broadly expressed (in samples other than the Azacytidine treated samples)
pandora_broad_ids <- active_aza %>% filter(type == "Pandoravirus", expression_pattern == "Broad") %>% .$gene_id


#######################################################################
### Obtaining gene ontologies for genes in endogenised regions
#######################################################################

# eggnog mapper (http://eggnog-mapper.embl.de/) was used to annotate the protein fasta file 

geneID2GO <- readMappings(file = "ReferenceGenome/Amoebidium.chr.v3.pep.Eggnog.GOs")
geneUniverse <- names(geneID2GO)

# read the gene IDs from various repetitive types
adinto_ids <- scan(file = "DifferentialExpression/Adintovirus.genes.ids",what = "character")
pandora_ids <- scan(file = "DifferentialExpression/Pandoravirus.genes.ids",what = "character")
plavaka_ids <- scan(file = "DifferentialExpression/Plavaka.genes.ids",what = "character")
te_genes_ids <- scan(file = "DifferentialExpression/TE.ORFs.ids",what = "character")

# function to obtain GOs

get_GOs <- function(geneNames, GOlist){
  # total number of genes with GOs
  geneUniverse <- names(GOlist)
  # logical for Genes in list
  list_for_GO <- factor(as.integer(geneUniverse %in% geneNames))
  names(list_for_GO) <- geneUniverse
  GOdata_BP <- new("topGOdata", description="cluster_1_program",
                   ontology="BP", allGenes=list_for_GO,
                   annot = annFUN.gene2GO, gene2GO = GOlist) 
  resultFisher_BP <- runTest(GOdata_BP,
                             algorithm="classic", statistic="fisher") 
  results_BP <- GenTable(GOdata_BP, classicFisher = resultFisher_BP, orderBy = "resultFisher", ranksOf = "classicFisher", topNodes = 30)
  results_BP$Ontology <- "BP"
  GOdata_MF <- new("topGOdata", description="cluster_1_program",
                   ontology="MF", allGenes=list_for_GO,
                   annot = annFUN.gene2GO, gene2GO = GOlist) 
  resultFisher_MF <- runTest(GOdata_MF,
                             algorithm="classic", statistic="fisher") 
  results_MF <- GenTable(GOdata_MF, classicFisher = resultFisher_MF, orderBy = "resultFisher", ranksOf = "classicFisher", topNodes = 30)
  results_MF$Ontology <- "MF"
  
  df <- rbind(results_BP, results_MF)
  return(df)
}

adinto_GOs <- get_GOs(geneNames = adinto_ids, GOlist = geneID2GO)
pandora_GOs <- get_GOs(geneNames = pandora_ids, GOlist = geneID2GO)
pandora_expressed_GOs <- get_GOs(geneNames = active_aza %>% dplyr::filter(type == "Pandoravirus", expression_pattern == "Broad") %>% .$gene_id, GOlist = geneID2GO)
pandora_reactivated_GOs <- get_GOs(geneNames = active_aza %>% dplyr::filter(type == "Pandoravirus", expression_pattern %in% c("Aza_only","Aza_most")) %>% .$gene_id, GOlist = geneID2GO)
plavaka_GOs <- get_GOs(geneNames = plavaka_ids, GOlist = geneID2GO)

write.table(x = pandora_GOs, file = "Pandora.GOs.tsv", 
            quote = FALSE, sep = "\t",
            row.names = FALSE, col.names = TRUE)

# visualize GO enrichments

plot_GOs <- function(df, name, GOnum = 10){
  # fix order
  
  df <- df %>% group_by(Ontology) %>% top_n(GOnum)
  
  df <- df %>% filter(as.numeric(classicFisher) < 0.05)
  
  df$id <- paste0(df$GO.ID,":",df$Term)
  
  df$id <- factor(df$id, levels = rev(df$id))
  
  # make plot
  ggbar <- ggplot(df, aes(x = id, y = -log10(as.numeric(classicFisher)), fill = Ontology )) +
    stat_summary(geom = "bar", fun = mean, position = "dodge") +
    xlab("Molecular Function // Biological process") +
    ylab("-log10(p-value)") +
    ggtitle(name) +
    theme_classic() +
    theme(
      legend.position='none',
      legend.background=element_rect(),
      plot.title=element_text(angle=0, size=12, face="bold", vjust=1),
      axis.text.x=element_text(angle=0, size=8, hjust=1.10),
      axis.text.y=element_text(angle=0, size=9, vjust=0.5),
      axis.title=element_text(size=12),
      legend.key=element_blank(),     #removes the border
      legend.key.size=unit(1, "cm"),      #Sets overall area/size of the legend
      legend.text=element_text(size=12),  #Text size
      title=element_text(size=12)) +
    guides(colour=guide_legend(override.aes=list(size=2.5))) +
    coord_flip()
  
  return(ggbar)
  
}

gg_pandora_GOs <- plot_GOs(df = pandora_GOs_slimmed, name = "Pandoravirus derived genes")
gg_pandora_expressed_GOs <- plot_GOs(df = pandora_expressed_GOs_slimmed, name = "Pandoravirus broadly expressed genes")


ggsave(gg_pandora_GOs, 
       filename = "GOs_pandoravirus.pdf",
       height = 9, width = 7)

ggsave(gg_pandora_expressed_GOs, 
       filename = "GOs_pandoravirus_broadly_expressed.pdf",
       height = 9, width = 7)

#######################################################################
### Obtaining Pfam domain enrichments for genes in endogenised regions
#######################################################################

# The PfamA annotation was obtained using hmmscan with this command:
# hmmscan --cut_ga --noali --acc --domtblout Amoebidium.chr.v3.pep.PfamA.domtblout Pfam-A.hmm Amoebidium.chr.v3.pep.fasta

pfam_numbers <- fread("ReferenceGenome/Amoebidium.chr.v3.pep.PfamA.domtblout.gz", sep = "\t", fill = TRUE, col.names = c("Pfam","PFID","gene_id","evalue")) %>% .[.$evalue < 0.001,]
viralrecall_numbers <- fread("ReferenceGenome/Amoebidium.chr.v3.pep.ViralRecallMarkers.domtblout.tsv", sep = "\t", fill = TRUE, col.names = c("Pfam","PFID","gene_id","evalue")) %>% .[.$evalue < 0.001,]

pfam_and_viral_numbers <- rbind( pfam_numbers, viralrecall_numbers )


pfamEnrichments <- function(ids, pfam_df = pfam_numbers){
  a <- pfam_df %>% filter(gene_id %in% ids) %>% group_by(Pfam) %>% tally() %>% dplyr::rename(test = n)
  b <- pfam_df %>% filter(!gene_id %in% ids) %>% group_by(Pfam) %>% tally() %>% dplyr::rename(background = n)
  
  df <- left_join(a,b) %>% mutate(background = ifelse(is.na(background), 0, background)) %>%
    mutate(tot1 = sum(test), tot2 = sum(background))
    
  df <- df %>% rowwise() %>%
    mutate(fishp=fisher.test(matrix(c(test,tot1-test,background,tot2-background ),ncol=2),alternative='greater')$p.value,
           odds=fisher.test(matrix(c(test,tot1-test,background,tot2-background ),ncol=2),alternative='greater')$estimate ) %>% arrange(odds)
  
  return(df)
  
}



pfamEnrichments(ids = adinto_ids)
pfamEnrichments(ids = plavaka_ids)


PfamFunctions <- read.xlsx("~/Dropbox/Amoebidium_genome/Chromosome_scale/ViralGenomeAnalysis/Repeats.xlsx", sheet = 6)

pandora_pfams <- pfamEnrichments(ids = pandora_ids, pfam_df = pfam_and_viral_numbers) %>% arrange(test) %>% filter(test > 5) %>%
  filter(!Pfam %in% c("Pox_VLTF3","Capsid_N","mcp","RVT_2","rve","gag_pre-integrs")) 
pandora_pfams <- left_join(pandora_pfams, PfamFunctions) %>% mutate(Function = ifelse(is.na(Function),"Other",Function))
pandora_pfams$Pfam <- factor(pandora_pfams$Pfam, levels= c(pandora_pfams$Pfam))


gg_pandora_pfam <- ggplot(pandora_pfams, aes(x = Pfam, y = -log10(fishp))) + 
  stat_summary(geom = "bar", fun = mean, position = "dodge") +
  xlab("Pfam domain") +
  ylab("-log10(p-value)") +
  ggtitle("GEVE domain enrichments") +
  theme_classic() + coord_flip()

gg_pandora_pfam2 <- ggplot(pandora_pfams, aes(x = Pfam, y = test, color = Function, size = -log10(fishp))) +
  geom_point(alpha = 0.7) +
  xlab("Pfam domain") +
  ylab("Pfam containing genes") +
  ggtitle("GEVE domain enrichments") + 
  theme_bw() + coord_flip()

ggsave(gg_pandora_pfam, 
       filename = "PFAMs_pandoravirus.pdf",
       height = 9, width = 7)

ggsave(gg_pandora_pfam2, 
       filename = "PFAMs_pandoravirus_dots.pdf",
       height = 4, width = 5)

pandora_expressed_ids <- active_aza %>% dplyr::filter(type == "Pandoravirus", expression_pattern == "Broad") %>% .$gene_id
pandora_Aza_expressed_ids <- active_aza %>% dplyr::filter(type == "Pandoravirus", expression_pattern %in% c("Aza_only","Aza_most") ) %>% .$gene_id

pfam_numbers %>% filter(gene_id %in% pandora_expressed_ids) %>% inner_join(., active_aza)
reactivated_domains <- pfam_numbers %>% filter(gene_id %in% pandora_Aza_expressed_ids) %>% inner_join(., active_aza) %>%
  .$Pfam %>% table(.) %>% data.frame() %>% arrange(desc(Freq)) %>% dplyr::rename(Domain = quote(.))
reactivated_domains$Domain <- factor(reactivated_domains$Domain, levels = rev(reactivated_domains$Domain))

gg_pandora_pfam_reactivated <- ggplot(reactivated_domains, aes(x = Domain, y = Freq)) + 
  stat_summary(geom = "bar", fun = mean, position = "dodge") +
  xlab("Pfam domain") +
  ylab("Domain containing genes") +
  ggtitle("Aza Reactivated GEVE genes") +
  theme_classic() + coord_flip()

ggsave(gg_pandora_pfam_reactivated, 
       filename = "PFAMs_reactivated_pandoravirus.pdf",
       height = 9, width = 7)
##################################################################
########## DEseq genes

DEseq_Aza_vs_DMSO <-fread("A_par_drug_treatment_DE_genes_table.tsv")

DEseq_Aza_vs_DMSO <- DEseq_Aza_vs_DMSO %>% filter(!gene_id %in% c(adinto_ids, plavaka_ids, te_genes_ids, pandora_ids))

upregulated_ids <- DEseq_Aza_vs_DMSO %>% filter(DMSO_vs_Aza_padj < 0.01, DMSO_vs_Aza_log2FC > 0) %>% .$gene_id 
downregulated_ids <- DEseq_Aza_vs_DMSO %>% filter(DMSO_vs_Aza_padj < 0.01, DMSO_vs_Aza_log2FC < 0) %>% .$gene_id 

upregulated_GOs <- get_GOs(geneNames = upregulated_ids, GOlist = geneID2GO)
downregulated_GOs <- get_GOs(geneNames = downregulated_ids, GOlist = geneID2GO)

upregulated_GOs_slimmed <- upregulated_GOs %>% filter(GO.ID %in% slimGOs)
downregulated_GOs_slimmed <- downregulated_GOs %>% filter(GO.ID %in% slimGOs)

write.table(x = upregulated_GOs, file = "UpregulatedAza.GOs.tsv", 
            quote = FALSE, sep = "\t",
            row.names = FALSE, col.names = TRUE)
write.table(x = downregulated_GOs, file = "DownregulatedAza.GOs.tsv", 
            quote = FALSE, sep = "\t",
            row.names = FALSE, col.names = TRUE)


gg_upregulated_GOs <- plot_GOs(df = upregulated_GOs, name = "Upregulated 5-Aza genes")
gg_downregulated_GOs <- plot_GOs(df = downregulated_GOs, name = "Downregulated 5-Aza genes")
ggsave(gg_upregulated_GOs, 
       filename = "GOs_Apar_upregulatedAza_GOs.pdf",
       height = 9, width = 7)
ggsave(gg_downregulated_GOs, 
       filename = "GOs_Apar_downregulatedAza_GOs.pdf",
       height = 9, width = 7)

upregulated_pfams <- pfamEnrichments(ids = upregulated_ids) %>% arrange(fishp)
upregulated_pfams$padj <- p.adjust(upregulated_pfams$fishp, method = "bonferroni")
upregulated_pfams <- upregulated_pfams %>% filter(padj < 0.05)
upregulated_pfams$Pfam <- factor(upregulated_pfams$Pfam, levels= rev(c(upregulated_pfams$Pfam)))
gg_upregulated_pfam <- ggplot(upregulated_pfams, aes(x = Pfam, y = -log10(fishp))) + 
  stat_summary(geom = "bar", fun = mean, position = "dodge") +
  xlab("Pfam domain") +
  ylab("-log10(p-value)") +
  ggtitle("5-Aza upregulated domain enrichment") +
  theme_classic() + coord_flip()

downregulated_pfams <- pfamEnrichments(ids = downregulated_ids) %>% arrange(fishp)
downregulated_pfams$padj <- p.adjust(downregulated_pfams$fishp, method = "bonferroni")
downregulated_pfams <- downregulated_pfams %>% filter(padj < 0.05)
downregulated_pfams$Pfam <- factor(downregulated_pfams$Pfam, levels= rev(c(downregulated_pfams$Pfam)))
gg_downregulated_pfam <- ggplot(downregulated_pfams, aes(x = Pfam, y = -log10(fishp))) + 
  stat_summary(geom = "bar", fun = mean, position = "dodge") +
  xlab("Pfam domain") +
  ylab("-log10(p-value)") +
  ggtitle("5-Aza downregulated domain enrichment") +
  theme_classic() + coord_flip()



ggsave(gg_upregulated_pfam, 
       filename = "PFAMs_upregulated_Aza.pdf",
       height = 9, width = 7)
ggsave(gg_downregulated_pfam, 
       filename = "PFAMs_downregulated_Aza.pdf",
       height = 9, width = 7)


eggnog <- fread("../Eggnog.annotations.tsv", fill = TRUE)
eggnog$query <- str_remove(eggnog$query, pattern = "\\.1")

# response to host defenses
eggnog %>% filter(grepl(GOs,pattern = "0052200"), query %in% upregulated_ids )
# response to host
eggnog %>% filter(grepl(GOs,pattern = "0075136"), query %in% upregulated_ids )
# response to other organism
eggnog %>% filter(grepl(GOs,pattern = "0051707"), query %in% upregulated_ids )
# cgutub based cuticle
eggnog %>% filter(grepl(GOs,pattern = "0018990"), query %in% upregulated_ids )


argonautes <- rbind(pfam_numbers %>% filter(gene_id %in% upregulated_ids, Pfam == "Piwi"),
pfam_numbers %>% filter(gene_id %in% upregulated_ids, Pfam == "RdRP"))



piwi_to_plot <-left_join(argonautes, tpm) %>% dplyr::select(Pfam, gene_id, AZA1,AZA2,AZA3,DMSO1,DMSO2,DMSO3) %>%
  filter( AZA1+AZA2+AZA3 > 3 ) %>% melt() %>% mutate(variable2 = ifelse(grepl(variable, pattern = "AZA"),"AZA","DMSO" ))

piwi_to_plot$variable2 <- factor(piwi_to_plot$variable2, levels = c("DMSO","AZA"))
  
gg_piwi <- ggplot(piwi_to_plot, aes(x = gene_id, y = value, fill = variable2)) +
  geom_boxplot(aes(fill=variable2)) +
  geom_point(position=position_dodge(width=0.75),aes(group=variable2)) +
  theme_bw() +
  ylab("TPM") + xlab("") +
  facet_grid(.~Pfam, scales="free_x")

ggsave(gg_piwi, 
       filename = "PiwiRdRP_upregulated_Aza.pdf",
       height = 4, width = 4)



######################################
#### Pandoravirus orthogroups

GVorthogroups <- fread("~/Dropbox/Amoebidium_genome/Chromosome_scale/ViralGenomeAnalysis/Orthogroups.GeneCount.GiantVirus.tsv")

GVorthogroups$presentIn <- rowSums(GVorthogroups > 0) -2
GVorthogroups$singlecopy <- rowSums(GVorthogroups == 1) 

GVorthogroups %>% arrange(desc(presentIn)) %>% head()
GVorthogroups %>% arrange(desc(singlecopy)) %>% head()

library(UpSetR)

movies <- read.csv(system.file("extdata", "movies.csv", package = "UpSetR"), 
                   header = T, sep = ";")

 
GVorthogroups_true <- 1*(GVorthogroups > 0) %>% data.frame() %>% dplyr::select(-Total,-presentIn, -singlecopy)
GVorthogroups_true$V1 <- GVorthogroups$V1
GVorthogroups_true <- GVorthogroups_true[GVorthogroups$presentIn > 2,]

pdf(file="~/Dropbox/Amoebidium_genome/Chromosome_scale/ViralGenomeAnalysis/UpsetPandoravirus.pdf") # or other device
upset(GVorthogroups_true, nsets = 40, point.size = 2, line.size = 1, order.by = "freq",
      mainbar.y.label = "Shared orthogroups", sets.x.label = "Orthogroup #",
      mb.ratio = c(0.55, 0.45),)
dev.off()

upset(GVorthogroups_true, mb.ratio = c(0.55, 0.45), order.by = "freq")


###########################
### What are the most common conserved GV orthogroups?

GVorthogroups_list <- fread("~/Dropbox/Amoebidium_genome/Chromosome_scale/ViralGenomeAnalysis/Orthogroups.GiantVirus.txt", sep = ":", header = F)
GVorthogroups %>% arrange(desc(presentIn)) %>% head(n = 17)

OG_ids <- GVorthogroups_list %>% filter(V1 == "OG0000005") %>% 
  .$V2 %>% str_split(., pattern= " ") %>% 
  str_remove_all(., pattern = "GV\\d+_chr\\d+_") %>% 
  str_remove_all(., pattern = "GV\\d+_unplaced\\d+_") %>%
  str_remove_all(., pattern = "\\.1")

eggnog %>% filter(query %in% OG_ids)

######################################
#### Timecourse

library(ggrepel)
library(ggfortify)
library(DESeq2)
tpm_mtx %>% rownames()


allZero <- rowSums(tpm_mtx==0)==ncol(tpm_mtx)
gene_counts <- tpm_mtx[!allZero,]
rownames(gene_counts ) <- rownames(tpm_mtx)[!allZero] %>% as.character()
gene_counts <- round(gene_counts, digits = 0)


d.correlation <- as.dist(1 - cor(gene_counts,method=c("spearman")))
fit <- hclust(d.correlation, method="complete")
pdf("RNAseq_clustering_byDEgenes.pdf", height = 5)
plot(fit) # display dendogram
dev.off()

gene_counts <- gene_counts[,1:8] %>% data.frame()
rownames(gene_counts ) <- rownames(tpm_mtx)[!allZero]
colnames(gene_counts ) <- colnames(tpm_mtx)[1:8]

d.correlation <- as.dist(1 - cor(gene_counts,method=c("spearman")))
fit <- hclust(d.correlation, method="complete")
plot(fit) # display dendogram


#coldata <- data.frame(sample=colnames(gene_counts)) %>%
#  mutate(condition = ifelse(grepl("AZA", sample), "AZA",
# ifelse(grepl("DMSO", sample), "DMSO",
coldata <- data.frame(sample=colnames(gene_counts)) %>%
  mutate(condition = ifelse(grepl("5h", sample), "5h",
                                         ifelse(grepl("14h", sample), "14h",
                                                ifelse(grepl("20h", sample), "20h",
                                                       ifelse(grepl("33h", sample), "33h","Broad"))))) 

dds_all <- DESeqDataSetFromMatrix(countData = gene_counts,
                                  colData = coldata,
                                  design = ~ condition)
dds_all <- estimateSizeFactors(dds_all)


dds_all <- dds_all[rowSums(counts(dds_all)) >= 10,]
dds_all <- dds_all[rowSums(counts(dds_all)==0) <= ncol(counts(dds_all))*0.5,]

#### PCA
tpm_mtx_dev <- tpm_mtx[,1:8]

pr <- prcomp(t(tpm_mtx_dev),center = T,scale = F)

#original Sam's function, where "group" is a separate object
pc1 <- (summary(pr)$importance[2, 1] * 100) %>% round(digits = 2)
pc2 <- (summary(pr)$importance[2, 2] * 100) %>% round(digits = 2)
group_col <- c(rep("black", times = 2),rep("red", times = 2), rep("green", times = 2), rep("darkgreen", times = 2))
group <- colnames(counts_mRNA)
pp1 <- autoplot(pr, alpha=0) + 
  geom_point(aes(color=factor(group)),
             alpha=0.8, size=4, fill=group_col) +
  geom_text_repel(aes(label=group),point.padding = unit(0.75, "lines")) +
  scale_color_manual(values=group_col) +
  xlab(str_c("PC1 (", pc1, "%)")) +
  ylab(str_c("PC2 (", pc2, "%)")) +
  theme_bw() + ggtitle(label = "PCA DMRs mCG/CG") +
  theme(legend.position="none")

ggsave(pp1, filename = "PCA_RNA-seq.pdf",
       height = 7, width = 7)

###################
# 1st comparison: noDox vs Dox


get_diff_genes <- function(cond1, cond2, df_deq = dds_all, fdrlevel.de = 0.05){
  dds <- df_deq
  
  dds <- dds[, dds$condition == cond1 | 
               dds$condition == cond2 ]
  dds$condition <- droplevels(dds$condition)
  
  #filter low expressed genes (at least 10 counts in total, or more than 0 in at least 50% of the samples)
  dds <- dds[rowSums(counts(dds)) >= 10,]
  dds_good <- dds[rowSums(counts(dds)==0) <= ncol(counts(dds))*0.5,]
  
  dds_good <- estimateSizeFactors(dds_good)
  
  rowData(dds_good)$control_expr <- rowMeans(counts(dds_good, normalize = TRUE)[,dds_good$condition == cond1])
  rowData(dds_good)$condition_expr <- rowMeans(counts(dds_good, normalize = TRUE)[,dds_good$condition == cond2])
  
  #######
  #do the test
  dds_good <- DESeq(dds_good)
  res_good <- results(dds_good)
  
  res_good$control_expr <- rowData(dds_good)$control_expr
  res_good$condition_expr <- rowData(dds_good)$condition_expr
  
  dds_good <- dds_good[!is.na(res_good$padj)]
  res_good <- res_good %>% na.omit()
  
  res_good <- res_good[res_good$padj < fdrlevel.de,]
  res_good$status <- ifelse(res_good$control_expr < res_good$condition_expr,"upregulated","downregulated")
  res_good$comparison <- paste0(cond1,"_vs_",cond2)
  
  return(res_good)
  
}


First_deseq <- get_diff_genes(cond1 = "5h", cond2 = "14h")
First_deseq$status %>% table()
Second_deseq <- get_diff_genes(cond1 = "14h", cond2 = "20h")
Second_deseq$status %>% table()
Third_deseq <- get_diff_genes(cond1 = "20h", cond2 = "33h")
Third_deseq$status %>% table()
Fourth_deseq <- get_diff_genes(cond1 = "33h", cond2 = "5h")
Fourth_deseq$status %>% table()

de_genes <- unique(c(rownames(First_deseq) ,rownames(Second_deseq), rownames(Third_deseq), rownames(Fourth_deseq)))

first_down_GOs <- get_GOs(geneNames = rownames(First_deseq)[First_deseq$status == "downregulated"], GOlist = geneID2GO)
first_up_GOs <- get_GOs(geneNames = rownames(First_deseq)[First_deseq$status == "upregulated"], GOlist = geneID2GO)

plot_GOs(df = first_down_GOs, name = "5h vs 14h downregulated")
plot_GOs(df = first_up_GOs, name = "5h vs 14h upregulated")

third_down_GOs <- get_GOs(geneNames = rownames(Third_deseq)[First_deseq$status == "downregulated"], GOlist = geneID2GO)
third_up_GOs <- get_GOs(geneNames = rownames(Third_deseq)[First_deseq$status == "upregulated"], GOlist = geneID2GO)

plot_GOs(df = third_down_GOs, name = "20h vs 33h downregulated")
plot_GOs(df = third_up_GOs, name = "20h vs 33h upregulated")

eggnog <- fread("../Eggnog.annotations.tsv", fill = TRUE)
eggnog$query <- str_remove(eggnog$query, pattern = "\\.1")

# striated muscle cell differentiation
eggnog %>% filter(grepl(GOs,pattern = "0051146"), query %in% rownames(Third_deseq)[First_deseq$status == "upregulated"] )
# it is a Myosin MYH9
eggnog %>% filter(grepl(GOs,pattern = "0140253"), query %in% rownames(Third_deseq)[First_deseq$status == "upregulated"] )



library(Mfuzz)

tpm_dev <- tpm_mtx_dev %>% data.frame()
rownames(tpm_dev) <- rownames(tpm_mtx)
colnames(tpm_dev) <- colnames(tpm_mtx_dev)
tpm_dev <- tpm_dev[,c("5h1","5h2","14h1","14h2","20h1","20h2","33h1","33h2")]

tpm_dev_mean <- tpm_dev %>% mutate(first = (`5h1`+`5h2`)/2,
                                   second = (`14h1`+`14h2`)/2,
                                   third = (`20h1`+`20h2`)/2,
                                   fourth = (`33h1`+`33h2`)/2) %>%
  dplyr::select(first, second, third, fourth)
  
scale_df <- function(df){
  df_scaled <- as.matrix(df) %>%
    ExpressionSet() %>% standardise() %>% exprs() %>% data.frame()
  colnames(df_scaled) <- colnames(df)
  return(df_scaled)
}

scaled_tpm <- scale_df(tpm_dev_mean)

scaled_tpm <- scaled_tpm[complete.cases(scaled_tpm), ] %>% as.matrix()

scaled_tpm_var <- scaled_tpm[de_genes,]


kmeans_ordered_clusters <- function(df, knum = 5){
  set.seed(20)
  df <- kmeans(df, centers = knum, nstart = 20)
  df <- df$cluster %>% as.data.frame() %>% dplyr::rename(cluster = ".") 
  df$loci <- rownames(df)
  df <- arrange(df, cluster)
  rownames(df) <- df$loci
  return(df)
}

df_clusters <- kmeans_ordered_clusters( df = scaled_tpm_var, knum = 5 )

#merge cluster IDs with mCG values

add_cluster_id_to_DF_to_ggplot <- function(df,cluster_id_dic){
  sample_count <- ncol(df)
  df$gene_id <- rownames(df)
  df$cluster <- cluster_id_dic[df$gene_id,]$cluster %>% as.factor()
  
  df <- melt(df)
  counts <- table(df$cluster)/sample_count
  counts <- as.data.frame(counts) %>% dplyr::rename(cluster = Var1) %>% melt()
  counts$value <- paste('Cluster ',counts$cluster,' (n = ', counts$value,')', sep = "")
  df$cluster2 <- factor(df$cluster, labels = as.character(counts$value))
  
  return(df)
}

df_clusters <- add_cluster_id_to_DF_to_ggplot(df = data.frame(scaled_tpm_var), cluster_id_dic = df_clusters )

#plot boxplot by DMR cluster

boxplot_DMRs <- ggplot(df_clusters, aes(y=value, x=variable )) +
  #geom_violin(alpha=0.1) +
  geom_boxplot(width=0.8, notch = TRUE, 
               outlier.size = 1, outlier.shape = 1, outlier.stroke = NA) +
  ylab("standarised TPM") + theme_bw() +
  xlab("stages") +
  theme( axis.text.x = element_text(angle = 45, hjust = 1)) + facet_grid(.~cluster2)

boxplot_DMRs <- boxplot_DMRs + facet_grid(. ~ cluster2) + ggtitle(label = "DSS DMRs ∆0.4 by kmeans clusters")

cluster1_GOs <- get_GOs(geneNames = df_clusters$gene_id[df_clusters$cluster == 1] %>% unique(), GOlist = geneID2GO)
cluster2_GOs <- get_GOs(geneNames = df_clusters$gene_id[df_clusters$cluster == 2] %>% unique(), GOlist = geneID2GO)
cluster3_GOs <- get_GOs(geneNames = df_clusters$gene_id[df_clusters$cluster == 3] %>% unique(), GOlist = geneID2GO)
cluster4_GOs <- get_GOs(geneNames = df_clusters$gene_id[df_clusters$cluster == 4] %>% unique(), GOlist = geneID2GO)
cluster5_GOs <- get_GOs(geneNames = df_clusters$gene_id[df_clusters$cluster == 5] %>% unique(), GOlist = geneID2GO)

plot_GOs(df = cluster1_GOs, name = "Cluster 1 = Up in 33h")
plot_GOs(df = cluster2_GOs, name = "Cluster 2 = Up in 14h to 33h")
plot_GOs(df = cluster3_GOs, name = "Cluster 3 = Downwards trend")
plot_GOs(df = cluster4_GOs, name = "Cluster 4 = Upwards trend")
plot_GOs(df = cluster5_GOs, name = "Cluster 5 = High 5h")

#profilin
tpm_dev_mean["APARG00000007788",] %>% boxplot()
#cofilin
tpm_dev_mean["APARG00000021960",] %>% boxplot()
#septins
tpm_dev_mean[c("APARG00000007402",
               "APARG00000009117",
               "APARG00000009250",
               "APARG00000010351",
               "APARG00000010751",
               "APARG00000011291",
               "APARG00000015424",
               "APARG00000022020"),] %>% boxplot()
# Talin
tpm_dev_mean[c("APARG00000012711",
               "APARG00000020228"),] %>% boxplot()

omaya_genes <- fread("~/Downloads/OmayaGenes.tsv", header = F, col.names = c("Name","SarcID"))
orthologues <- fread("~/Downloads/Amoebidium_parasiticum__v__Sphaeroforma_arctica.tsv", header = T, col.names = c("OG","gene_id","SarcID"))
inner_join(orthologues, omaya_genes) %>% mutate(ID= str_remove(gene_id, pattern = "\\.1")) %>% 
  dplyr::select(ID,Name)

omaya_apar <- fread("~/Downloads/Omaya_Apar.tsv")

Omaya_to_plot <- tpm_dev_mean %>% mutate(ID = rownames(.)) %>% inner_join(.,omaya_apar) %>% 
  melt() %>% mutate(hours = ifelse(variable == "first", "5h",ifelse(variable == "second","14h",ifelse(variable == "third", "20h","33h"))))
Omaya_to_plot$hours <- factor(Omaya_to_plot$hours, levels = c("5h","14h","20h","33h"))

library(cowplot)
ggMyosin <- ggplot(Omaya_to_plot %>% filter(grepl(Name, pattern = "Myosin")), aes(x = hours , y = value, group = ID, colour = Name)) + geom_line() + geom_point() + theme_bw() + ylab("TPM") + xlab("")
ggRho <- ggplot(Omaya_to_plot %>% filter(grepl(Name, pattern = "Rho")), aes(x = hours , y = value, group = ID, colour = Name)) + geom_line() + geom_point() + theme_bw() + ylab("TPM") + xlab("")
ggProfilin <- ggplot(Omaya_to_plot %>% filter(grepl(Name, pattern = "Profilin|Cofilin")), aes(x = hours , y = value, group = ID, colour = Name)) + geom_line() + geom_point() + theme_bw() + ylab("TPM") + xlab("")
ggKinesin <- ggplot(Omaya_to_plot %>% filter(grepl(Name, pattern = "Kinesin")), aes(x = hours , y = value, group = ID, colour = Name)) + geom_line() + geom_point() + theme_bw() + ylab("TPM") + xlab("")
ggArpFormin <- ggplot(Omaya_to_plot %>% filter(grepl(Name, pattern = "Arp|Formin")), aes(x = hours , y = value, group = ID, colour = Name)) + geom_line() + geom_point() + theme_bw() + ylab("TPM") + xlab("")
ggTubulin <- ggplot(Omaya_to_plot %>% filter(grepl(Name, pattern = "tubulin")), aes(x = hours , y = value, group = ID, colour = Name)) + geom_line() + geom_point() + theme_bw() + ylab("TPM") + xlab("")
ggSeptin <- ggplot(Omaya_to_plot %>% filter(grepl(Name, pattern = "Septin")), aes(x = hours , y = value, group = ID, colour = Name)) + geom_line() + geom_point() + theme_bw() + ylab("TPM") + xlab("")
ggAdhesion <- ggplot(Omaya_to_plot %>% filter(grepl(Name, pattern = "Pinch|Parvin|Paxilin|Talin|Vinculin")), aes(x = hours , y = value, group = ID, colour = Name)) + geom_line() + geom_point() + theme_bw() + ylab("TPM") + xlab("")

ggsave(ggMyosin, filename = "Myosin_plot.pdf", height =3 , width = 4)
ggsave(ggProfilin, filename = "ProCofilinplot.pdf", height = 3 , width = 4)
ggsave(ggArpFormin, filename = "ArpFormin_plot.pdf", height = 3 , width = 4)
ggsave(ggTubulin, filename = "Tubulin_plot.pdf", height = 3 , width = 4)
ggsave(ggSeptin, filename = "Septin_plot.pdf", height = 3 , width = 4)
ggsave(ggAdhesion, filename = "Adheesion_plot.pdf", height = 3 , width = 4)

################################################
### Select longer TEs according to family size


repeat_modeler <- fread("~/Dropbox/Amoebidium_genome/Chromosome_scale/RepeatSmallRNA/RepeatModelar_EarlGrey_combined_library.fasta.fai") %>% .[,c(1,2)]
colnames(repeat_modeler) <- c("TEid","length")

inserts <- fread("~/Dropbox/Amoebidium_genome/Chromosome_scale/RepeatSmallRNA/Apar_genome.fasta.RM.500.noGiant.bed", col.names = c("chr","start","stop","class","id","strand"))

inserts  <-inserts %>% mutate(TEid = paste0(id,"#",class), width = stop - start, locus = paste0(chr,":",start,"-",stop))

inserts <- inner_join(inserts,repeat_modeler ) %>% mutate(perc_original = 100*width/length )

inserts_70 <- inserts %>% filter(perc_original > 70) 

write.table(inserts_70, file = "~/Dropbox/Amoebidium_genome/Chromosome_scale/RepeatSmallRNA/Longer_repeats70perc_1000bp.bed", quote = F, row.names = F, col.names = F, sep = "\t")

hist(inserts_70$width, breaks = 1000, xlim = c(0,5000))


################################################
### Plot genome size / repeats / global mCG

library(openxlsx)
library(ggplot2)
library(dplyr)

comparative_analysis <- read.xlsx("~/Dropbox/CGmaps/Book_of_mCGs.xlsx",sheet = 2) %>%
  dplyr::rename(repeat_perc = `Repeat.%`, genome_size = Genome_Size, species = Name ) %>%
  mutate(repeats = genome_size * repeat_perc / 100) %>% 
  mutate(rest = genome_size - repeats)
  

comparative_analysis$species <- factor(comparative_analysis$species, levels = rev(comparative_analysis$species))

df <- comparative_analysis %>% dplyr::select(species, repeats, rest) %>% melt()
df$variable <- factor(df$variable, levels = c("rest","repeats"))  

  # Create a ggplot barplot
gg_genome_sizes <- ggplot(df, aes(x = species, y = value, colour = variable, fill = variable)) +
    geom_bar(stat = "identity") + scale_y_continuous(limits = c(0, 400), breaks = seq(0, 400, by = 50)) +
    theme_bw() + coord_flip()

gg_mCG_species <- ggplot(comparative_analysis, aes(x = species, y = mCG_level)) + geom_bar(stat = "identity") + 
  ylim(c(0,100)) + theme_bw() + coord_flip()


gg_both_plots <- plot_grid(gg_genome_sizes, gg_mCG_species, ncol = 2)

ggsave(gg_both_plots, filename = "~/Dropbox/Amoebidium_genome/Chromosome_scale/ViralGenomeAnalysis/GenomeSize_GlobalmCGplot.pdf", height = 3.5, width = 8)


################################################
### Plot Holozoan BUSCO and Repeat Landscapes

library(openxlsx)
library(ggplot2)
library(dplyr)

holozoanTEs <- read.xlsx("~/Dropbox/Amoebidium_genome/Chromosome_scale/ViralGenomeAnalysis/Repeats.xlsx")

holozoanTEs$Species <- factor(holozoanTEs$Species, levels = c("Amoebidium","Ichthyophonus","Sphaeroforma","Creolimax","Pirum","Abeoforma","Chromosphaera",
                                                              "Corallochytrium","Capsaspora","Salpingoeca")) 

holozoanTEs <- melt(holozoanTEs) %>% filter(variable != "Total.repeats")

gg_TE_distribution <- ggplot(holozoanTEs, aes(x = Species, y = value, fill = variable)) +
  geom_bar(stat = "identity") +
  theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + facet_grid(variable~.)

ggsave(gg_TE_distribution, filename = "~/Dropbox/Amoebidium_genome/Chromosome_scale/ViralGenomeAnalysis/HolozoanTEs.pdf", height = 3.5, width = 4)


holozoanBUSCOs <- read.xlsx("~/Dropbox/Amoebidium_genome/Chromosome_scale/ViralGenomeAnalysis/Repeats.xlsx", sheet = 2)

holozoanBUSCOs$Species <- factor(holozoanBUSCOs$Species, levels = rev(c("Amoebidium","Ichthyophonus","Sphaeroforma","Creolimax","Pirum","Abeoforma","Chromosphaera",
                                                              "Corallochytrium","Capsaspora","Salpingoeca","Monosiga")) )
holozoanBUSCOs <- melt(holozoanBUSCOs)

gg_BUSCO_distribution <- ggplot(holozoanBUSCOs, aes(x = Species, y = value, fill = variable)) +
  geom_bar(position="fill", stat="identity") +
  theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  coord_flip()
ggsave(gg_BUSCO_distribution, filename = "~/Dropbox/Amoebidium_genome/Chromosome_scale/ViralGenomeAnalysis/HolozoanBUSCOs.pdf", height = 3.5, width = 4)


################################################
### Plot insertion length

GiantRepeatsBed <- fread("~/Dropbox/Amoebidium_genome/Chromosome_scale/GiantRepeats.bed", col.names = c("chr","start","end","Type")) %>%
  mutate(length = (end - start)/ 1000)

GiantRepeatsBed$Type <- factor(GiantRepeatsBed$Type, levels = c("Pandoravirus","Adintovirus","Plavaka"))

gg_InsertionDistribution <- ggplot(GiantRepeatsBed %>% filter(grepl(chr, pattern = "chr" )), aes(x = Type, y = length)) + geom_jitter(size = 0.5) + geom_boxplot(outlier.shape = NA) +
  scale_y_sqrt(breaks = c(10, 25, 50, 100, 200,300,400)) + theme_bw()
ggsave(gg_InsertionDistribution, filename = "~/Dropbox/Amoebidium_genome/Chromosome_scale/ViralGenomeAnalysis/Insertions.pdf", height = 2.5, width = 1.5)

GiantRepeatsBed %>% filter(grepl(chr, pattern = "chr" )) %>% .$Type %>% table()

fai <- fread("~/Dropbox/Amoebidium_genome/Chromosome_scale/Apar_genome.fasta.fai", col.names = c("chr","length","cumulative","stuff","stuff2"))  

chr_cumulative <- fai$cumulative[fai$chr == "chr18"]

insertion_to_plot <- GiantRepeatsBed %>% group_by(Type) %>% summarize(100*sum(end - start)/chr_cumulative) %>% mutate(genome = "Aapp")

colnames(insertion_to_plot) <- c("Type","GenomePerc","genome")

perc_giant_gg <- ggplot(insertion_to_plot, aes(x = genome,  y = GenomePerc, fill = Type)) +
  geom_bar(position="stack", stat="identity") + ylim(c(0,20)) + theme_classic()

gene_total <-  fread("~/Dropbox/Amoebidium_genome/Chromosome_scale/ViralGenomeAnalysis/Amoebidium.chr.v3.annot.gene.bed") %>% nrow()

giant_gene_num <- data.frame( type = c("pandora","adinto","plavaka","rest"), 
            gene_num = c(length(pandora_ids), length(adinto_ids), length(plavaka_ids), 
                         gene_total - length(pandora_ids) - length(adinto_ids) - length(plavaka_ids))) %>%
  mutate(perc = 100 * gene_num / gene_total, genome = "App")

giant_gene_num$type <- factor(giant_gene_num$type, levels = c("pandora","adinto","plavaka","rest"))

perc_giantGenes_gg <- ggplot(giant_gene_num, aes(x = genome,  y = perc, fill = type)) +
  geom_bar(position="stack", stat="identity") + ylim(c(0,20)) + theme_classic()

gg_giant_combo <- plot_grid(perc_giant_gg, perc_giantGenes_gg, ncol = 2)

ggsave(gg_giant_combo, filename = "~/Dropbox/Amoebidium_genome/Chromosome_scale/ViralGenomeAnalysis/GiantGenomePercs.pdf", height = 2.5, width = 4)


################################################
### Compare Amoebidium genomes

library(openxlsx)
library(ggplot2)
library(dplyr)

IsolatesDF <- read.xlsx("~/Dropbox/Amoebidium_genome/Chromosome_scale/ViralGenomeAnalysis/Repeats.xlsx", sheet = 4)

IsolatesDF$Isolate <- factor(IsolatesDF$Isolate, levels = IsolatesDF$Isolate) 

Isolate_repeats <- IsolatesDF %>% dplyr::select(Isolate, LTR, LINE, DNA, Repeats_perc, Assembly_size) %>%
  mutate(rest_repeats = Repeats_perc - (LTR + LINE + DNA) ) %>% 
  mutate(rest_genome = Assembly_size - Repeats_perc*Assembly_size/100,
         LTR = Assembly_size*(LTR/100),LINE = Assembly_size*(LINE/100),DNA = Assembly_size*(DNA/100),
         rest_repeats = Assembly_size*(rest_repeats/100)) %>% 
  dplyr::select(-Assembly_size, -Repeats_perc) %>% melt()

Isolate_repeats$variable <- factor(Isolate_repeats$variable, levels = rev(c("LTR","LINE","DNA","rest_repeats","rest_genome")))

gg_TE_isolates <- ggplot(Isolate_repeats, aes(x = Isolate, y = value, colour = variable, fill = variable)) +
  geom_bar(stat = "identity") +theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggsave(gg_TE_isolates, filename = "~/Dropbox/Amoebidium_genome/Chromosome_scale/ViralGenomeAnalysis/IsolatesTEs.pdf", height = 3, width = 2.6)


IsolateBUSCOs <- read.xlsx("~/Dropbox/Amoebidium_genome/Chromosome_scale/ViralGenomeAnalysis/Repeats.xlsx", sheet = 2) %>%
  dplyr::filter(Species %in% c("Amoebidium","Isolate_9181","Isolate_9257")) %>% dplyr::select(-X6)

IsolateBUSCOs <- melt(IsolateBUSCOs)

gg_BUSCO_isolates <- ggplot(IsolateBUSCOs, aes(x = Species, y = value, fill = variable)) +
  geom_bar(position="fill", stat="identity") +
  theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 
ggsave(gg_BUSCO_isolates, filename = "~/Dropbox/Amoebidium_genome/Chromosome_scale/ViralGenomeAnalysis/IsolatesBUSCOs.pdf", height = 3.5, width = 2.6)


################################################
### Compare Amoebidium isolates mCG %

library(openxlsx)
library(ggplot2)
library(dplyr)

Isolates_mCG <- read.xlsx("~/Dropbox/Amoebidium_genome/Chromosome_scale/ViralGenomeAnalysis/Repeats.xlsx", sheet = 3)

Isolates_mCG <- Isolates_mCG %>% dplyr::filter(Context %in% c("WmCGW","CGC","mCG"))

gg_mCG_isolates <- ggplot(Isolates_mCG, aes(x = Isolate, y = Percentage, colour = Context )) + geom_point(size = 3) + 
  ylim(c(0,100)) + theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggsave(gg_mCG_isolates, filename = "~/Dropbox/Amoebidium_genome/Chromosome_scale/ViralGenomeAnalysis/Isolates_mCGs.pdf", height = 3.5, width = 2.6)
