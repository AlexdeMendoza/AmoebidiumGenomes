#!/bin/bash
#$ -m bea
#$ -cwd
#$ -pe smp 1
#$ -l h_vmem=16G
#$ -j y
#$ -l h_rt=1:0:0

cd /data/SBCS-ademendoza/02-lukesarre/amoebidium_paper_trimmed/deeptools

module load anaconda3
source activate deeptools-env

cd CGmaps/

# . ~/bin/FromCGmap_to_BigWig_local.sh BS_cgc.CGmap.gz BS_cgc.CGmap.gz Apar_genome.fasta.fai
# . ~/bin/FromCGmap_to_BigWig_local.sh BS_non_cgc.CGmap.gz BS_non_cgc.CGmap.gz Apar_genome.fasta.fai

cd ../

computeMatrix scale-regions -out repeat_families_with_CGC_CG.dtmatrix -S ./CGmaps/BS_cgc.CGmap.gz.CG.level.bigwig -R ./BEDs/DNA.bed ./BEDs/LINE.bed ./BEDs/LTR.bed -a 3000 -b 3000 -m 3000 -bs 100 -p 1 --transcriptID "mRNA"

computeMatrix scale-regions -out repeat_families_without_CGC_CG.dtmatrix -S ./CGmaps/BS_non_cgc.CGmap.gz.CG.level.bigwig -R ./BEDs/DNA.bed ./BEDs/LINE.bed ./BEDs/LTR.bed -a 3000 -b 3000 -m 3000 -bs 100 -p 1 --transcriptID "mRNA"

plotHeatmap -m repeat_families_with_CGC_CG.dtmatrix -out repeat_families_with_CGC_CG.heatmap.pdf --colorMap Reds --regionsLabel "DNA" "LINE" "LTR" -max 1 --interpolationMethod nearest

plotHeatmap -m repeat_families_without_CGC_CG.dtmatrix -out repeat_families_without_CGC_CG.heatmap.pdf --colorMap Reds --regionsLabel "DNA" "LINE" "LTR" -max 1 --interpolationMethod nearest

computeMatrix scale-regions -out expression_genes_with_CGC_CG.dtmatrix -S ./CGmaps/BS_cgc.CGmap.gz.CG.level.bigwig -R ./BEDs/unexpressed_genes.bed ./BEDs/gene_expression_decile_3.bed ./BEDs/gene_expression_decile_6.bed ./BEDs/gene_expression_decile_10.bed -a 3000 -b 3000 -m 3000 -bs 100 -p 1 --transcriptID "mRNA"

computeMatrix scale-regions -out expression_genes_without_CGC_CG.dtmatrix -S ./CGmaps/BS_non_cgc.CGmap.gz.CG.level.bigwig -R ./BEDs/unexpressed_genes.bed ./BEDs/gene_expression_decile_3.bed ./BEDs/gene_expression_decile_6.bed ./BEDs/gene_expression_decile_10.bed -a 3000 -b 3000 -m 3000 -bs 100 -p 1 --transcriptID "mRNA"

plotHeatmap -m expression_genes_with_CGC_CG.dtmatrix -out expression_genes_with_CGC_CG.heatmap.pdf --colorMap Reds --regionsLabel "Unexpressed" "3rd" "6th" "10th" -max 100 --interpolationMethod nearest
plotHeatmap -m expression_genes_without_CGC_CG.dtmatrix -out expression_genes_without_CGC_CG.heatmap.pdf --colorMap Reds --regionsLabel "Unexpressed" "3rd" "6th" "10th" -max 100 --interpolationMethod nearest
