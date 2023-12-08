##### This script takes a CGmap and makes 3 bigwigs: 
##### one for mCG, one for mCH, and one for coverage on Cs

PATH=$PATH:/data/home/btx832/src/UCSC

# you will need the CGmap as first variable
# the name of the output file as second variable
# the name of the "genome index" as third variable

# the genome index you can obtain from "samtools faidx genome.fasta"
# and the resulting "genome.fasta.fai" would be the genome index

CGmap=$1
prefix=$2
chromFile=$3

### make bedgraph of CG methylation level (0-1)
#Correct for change from 1-based to 0-based
#and remove unplaced contigs
# Seperate file for each strand

### Methylation level

#### CG context
zcat "$CGmap" | awk -v OFS="\t" \
'$1 != "M77789.2" && $1 != "NC_001416.1" && $5 == "CG" {print $1,$3-1,$3,$6}' \
|sort -k1,1 -k2,2n > "$prefix".CG.level.bedgraph &&

#### CH context
zcat "$CGmap" | awk -v OFS="\t" \
'$1 != "M77789.2" && $1 != "NC_001416.1" && $4 ~/CH/ {print $1,$3-1,$3,$6}' \
|sort -k1,1 -k2,2n > "$prefix".CH.level.bedgraph &&

### Cytosine coverage
zcat "$CGmap" | awk -v OFS="\t" \
'{ if ( $1 != "M77789.2" && $1 != "NC_001416.1" ){ if ($2 == "G") {print $1,$3-1,$3,-$8}else if ($2 == "C") {print $1,$3-1,$3,$8}}}' \
|sort -k1,1 -k2,2n > "$prefix".C.coverage.bedgraph

bedSort "$prefix".CG.level.bedgraph "$prefix".CG.level.sorted.bedgraph
bedSort "$prefix".CH.level.bedgraph "$prefix".CH.level.sorted.bedgraph
bedSort "$prefix".C.coverage.bedgraph "$prefix".C.coverage.sorted.bedgraph

# bedgraph to bigWig
bedGraphToBigWig "$prefix".CG.level.sorted.bedgraph $chromFile "$prefix".CG.level.bigwig
bedGraphToBigWig "$prefix".CH.level.sorted.bedgraph $chromFile "$prefix".CH.level.bigwig
bedGraphToBigWig "$prefix".C.coverage.sorted.bedgraph $chromFile "$prefix".C.coverage.bigwig
