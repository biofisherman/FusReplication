#!/bin/bash
## enrich plots for Z score normalized data
BIN=20000
WIN=5000 # or 500,000
R=FUS.bed ## or hg38GencodeV32.bed for all genes
computeMatrix reference-point \
 --referencePoint TSS \
 -S U2OS_RT_R2-X_z_score_normal.bedgraph.bw \
    Clone110_RT_R2-X_z_score_normal.bedgraph.bw \
    FUSClone110_RT_R2-X_z_score_normal.bedgraph.bw \
 -R $R \
 -a $WIN -b $WIN \
 --skipZeros \
 --missingDataAsZero \
 -p 30 \
 -out results/R2/matrix_RT_TSS_$R"_"$BIN"_WS_"$WIN.gz 
 # --outFileNameMatrix results/R2/matrix_RT_TSS_$R"_"$BIN"_WS_"$WIN.tab \
 # --outFileSortedRegions results/R2/regions_RT_TSS_$R"_"$BIN"_WS_"$WIN.bed
## plot heatmap
plotHeatmap \
 -m results/R2/matrix_RT_TSS_$R"_"$BIN"_WS_"$WIN.gz \
 -out figures/R2/RT_TSS_$R"_"$BIN"_WS_"$WIN.pdf \
 --colorMap YlGnBu \
 --heatmapHeight 20 \
 --plotTitle 'Replication timing Signal in TSS'
# plot profile
plotProfile \
     -m results/R2/matrix_RT_TSS_$R"_"$BIN"_WS_"$WIN.gz \
     -out figures/R2/Profile_RT_TSS_$R"_"$BIN"_WS_"$WIN.pdf \
     --perGroup \
     --plotTitle 'Replication timing Signal in TSS'