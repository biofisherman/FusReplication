#!/bin/bash
## enrich plots
WIN=500000 # or 5000
R=LRD_FAD.bed ## or other domain related bed files used in figures
## generate matrix
computeMatrix reference-point \
 -S $trdir/W_2.bw \
    $trdir/C_2.bw \
    $trdir/F_2.bw \
 -R $R \
 --referencePoint center \
 -a $WIN -b $WIN \
 --skipZeros \
 --missingDataAsZero \
 -p 30 \
 -out $dir/results/matrix_transcription"_"$R"_"$WIN.gz 
## plot heatmap
plotHeatmap \
 -m $dir/results/matrix_transcription"_"$R"_"$WIN.gz \
 -out $dir/figures/RD/transcription"_"$R"_"$WIN.pdf \
 --colorMap YlGnBu \
 --heatmapHeight 20 \
 --plotTitle 'Transcription Signal in $R'
# plot profile
plotProfile \
     -m $dir/results/matrix_transcription"_"$R"_"$WIN.gz \
     -out $dir/figures/RD/profile_transcription"_"$R"_"$WIN.pdf \
     --perGroup \
     --plotTitle 'Transcription Signal in $R'
