#!/bin/bash
# build FUS regulated genes annotation bed file
# FUS.csv: the FUS specific regulated genes from DESeq2
# batch gene list to transcripts list by biomart() Gene stable ID(s) to Transcript stable ID. The file was saved as "mart_export_fus.txt"
## the total genes are 708, two Pseudogene was removed.
# awk -F "\"*,\"*" 'NR>1{print $2}' FUS.csv > FUS.txt # $2: column from FUS.csv; NR>1: Record Number (NR) is greater than 1 to remove head.
awk -F "\t" 'NR>1{print $2}' mart_export_fus.txt > FUSTransID.txt
# awk -F "\t" 'NR>1{print $1}' mart_export_fus.txt |sort -u > FUSGeneID.txt
grep -F -f FUSTransID.txt hg38GencodeV32.bed > FUS.bed