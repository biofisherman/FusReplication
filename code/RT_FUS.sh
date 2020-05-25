#!/bin/bash
## Replication timing analysis
## directory
# ulimit -n 100000
DIR=/home/weiyan/Desktop/ubuntu_share/RT_transcription
FQDIR=/home/weiyan/Desktop/ubuntu_share/RT_transcription/data/fq
RAW=/home/weiyan/Desktop/ubuntu_share/RT_transcription/data/fq/raw/R1

### trimmed reads directory
mkdir -p $FQDIR/trimmed/R1
mkdir -p $FQDIR/trimmed/R1/report
TRIM=$FQDIR/trimmed/R1

#reference genome.
REF=$DIR/ref/Homo_sapiens.GRCh38.dna.primary_assembly.fa

### index directory
# mkdir -p $DIR/index
# mkdir -p $DIR/index/hg38
INX=$DIR/index/hg38

### aligned reads directory
mkdir -p $DIR/data/aligned/R1
ALG=$DIR/data/aligned/R1

### Signal files directory
mkdir -p $DIR/data/signal/R1
SIG=$DIR/data/signal/R1

### setting for signal
SAM=sample.txt
BIN=20000
SLE=100000
#Generate the root names
parallel -j 1 echo {1}_{2} ::: U2OS Clone110 FUSClone110 ::: G1 S > names.txt
NAMEs=names.txt

# ## Raw reads QC by fastp(version 0.20.1)
cat $NAMEs | parallel --jobs 6 "fastp -i $RAW/{}.fastq.gz -o $TRIM/{}.fastq -w 30 -j $TRIM/report/{}".json" -h $TRIM/report/{}".html""

## Map trimmed reads with bowtie2

### build index
bowtie2-build $REF $INX --threads 30

### mapping
cat $NAMEs | parallel --jobs 1 "bowtie2 -x $INX --threads 30 -U $TRIM/{}.fastq | samtools view -bS -q 20 | samtools sort - > $ALG/{}.bam"

### index
cat $NAMEs | parallel --jobs 6 "samtools index $ALG/{}.bam"

## calculate log2 ratios between S/G1

## download hg38 blacklist "ENCFF419RSJ.bed" on encode website

bedtools merge -i ENCFF419RSJ.bed > hg38Blacklist.bed

cat $SAM | parallel --jobs 1 " \
bamCompare -b1 $ALG/{}_S.bam -b2 $ALG/{}_G1.bam \
           --operation log2 \
           --skipZeroOverZero \
           --scaleFactorsMethod None \
           --binSize $BIN \
           --blackListFileName $DIR/ref/hg38Blacklist.bed \
           --effectiveGenomeSize 2805636331 \
           --normalizeUsing CPM \
           --ignoreDuplicates \
           --minMappingQuality 20 \
           -p 30 \
           --smoothLength $SLE \
           --outFileFormat bigwig \
           -o $SIG/{}_$BIN.bw
"