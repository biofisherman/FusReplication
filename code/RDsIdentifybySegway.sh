# Segway environment:
 
conda create --name segway python=2 segway=2.0.2 genomedata=1.4.1 optplus=0.1 optbuild=0.1 autolog=0.1

source activate segway

# loading data:

## R1
genomedata-load --sizes --sequence hg38.chrom.sizes.txt --track U2OS_R2=U2OS_RT_R1-X_Loess_smoothing.bedGraph U2OS_RT.genomedata

genomedata-load --sizes --sequence hg38.chrom.sizes.txt --track Clone110_R2=Clone110_RT_R1-X_Loess_smoothing.bedGraph Clone110_RT.genomedata

genomedata-load --sizes --sequence hg38.chrom.sizes.txt --track FUSClone110_R2=FUSClone110_RT_R1-X_Loess_smoothing.bedGraph FUSClone110_RT.genomedata

## R2
genomedata-load --sizes --sequence hg38.chrom.sizes.txt --track U2OS_R2=U2OS_RT_R2-X_Loess_smoothing.bedGraph U2OS_RT1.genomedata

genomedata-load --sizes --sequence hg38.chrom.sizes.txt --track Clone110_R2=Clone110_RT_R2-X_Loess_smoothing.bedGraph Clone110_RT1.genomedata

genomedata-load --sizes --sequence hg38.chrom.sizes.txt --track FUSClone110_R2=FUSClone110_RT_R2-X_Loess_smoothing.bedGraph FUSClone110_RT1.genomedata

# Train and identify segments

## for U2OS_3lables

### train
segway --resolution=10000 --num-instances=25 --minibatch-fraction=0.01 --num-labels=3 --max-train-rounds=20 --exclude-coords=hg38.blacklist.bed train U2OS_RT.genomedata U2OS_RT1.genomedata train_results
#running locally 9531: emt24.19.bundle.train_results.2175879ec4df11e8a862ab7fb22a3465 ()

### Identify
segway --exclude-coords="hg38.blacklist.bed" --bigBed=annotate_results/segway.layered.bb identify U2OS_RT.genomedata U2OS_RT1.genomedata train_results annotate_results

## for Clone110_3lables
### train
segway --resolution=10000 --num-instances=25 --minibatch-fraction=0.01 --num-labels=3 --max-train-rounds=20 --exclude-coords=hg38.blacklist.bed train Clone110_RT.genomedata Clone110_RT1.genomedata train_results

### Identify
segway --exclude-coords="hg38.blacklist.bed" --bigBed=annotate_results/segway.layered.bb identify Clone110_RT.genomedata Clone110_RT1.genomedata train_results annotate_results

## For FUSClone110_3_Lables
### train
segway --resolution=10000 --num-instances=25 --minibatch-fraction=0.01 --num-labels=3 --max-train-rounds=20 --exclude-coords=hg38.blacklist.bed train FUSClone110_RT.genomedata FUSClone110_RT1.genomedata train_results

### Identify
segway --exclude-coords="hg38.blacklist.bed" --bigBed=annotate_results/segway.layered.bb identify FUSClone110_RT.genomedata FUSClone110_RT1.genomedata train_results annotate_results

## Splite bed files to ERD, LRD and MRD:
### U2OS
awk '4=="0"' U2OS_segway.bed > U2OS_ERD.bed
awk '4=="1"' U2OS_segway.bed > U2OS_LRD.bed
awk '4=="2"' U2OS_segway.bed > U2OS_MRD.bed

### Clone110
awk '4=="0"' Clone110_segway.bed > Clone110_ERD.bed
awk '4=="1"' Clone110_segway.bed > Clone110_LRD.bed
awk '4=="2"' Clone110_segway.bed > Clone110_MRD.bed

### FUSClone110
awk '4=="0"' FUSClone110_segway.bed > FUSClone110_MRD.bed
awk '4=="1"' FUSClone110_segway.bed > FUSClone110_ERD.bed
awk '4=="2"' FUSClone110_segway.bed > FUSClone110_LRD.bed

#### Change labels of FUSClone110 to make it consistent with U2OS and Clone110

awk -F "\t" -v OFS="\t" '{4="0"; 9="27,158,119" ; print 0}' FUSClone110_ERD.bed > FUSClone110_ERD_1.bed
awk -F "\t" -v OFS="\t" '{4="1"; 9="217,95,2" ; print 0}' FUSClone110_LRD.bed > FUSClone110_LRD_1.bed
awk -F "\t" -v OFS="\t" '{4="2"; 9="117,112,179" ; print 0}' FUSClone110_MRD.bed > FUSClone110_MRD_1.bed


# ideantify FADs

## FAD-ERDs(ERD_lost.bed)
### 1. Find overlapped ERD domains in U2OS and FUSClone110
bedtools intersect -a U2OS_ERD.bed -b FUSClone110_ERD_1.bed -sorted > U2OS_FUSClone110_overlap.bed
### 2. find lost domains in U2OS_FUSClone110_overlap.bed
bedtools subtract -a U2OS_FUSClone110_overlap.bed -b Clone110_ERD.bed > ERD_lost.bed

#For FAD-LRD:(LRD_lost.bed)
### 1. Find overlapped LRD domains in U2OS and FUSClone110
bedtools intersect -a U2OS_LRD.bed -b FUSClone110_LRD_1.bed -sorted > U2OS_FUSClone110_LRD_overlap.bed
### 2. find lost domains in U2OS_FUSClone110_overlap.bed
bedtools subtract -a U2OS_FUSClone110_LRD_overlap.bed -b Clone110_LRD.bed > LRD_lost.bed

#For FAD-MRD:(MRD_lost.bed)
### 1. Find overlapped MRD domains in U2OS and FUSClone110
bedtools intersect -a U2OS_MRD.bed -b FUSClone110_MRD_1.bed -sorted > U2OS_FUSClone110_MRD_overlap.bed
### 2. find lost domains in U2OS_FUSClone110_overlap.bed
bedtools subtract -a U2OS_FUSClone110_MRD_overlap.bed -b Clone110_MRD.bed > MRD_lost.bed


# extract gene annotation information in FADs

## 1. annotation informatin downloaded from UCSC table(all Gencode v28)

## 2. Add regulation region to gene.(3kb on both directions)
$ bedtools slop -i hg38_genecodeV28_gene.bed -g hg38.txt -b 3000 > hg38_genecodV28_gene_3k    
## 3. Find genes in RTDs:
bedtools intersect -a hg38_genecodV28_gene_3k.bed -b ERD_lost.bed  > ERD_lost_gene.bed
bedtools intersect -a hg38_genecodV28_gene_3k.bed -b MRD_lost.bed  > MRD_lost_gene.bed
bedtools intersect -a hg38_genecodV28_gene_3k.bed -b LRD_lost.bed  > LRD_lost_gene.bed
	
# cut column
cut -f 4 ERD_lost_gene_3k.bed > ERD_lost_genes.txt
cut -d. -f 1 ERD_lost_genes.txt > ERD_lost_genes_Ref.txt

cut -f 4 MRD_lost_gene_3k.bed > MRD_lost_genes.txt
cut -d. -f 1 MRD_lost_genes.txt > MRD_lost_genes_Ref.txt

cut -f 4 LRD_lost_gene_3k.bed > LRD_lost_genes.txt
cut -d. -f 1 LRD_lost_genes.txt > LRD_lost_genes_Ref.txt


