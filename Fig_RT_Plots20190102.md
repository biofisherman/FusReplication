RT\_Plots\_Jia\_01022019
================
Weiyan
1/2/2019

Note: all the following data sets were generated by Lu Liu through R
anlaysis.

# 1\. Data input

``` r
setwd("/Users/weiyanjia/SynologyDrive/Randal\ S.\ Tibbetts/FUS\ project/Bioinfo_analysis/FUS_Replication\ timing/Replication\ timing_08312018/01_02_normalization_smoothing/plots_Jia/01022019")
library(RColorBrewer)
RT_Loess <- read.delim("all_sample_r1_r2_smooth.txt", header = TRUE, sep ="\t", stringsAsFactors = FALSE)
RT_Loess[, c(2:3)] <- sapply(RT_Loess[, c(2:3)], as.integer)
head(RT_Loess)
```

    ##    chr  start    end Clone110_R1 FUSClone110_R1  U2OS_R1 Clone110_R2
    ## 1 chr1 800000 820000    1.003375       1.208018 1.406038   0.9013547
    ## 2 chr1 820000 840000    1.022251       1.196786 1.414605   0.9081274
    ## 3 chr1 840000 860000    1.040475       1.186075 1.422639   0.9143926
    ## 4 chr1 860000 880000    1.058048       1.175888 1.430135   0.9201442
    ## 5 chr1 880000 900000    1.074967       1.166228 1.437088   0.9253764
    ## 6 chr1 900000 920000    1.091232       1.157098 1.443494   0.9300832
    ##   FUSClone110_R2   U2OS_R2
    ## 1      0.9879291 0.9895395
    ## 2      1.0126643 1.0077907
    ## 3      1.0363650 1.0255111
    ## 4      1.0590256 1.0427003
    ## 5      1.0806404 1.0593578
    ## 6      1.1012037 1.0754831

``` r
chr_list<-c()
for (i in 1:22){
  chr_list[i] <- paste0("chr",i,sep="")
  print(chr_list)
}
```

    ## [1] "chr1"
    ## [1] "chr1" "chr2"
    ## [1] "chr1" "chr2" "chr3"
    ## [1] "chr1" "chr2" "chr3" "chr4"
    ## [1] "chr1" "chr2" "chr3" "chr4" "chr5"
    ## [1] "chr1" "chr2" "chr3" "chr4" "chr5" "chr6"
    ## [1] "chr1" "chr2" "chr3" "chr4" "chr5" "chr6" "chr7"
    ## [1] "chr1" "chr2" "chr3" "chr4" "chr5" "chr6" "chr7" "chr8"
    ## [1] "chr1" "chr2" "chr3" "chr4" "chr5" "chr6" "chr7" "chr8" "chr9"
    ##  [1] "chr1"  "chr2"  "chr3"  "chr4"  "chr5"  "chr6"  "chr7"  "chr8"  "chr9" 
    ## [10] "chr10"
    ##  [1] "chr1"  "chr2"  "chr3"  "chr4"  "chr5"  "chr6"  "chr7"  "chr8"  "chr9" 
    ## [10] "chr10" "chr11"
    ##  [1] "chr1"  "chr2"  "chr3"  "chr4"  "chr5"  "chr6"  "chr7"  "chr8"  "chr9" 
    ## [10] "chr10" "chr11" "chr12"
    ##  [1] "chr1"  "chr2"  "chr3"  "chr4"  "chr5"  "chr6"  "chr7"  "chr8"  "chr9" 
    ## [10] "chr10" "chr11" "chr12" "chr13"
    ##  [1] "chr1"  "chr2"  "chr3"  "chr4"  "chr5"  "chr6"  "chr7"  "chr8"  "chr9" 
    ## [10] "chr10" "chr11" "chr12" "chr13" "chr14"
    ##  [1] "chr1"  "chr2"  "chr3"  "chr4"  "chr5"  "chr6"  "chr7"  "chr8"  "chr9" 
    ## [10] "chr10" "chr11" "chr12" "chr13" "chr14" "chr15"
    ##  [1] "chr1"  "chr2"  "chr3"  "chr4"  "chr5"  "chr6"  "chr7"  "chr8"  "chr9" 
    ## [10] "chr10" "chr11" "chr12" "chr13" "chr14" "chr15" "chr16"
    ##  [1] "chr1"  "chr2"  "chr3"  "chr4"  "chr5"  "chr6"  "chr7"  "chr8"  "chr9" 
    ## [10] "chr10" "chr11" "chr12" "chr13" "chr14" "chr15" "chr16" "chr17"
    ##  [1] "chr1"  "chr2"  "chr3"  "chr4"  "chr5"  "chr6"  "chr7"  "chr8"  "chr9" 
    ## [10] "chr10" "chr11" "chr12" "chr13" "chr14" "chr15" "chr16" "chr17" "chr18"
    ##  [1] "chr1"  "chr2"  "chr3"  "chr4"  "chr5"  "chr6"  "chr7"  "chr8"  "chr9" 
    ## [10] "chr10" "chr11" "chr12" "chr13" "chr14" "chr15" "chr16" "chr17" "chr18"
    ## [19] "chr19"
    ##  [1] "chr1"  "chr2"  "chr3"  "chr4"  "chr5"  "chr6"  "chr7"  "chr8"  "chr9" 
    ## [10] "chr10" "chr11" "chr12" "chr13" "chr14" "chr15" "chr16" "chr17" "chr18"
    ## [19] "chr19" "chr20"
    ##  [1] "chr1"  "chr2"  "chr3"  "chr4"  "chr5"  "chr6"  "chr7"  "chr8"  "chr9" 
    ## [10] "chr10" "chr11" "chr12" "chr13" "chr14" "chr15" "chr16" "chr17" "chr18"
    ## [19] "chr19" "chr20" "chr21"
    ##  [1] "chr1"  "chr2"  "chr3"  "chr4"  "chr5"  "chr6"  "chr7"  "chr8"  "chr9" 
    ## [10] "chr10" "chr11" "chr12" "chr13" "chr14" "chr15" "chr16" "chr17" "chr18"
    ## [19] "chr19" "chr20" "chr21" "chr22"

``` r
chr_list<-c(chr_list,"chrX")

RT_Loess$chr <- factor(RT_Loess$chr, levels= chr_list)

# write.table(RT_Loess,"RT_Loess_x.txt", sep="\t", row.names=FALSE, quote=FALSE, col.names = TRUE)
```

# 2\. boxplots of timing values

``` r
library(reshape2)
dat.m <- melt(RT_Loess, id.vars='chr', measure.vars =c("U2OS_R1", "Clone110_R1","FUSClone110_R1","U2OS_R2","Clone110_R2","FUSClone110_R2") )

library(ggpubr)
```

    ## Loading required package: ggplot2

``` r
ggboxplot(dat.m, x="variable", y="value",
         ylab = "Replication Timing", 
         ylim=c(-4,4),
         xlab ="",
         legend.title ="",
         x.text.angle = 45,
         color = "variable", palette = "aaas"
         )+
  border(color ="black", size =1, linetype = NULL)
```

![](Fig_RT_Plots20190102_files/figure-gfm/unnamed-chunk-2-1.png)<!-- -->

``` r
palette = "aaas"
```

# 3\. replication timing profile

## 3.1 RT-Chr18-R1

``` r
sub_chr <- subset(RT_Loess, RT_Loess$chr == "chr18")
head(sub_chr)
```

    ##         chr  start    end Clone110_R1 FUSClone110_R1     U2OS_R1 Clone110_R2
    ## 50006 chr18  80000 100000  -0.9158458     -0.4703737 -0.20717267  -1.3131937
    ## 50007 chr18 100000 120000  -0.8365058     -0.3986828 -0.14072280  -1.2315270
    ## 50008 chr18 120000 140000  -0.7589561     -0.3293648 -0.07661734  -1.1511733
    ## 50009 chr18 140000 160000  -0.6831967     -0.2624348 -0.01486341  -1.0721301
    ## 50010 chr18 160000 180000  -0.6092274     -0.1979082  0.04453186  -0.9943948
    ## 50011 chr18 180000 200000  -0.5370481     -0.1358000  0.10156135  -0.9179649
    ##       FUSClone110_R2     U2OS_R2
    ## 50006     -0.8781096 -0.06621107
    ## 50007     -0.7956762 -0.01103857
    ## 50008     -0.7154465  0.04214985
    ## 50009     -0.6374241  0.09334109
    ## 50010     -0.5616127  0.14252203
    ## 50011     -0.4880161  0.18967957

``` r
class(sub_chr)
```

    ## [1] "data.frame"

``` r
ggline(sub_chr, x="start", y=c("U2OS_R1", "Clone110_R1","FUSClone110_R1"),
       merge = TRUE,
         ylab = "Replication Timing", 
         ylim=c(-2,2),
         xlim = c(2e+07,5e+07),
         xlab ="Coordinate",
         plot_type = "l",
          palette = "aaas") +
  geom_hline(yintercept =0, linetype =2) +
  xscale("none", .format = TRUE)+
  border(color ="black", size =1, linetype = NULL)
```

![](Fig_RT_Plots20190102_files/figure-gfm/unnamed-chunk-3-1.png)<!-- -->

\#\#3.2
RT-Chr18-R2

``` r
ggline(sub_chr, x="start", y=c("U2OS_R2", "Clone110_R2","FUSClone110_R2"),
       merge = TRUE,
         ylab = "Replication Timing", 
         ylim=c(-2,2),
         xlim = c(2e+07,5e+07),
         xlab ="Coordinate",
         plot_type = "l",
          palette = "aaas") +
  geom_hline(yintercept =0, linetype =2) +
  xscale("none", .format = TRUE)+
  border(color ="black", size =1, linetype = NULL)
```

![](Fig_RT_Plots20190102_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

# 4\. RT Density distribution plot

## 4.1 Density distribution of R1\_all\_chr

``` r
ggdensity(RT_Loess, x = c("U2OS_R1", "Clone110_R1","FUSClone110_R1"), 
          y="..density..",
          xlim = c(-3,3),
          xlab = "Replication Timing",
          merge = TRUE,
          add = "median",                  # Add median line. 
         palette = "aaas"
         )+
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0))
```

![](Fig_RT_Plots20190102_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

## 4.2 Density distribution of R2\_all\_chr

``` r
ggdensity(RT_Loess, x = c("U2OS_R2", "Clone110_R2","FUSClone110_R2"), 
          y="..density..",
          xlim = c(-3,3),
          xlab = "Replication Timing",
          merge = TRUE,
          add = "median",                  # Add median line. 
         palette = "aaas"
         )+
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0))
```

![](Fig_RT_Plots20190102_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->

## 4.3 Density distribution of R2\_single\_chr2

``` r
sub_chr <- subset(RT_Loess, RT_Loess$chr == "chr2")
ggdensity(sub_chr, x = c("U2OS_R2", "Clone110_R2","FUSClone110_R2"), 
          y="..density..",
          xlim = c(-3,3),
          xlab = "Replication Timing",
          merge = TRUE,
          add = "median",                  # Add median line. 
          # rug = TRUE , 
         palette = "aaas"
         )+
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0))
```

![](Fig_RT_Plots20190102_files/figure-gfm/unnamed-chunk-7-1.png)<!-- -->

## 4.4 Density distribution of R2\_single\_chr5

``` r
sub_chr <- subset(RT_Loess, RT_Loess$chr == "chr5")
ggdensity(sub_chr, x = c("U2OS_R2", "Clone110_R2","FUSClone110_R2"), 
          y="..density..",
          xlim = c(-3,3),
          xlab = "Replication Timing",
          merge = TRUE,
          add = "median",                  # Add median line. 
         palette = "aaas"
         )+
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0))
```

![](Fig_RT_Plots20190102_files/figure-gfm/unnamed-chunk-8-1.png)<!-- -->

## 4.5 Density distribution of R2\_single\_chr20

``` r
sub_chr <- subset(RT_Loess, RT_Loess$chr == "chr20")
ggdensity(sub_chr, x = c("U2OS_R2", "Clone110_R2","FUSClone110_R2"), 
          y="..density..",
          xlim = c(-3,3),
          xlab = "Replication Timing",
          merge = TRUE,
          add = "median",                  # Add median line. 
          # rug = TRUE , 
         palette = "aaas"
         )+
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0))
```

![](Fig_RT_Plots20190102_files/figure-gfm/unnamed-chunk-9-1.png)<!-- -->

## 4.3 Density distribution of R1\_single\_chr2

``` r
sub_chr <- subset(RT_Loess, RT_Loess$chr == "chr2")
ggdensity(sub_chr, x = c("U2OS_R1", "Clone110_R1","FUSClone110_R1"), 
          y="..density..",
          xlim = c(-3,3),
          xlab = "Replication Timing",
          merge = TRUE,
          add = "median",                  # Add median line. 
          # rug = TRUE , 
         palette = "aaas"
         )+
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0))
```

![](Fig_RT_Plots20190102_files/figure-gfm/unnamed-chunk-10-1.png)<!-- -->

## 4.4 Density distribution of R1\_single\_chr5

``` r
sub_chr <- subset(RT_Loess, RT_Loess$chr == "chr5")
ggdensity(sub_chr, x = c("U2OS_R1", "Clone110_R1","FUSClone110_R1"), 
          y="..density..",
          xlim = c(-3,3),
          xlab = "Replication Timing",
          merge = TRUE,
          add = "median",                  # Add median line. 
         palette = "aaas"
         )+
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0))
```

![](Fig_RT_Plots20190102_files/figure-gfm/unnamed-chunk-11-1.png)<!-- -->

## 4.5 Density distribution of R1\_single\_chr20

``` r
sub_chr <- subset(RT_Loess, RT_Loess$chr == "chr20")
ggdensity(sub_chr, x = c("U2OS_R1", "Clone110_R1","FUSClone110_R1"), 
          y="..density..",
          xlim = c(-3,3),
          xlab = "Replication Timing",
          merge = TRUE,
          add = "median",                  # Add median line. 
          # rug = TRUE , 
         palette = "aaas"
         )+
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0))
```

![](Fig_RT_Plots20190102_files/figure-gfm/unnamed-chunk-12-1.png)<!-- -->

# 5\. correlection-heatmap plot

## 5.1 heatmap

``` r
library(ComplexHeatmap)
```

    ## Loading required package: grid

    ## ========================================
    ## ComplexHeatmap version 2.0.0
    ## Bioconductor page: http://bioconductor.org/packages/ComplexHeatmap/
    ## Github page: https://github.com/jokergoo/ComplexHeatmap
    ## Documentation: http://jokergoo.github.io/ComplexHeatmap-reference
    ## 
    ## If you use it in published research, please cite:
    ## Gu, Z. Complex heatmaps reveal patterns and correlations in multidimensional 
    ##   genomic data. Bioinformatics 2016.
    ## ========================================

``` r
library(circlize)
```

    ## ========================================
    ## circlize version 0.4.10
    ## CRAN page: https://cran.r-project.org/package=circlize
    ## Github page: https://github.com/jokergoo/circlize
    ## Documentation: https://jokergoo.github.io/circlize_book/book/
    ## 
    ## If you use it in published research, please cite:
    ## Gu, Z. circlize implements and enhances circular visualization
    ##   in R. Bioinformatics 2014.
    ## 
    ## This message can be suppressed by:
    ##   suppressPackageStartupMessages(library(circlize))
    ## ========================================

``` r
corr_data <- RT_Loess[4:ncol(RT_Loess)]
corr_matrix <- round(cor(corr_data, method = "pearson"), 2)
col_fun = colorRamp2(c(0.6, 0.8, 1), c("#4DAF4A", "#FFD92F", "#E41A1C"))

Heatmap(corr_matrix, name = "", col = col_fun,
        row_names_side = "left",
        row_order = c("U2OS_R1","U2OS_R2","FUSClone110_R1","FUSClone110_R2",
                      "Clone110_R1","Clone110_R2"),
        column_order = c("U2OS_R1","U2OS_R2","FUSClone110_R1","FUSClone110_R2",
                         "Clone110_R1","Clone110_R2"),
    # column_order = c("Clone110_S_R2","Clone110_S_R1","FUSClone110_S_R2","FUSClone110_S_R1",
                       #"U2OS_S_R2","U2OS_S_R1"),
    # column_order = rev(order(colnames(corr_matrix))),
      clustering_distance_rows = "pearson",
      clustering_distance_columns = "pearson",
      cluster_rows = FALSE,
      cluster_columns = FALSE,
      cell_fun = function(j, i, x, y, width, height, fill) {
      grid.text(sprintf("%.2f", corr_matrix[i, j]), x, y, gp = gpar(fontsize = 10))
      },
      heatmap_legend_param = list(
      at = c(0.6,0.7, 0.8,0.9, 1.0),
     legend_height = unit(4.5, "inch"),
     title_position = "topleft")
    )
```

![](Fig_RT_Plots20190102_files/figure-gfm/unnamed-chunk-13-1.png)<!-- -->

## 5.2 PCA

``` r
library(ggfortify)
library(tidyverse)
```

    ## ── Attaching packages ───────────────────────────────────────────────────── tidyverse 1.3.0 ──

    ## ✓ tibble  3.0.1     ✓ dplyr   1.0.0
    ## ✓ tidyr   1.1.0     ✓ stringr 1.4.0
    ## ✓ readr   1.3.1     ✓ forcats 0.5.0
    ## ✓ purrr   0.3.4

    ## ── Conflicts ──────────────────────────────────────────────────────── tidyverse_conflicts() ──
    ## x dplyr::filter() masks stats::filter()
    ## x dplyr::lag()    masks stats::lag()

``` r
library(export)
```

``` r
pca_data <- RT_Loess[4:ncol(RT_Loess)]
# head(pca_data)

pca_data2<- tibble::rowid_to_column(corr_data, "ID")

pca_data2$ID <-as.factor(pca_data2$ID)
# head(pca_data2)

Tran_pca_data <- gather(pca_data2,sample,RT,Clone110_R1:U2OS_R2)
Tran_pca_data_S <- spread(Tran_pca_data,ID, RT)
# head(Tran_pca_data_S)

pca_res <- prcomp(Tran_pca_data_S[2:ncol(Tran_pca_data_S)], scale. = FALSE)
autoplot(pca_res, data = Tran_pca_data_S, colour = 'sample')
```

    ## Warning: `select_()` is deprecated as of dplyr 0.7.0.
    ## Please use `select()` instead.
    ## This warning is displayed once every 8 hours.
    ## Call `lifecycle::last_warnings()` to see where this warning was generated.

![](Fig_RT_Plots20190102_files/figure-gfm/unnamed-chunk-15-1.png)<!-- -->

``` r
graph2pdf(file="/Users/weiyanjia/SynologyDrive/Randal\ S.\ Tibbetts/FUS\ project/Bioinfo_analysis/FUS_Replication\ timing/Replication\ timing_08312018/01_02_normalization_smoothing/plots_Jia/01022019/PCA_RT.pdf", width=8, aspectr=sqrt(2),font = "Arial",bg = "transparent")
```

    ## Exported graph as /Users/weiyanjia/SynologyDrive/Randal S. Tibbetts/FUS project/Bioinfo_analysis/FUS_Replication timing/Replication timing_08312018/01_02_normalization_smoothing/plots_Jia/01022019/PCA_RT.pdf

``` r
Tran_pca_data_S_T<- Tran_pca_data_S%>%
                    column_to_rownames("sample")
pca_res <- prcomp(Tran_pca_data_S_T, scale. = FALSE)
var_explained <- pca_res$sdev^2/sum(pca_res$sdev^2)
var_explained[1:6]
```

    ## [1] 8.009698e-01 7.954648e-02 6.457223e-02 3.359288e-02 2.131863e-02
    ## [6] 2.133411e-30

``` r
names(pca_res)
```

    ## [1] "sdev"     "rotation" "center"   "scale"    "x"

``` r
pca_res$x %>% 
  as.data.frame %>%
  rownames_to_column("sample") %>%
  ggplot(aes(x=PC1,y=PC2)) + geom_point(aes(color=sample),size=4) +
  theme_bw(base_size=20) + 
  labs(x=paste0("PC1: ",round(var_explained[1]*100,1),"%"),
       y=paste0("PC2: ",round(var_explained[2]*100,1),"%")) +
  theme(legend.position="top")
```

![](Fig_RT_Plots20190102_files/figure-gfm/unnamed-chunk-16-1.png)<!-- -->

``` r
graph2pdf(file="/Users/weiyanjia/SynologyDrive/Randal\ S.\ Tibbetts/FUS\ project/Bioinfo_analysis/FUS_Replication\ timing/Replication\ timing_08312018/01_02_normalization_smoothing/plots_Jia/01022019/PCA_RT2.pdf", width=8, aspectr=sqrt(2),font = "Arial",bg = "transparent")
```

    ## Exported graph as /Users/weiyanjia/SynologyDrive/Randal S. Tibbetts/FUS project/Bioinfo_analysis/FUS_Replication timing/Replication timing_08312018/01_02_normalization_smoothing/plots_Jia/01022019/PCA_RT2.pdf

``` r
pca_res$x %>% 
  as.data.frame %>%
  rownames_to_column("sample") %>%
  ggplot(aes(x=PC1,y=PC2, label=sample, color=sample)) +
  geom_label(aes(fill = sample), colour = "white", fontface = "bold")+
  theme_bw(base_size=20) +
  labs(x=paste0("PC1: ",round(var_explained[1]*100,1),"%"),
       y=paste0("PC2: ",round(var_explained[2]*100,1),"%"))+
  theme(legend.position="top")
```

![](Fig_RT_Plots20190102_files/figure-gfm/unnamed-chunk-17-1.png)<!-- -->

``` r
graph2pdf(file="/Users/weiyanjia/SynologyDrive/Randal\ S.\ Tibbetts/FUS\ project/Bioinfo_analysis/FUS_Replication\ timing/Replication\ timing_08312018/01_02_normalization_smoothing/plots_Jia/01022019/PCA_RT3.pdf", width=8, aspectr=sqrt(2),font = "Arial",bg = "transparent")
```

    ## Exported graph as /Users/weiyanjia/SynologyDrive/Randal S. Tibbetts/FUS project/Bioinfo_analysis/FUS_Replication timing/Replication timing_08312018/01_02_normalization_smoothing/plots_Jia/01022019/PCA_RT3.pdf

``` r
pca_res_scale <- prcomp(Tran_pca_data_S[2:ncol(Tran_pca_data_S)], scale. = TRUE)
# pca_res$x

autoplot(pca_res_scale, data = Tran_pca_data_S, colour = 'sample')
```

![](Fig_RT_Plots20190102_files/figure-gfm/unnamed-chunk-18-1.png)<!-- -->

``` r
library(M3C)
# pca_data
# s <- colnames(pca_data[1,])
# 
# celltypes<-c("Clone110_R1","FUSClone110_R1","U2OS_R1","Clone110_R2","FUSClone110_R2","U2OS_R2")
# tsne(pca_data,labels=as.factor(celltypes))
pca(pca_data,legendtextsize = 10,axistextsize = 10,dotsize=2)
```

![](Fig_RT_Plots20190102_files/figure-gfm/unnamed-chunk-19-1.png)<!-- -->

``` r
pca(pca_data,labels=as.factor(colnames(pca_data[1,])),legendtextsize = 10,axistextsize = 10,dotsize=2)
```

![](Fig_RT_Plots20190102_files/figure-gfm/unnamed-chunk-19-2.png)<!-- -->

# 6 Histogram plots of RT difference

``` r
RT_Loess$diff_U2OS_Clone110_R1 <- RT_Loess$Clone110_R1 - RT_Loess$U2OS_R1
RT_Loess$diff_U2OS_FUSClone110_R1 <- RT_Loess$FUSClone110_R1 - RT_Loess$U2OS_R1
RT_Loess$diff_U2OS_Clone110_R2 <- RT_Loess$Clone110_R2 - RT_Loess$U2OS_R2
RT_Loess$diff_U2OS_FUSClone110_R2 <- RT_Loess$FUSClone110_R2 - RT_Loess$U2OS_R2

head(RT_Loess)
```

    ##    chr  start    end Clone110_R1 FUSClone110_R1  U2OS_R1 Clone110_R2
    ## 1 chr1 800000 820000    1.003375       1.208018 1.406038   0.9013547
    ## 2 chr1 820000 840000    1.022251       1.196786 1.414605   0.9081274
    ## 3 chr1 840000 860000    1.040475       1.186075 1.422639   0.9143926
    ## 4 chr1 860000 880000    1.058048       1.175888 1.430135   0.9201442
    ## 5 chr1 880000 900000    1.074967       1.166228 1.437088   0.9253764
    ## 6 chr1 900000 920000    1.091232       1.157098 1.443494   0.9300832
    ##   FUSClone110_R2   U2OS_R2 diff_U2OS_Clone110_R1 diff_U2OS_FUSClone110_R1
    ## 1      0.9879291 0.9895395            -0.4026630               -0.1980204
    ## 2      1.0126643 1.0077907            -0.3923539               -0.2178189
    ## 3      1.0363650 1.0255111            -0.3821632               -0.2365636
    ## 4      1.0590256 1.0427003            -0.3720869               -0.2542465
    ## 5      1.0806404 1.0593578            -0.3621212               -0.2708598
    ## 6      1.1012037 1.0754831            -0.3522625               -0.2863955
    ##   diff_U2OS_Clone110_R2 diff_U2OS_FUSClone110_R2
    ## 1           -0.08818479             -0.001610420
    ## 2           -0.09966328              0.004873531
    ## 3           -0.11111856              0.010853844
    ## 4           -0.12255610              0.016325284
    ## 5           -0.13398139              0.021282613
    ## 6           -0.14539989              0.025720594

## 6.1 R1

``` r
gghistogram(RT_Loess, x = c("diff_U2OS_FUSClone110_R1", "diff_U2OS_Clone110_R1"), 
          y="..count..",
          position = "dodge",
          xlim = c(-3,3),
          bins = 50,
          xlab = "Replication Timing",
          color = ".x.", 
          fill = ".x.",
          merge = TRUE,
          add = "median",                # Add median line. 
         palette = "aaas"
         )+
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0))
```

![](Fig_RT_Plots20190102_files/figure-gfm/unnamed-chunk-21-1.png)<!-- -->
\#\# 6.1
R2

``` r
gghistogram(RT_Loess, x = c("diff_U2OS_FUSClone110_R2","diff_U2OS_Clone110_R2"), 
          y="..count..",
          position = "dodge",
          xlim = c(-3,3),
          bins = 100,
          xlab = "Replication Timing",
          color = ".x.", fill = ".x.",
          merge = TRUE,
          add = "median",                  # Add median line. 
         palette = "aaas"
         )+
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) 
```

![](Fig_RT_Plots20190102_files/figure-gfm/unnamed-chunk-22-1.png)<!-- -->

# 7\. Plot whole chromosomes RT data

## 7.1 facet\_grid of RT

``` r
p <- ggline(RT_Loess, x="start", y=c("U2OS_R2", "Clone110_R2","FUSClone110_R2"),
       merge = TRUE,
       size= 0.5,
         ylab = "Replication Timing", 
         ylim=c(-3,3),
         xlab ="Coordinate",
         plot_type = "l",
          palette = "aaas") +
  geom_hline(yintercept =0, linetype =2) +
  xscale("none", .format = TRUE)+
  border(color ="black", size =0.5, linetype = NULL)+
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0), breaks = c(-2,-1, 0, 1, 2))

p + facet_grid(rows = vars(chr), scales = "free", space = "free", margins = FALSE, shrink = TRUE)+
theme(panel.spacing = unit(0.1, "lines"),
      panel.border = element_rect(linetype ="solid", fill = NA)) 
```

![](Fig_RT_Plots20190102_files/figure-gfm/unnamed-chunk-23-1.png)<!-- -->

## 7.2 RT plot of R2

``` r
library(ggsci)
## palette from package ggsci

palette1 <- pal_aaas("default")(10)
# palette1

palette2 <- pal_npg("nrc")(10)
# palette2

palette3 <- pal_lancet("lanonc")(9)
# palette3

big_palette <- c(palette1,palette2,palette3)
# big_palette

big_palette_clean <- big_palette[!duplicated(big_palette)]

# big_palette_clean

p <- ggline(RT_Loess, x="start", y="U2OS_R2",
       merge = TRUE,
       size = 0.5,
         ylab = "Replication Timing", 
         ylim=c(-3,3),
         #xlim = c(0.0e+00,2.6e+08),
         xlab ="Chromosome",
         plot_type = "l",
         color = "chr",
         palette = big_palette_clean) +
  geom_hline(yintercept =0, linetype =2) +
  xscale("none", .format = TRUE)+
  border(color ="black", size =0.5, linetype = NULL)+
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0), breaks = c(-2,-1, 0, 1, 2)) +
  theme (
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
    )

p + facet_grid(cols = vars(chr), scales = "free", space = "free", drop = TRUE, margins = FALSE, shrink = TRUE, switch = "x") +
   theme(panel.spacing = unit(0, "lines")) 
```

![](Fig_RT_Plots20190102_files/figure-gfm/unnamed-chunk-24-1.png)<!-- -->

``` r
sessionInfo()
```

    ## R version 3.6.3 (2020-02-29)
    ## Platform: x86_64-apple-darwin15.6.0 (64-bit)
    ## Running under: macOS Catalina 10.15.5
    ## 
    ## Matrix products: default
    ## BLAS:   /Library/Frameworks/R.framework/Versions/3.6/Resources/lib/libRblas.0.dylib
    ## LAPACK: /Library/Frameworks/R.framework/Versions/3.6/Resources/lib/libRlapack.dylib
    ## 
    ## locale:
    ## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
    ## 
    ## attached base packages:
    ## [1] parallel  grid      stats     graphics  grDevices utils     datasets 
    ## [8] methods   base     
    ## 
    ## other attached packages:
    ##  [1] ggsci_2.9            M3C_1.6.0            Biobase_2.44.0      
    ##  [4] BiocGenerics_0.30.0  export_0.2.2.9001    forcats_0.5.0       
    ##  [7] stringr_1.4.0        dplyr_1.0.0          purrr_0.3.4         
    ## [10] readr_1.3.1          tidyr_1.1.0          tibble_3.0.1        
    ## [13] tidyverse_1.3.0      ggfortify_0.4.10     circlize_0.4.10     
    ## [16] ComplexHeatmap_2.0.0 ggpubr_0.4.0         ggplot2_3.3.2       
    ## [19] reshape2_1.4.4       RColorBrewer_1.1-2  
    ## 
    ## loaded via a namespace (and not attached):
    ##   [1] readxl_1.3.1            uuid_0.1-4              snow_0.4-3             
    ##   [4] backports_1.1.8         systemfonts_0.2.3       NMF_0.22.0             
    ##   [7] plyr_1.8.6              splines_3.6.3           sigclust_1.1.0         
    ##  [10] crosstalk_1.1.0.1       gridBase_0.4-7          digest_0.6.25          
    ##  [13] foreach_1.5.0           htmltools_0.5.0         matrixcalc_1.0-3       
    ##  [16] viridis_0.5.1           fansi_0.4.1             magrittr_1.5           
    ##  [19] cluster_2.1.0           doParallel_1.0.15       openxlsx_4.1.5         
    ##  [22] modelr_0.1.8            officer_0.3.12          colorspace_1.4-1       
    ##  [25] blob_1.2.1              rvest_0.3.5             haven_2.3.1            
    ##  [28] xfun_0.15               crayon_1.3.4            jsonlite_1.7.0         
    ##  [31] survival_3.2-3          iterators_1.0.12        glue_1.4.1             
    ##  [34] registry_0.5-1          rvg_0.2.5               gtable_0.3.0           
    ##  [37] webshot_0.5.2           GetoptLong_1.0.2        car_3.0-8              
    ##  [40] shape_1.4.4             abind_1.4-5             scales_1.1.1           
    ##  [43] DBI_1.1.0               rngtools_1.5            bibtex_0.4.2.2         
    ##  [46] rstatix_0.6.0           miniUI_0.1.1.1          Rcpp_1.0.5             
    ##  [49] viridisLite_0.3.0       xtable_1.8-4            clue_0.3-57            
    ##  [52] foreign_0.8-75          htmlwidgets_1.5.1       httr_1.4.1             
    ##  [55] ellipsis_0.3.1          pkgconfig_2.0.3         farver_2.0.3           
    ##  [58] dbplyr_1.4.4            tidyselect_1.1.0        labeling_0.3           
    ##  [61] rlang_0.4.6             manipulateWidget_0.10.1 later_1.1.0.1          
    ##  [64] munsell_0.5.0           cellranger_1.1.0        tools_3.6.3            
    ##  [67] cli_2.0.2               generics_0.0.2          broom_0.5.6            
    ##  [70] evaluate_0.14           fastmap_1.0.1           yaml_2.2.1             
    ##  [73] knitr_1.29              fs_1.4.2                zip_2.0.4              
    ##  [76] rgl_0.100.54            dendextend_1.13.4       nlme_3.1-148           
    ##  [79] mime_0.9                xml2_1.3.2              compiler_3.6.3         
    ##  [82] rstudioapi_0.11         curl_4.3                png_0.1-7              
    ##  [85] ggsignif_0.6.0          reprex_0.3.0            stringi_1.4.6          
    ##  [88] gdtools_0.2.2           stargazer_5.2.2         lattice_0.20-41        
    ##  [91] Matrix_1.2-18           vctrs_0.3.1             pillar_1.4.4           
    ##  [94] lifecycle_0.2.0         GlobalOptions_0.1.2     corpcor_1.6.9          
    ##  [97] data.table_1.12.8       flextable_0.5.10        httpuv_1.5.4           
    ## [100] R6_2.4.1                promises_1.1.1          gridExtra_2.3          
    ## [103] rio_0.5.16              codetools_0.2-16        assertthat_0.2.1       
    ## [106] pkgmaker_0.31.1         rjson_0.2.20            withr_2.2.0            
    ## [109] hms_0.5.3               doSNOW_1.0.18           rmarkdown_2.3          
    ## [112] carData_3.0-4           Rtsne_0.15              shiny_1.5.0            
    ## [115] lubridate_1.7.9         base64enc_0.1-3
