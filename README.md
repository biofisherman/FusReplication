# Fused in sarcoma regulates DNA replication timing and progression
## Scripts and resource data for figures:
* __Fig. 3 and Sup. Fig. 5C to E:__ [GSEA anlysis of FUS RNA-Seq](/Fig_GSEA_FUS.md)
  * [Rmd file](code/Fig_GSEA_FUS.Rmd).
  * [MSigDB for biological processing(c5.bp.v6.2.symbols.gmt)](data/).
  * [raw data of res_WTvsKO_shrunken_tb_entrez.csv](data/res_WTvsKO_shrunken_tb_entrez.csv).
  * [raw data of normalized_counts.csv](data/normalized_counts.csv).
  * [raw data of ann.sample.csv](data/ann.sample.csv).
  * [raw data of FUS_Seq_SSIV_All.csv](data/FUS_Seq_SSIV_All.csv).
* __Fig. 6A and Sup. Fig. 7:__ [Analysis of FUS MS/MS data](/fig_FUS_MS.md)
  * [Rmd file](code/fig_FUS_MS.Rmd).
  * [MSigDB for biological processing(c5.bp.v6.2.symbols.gmt)](data/).
  * [raw data of AllResultsOfQuantification_FUS.csv](data/AllResultsOfQuantification_FUS.csv).
  * [raw data of PSMs_FUS.csv](data/PSMs_FUS.csv).
* __Fig. 6E:__ [Quantification results of PLA_FUS](/Fig_PLA_FUS.md)
  * [Rmd file](code/Fig_PLA_FUS.Rmd).
  * [CellProfiler pipeline(PLA_Nuclear.cppipe)](code/).
  * [raw data of POLD1_U2OS_CSK_Nuclei.csv](data/POLD1_U2OS_CSK_Nuclei.csv).
  * [raw data of POLD1_Clone110_CSK_Nuclei.csv](data/POLD1_Clone110_CSK_Nuclei.csv).
  * [raw data of PCNA_U2OS_CSK_Nuclei.csv](data/PCNA_U2OS_CSK_Nuclei.csv).
  * [raw data of PCNA_Clone110_CSK_Nuclei.csv](data/PCNA_Clone110_CSK_Nuclei.csv).
  * [raw data of FEN1_U2OS_CSK_Nuclei.csv](data/FEN1_U2OS_CSK_Nuclei.csv).
  * [raw data of FEN1_Clone110_CSK_Nuclei.csv](data/FEN1_Clone110_CSK_Nuclei.csv).
* __Fig. 7:__ [FUS regulates DNA replication timing in U-2 OS cells](/barplot_RT_FUS-KO.md)
  * [Rmd file](code/barplot_RT_FUS-KO.Rmd).
  * [raw data of 7B](data/RT_EdU_FUS-KO.csv).
  * [raw data of 7D](data/RT_BrdU_DoubleThymidine_FUS.csv).
* __Fig. 9D:__ [Doughnut Pie chart of FADs coverage](/PieChart_FADs_Coverage.md)
  * [Rmd file](code/PieChart_FADs_Coverage.Rmd).
  * [raw data of 9D](data/FADs_coverage.csv).
* __Sup. Fig. 2B:__ [IRIF of BRCA1 in FUS knockout cells](/fig_IRIF_BRCA1_FUS-KO.md)
  * [Rmd file](code/fig_IRIF_BRCA1_FUS-KO.Rmd).
  * [raw data](data/IRIF_BRCA1_mock_15min.csv).
* __Sup. Fig. 2D:__ [IRIF of BRCA1 in CyclinA postive cells](/fig_IRIF_BRCA1_CyclinA.md)
  * [Rmd file](code/fig_IRIF_BRCA1_CyclinA.Rmd).
  * [raw data](data/IRIF_BRCA1_CyclinA.csv).
* __Sup. Fig. 5A and B:__ [RNA SEQ analysis of FUS](/Fig_FUS_RnaSeqDESeq2.md)
  * [Rmd file](code/Fig_FUS_RnaSeqDESeq2.Rmd).
  * [Count Matrix by FeatureCounts](data/fus_featurecounts.txt.Rmatrix.txt).
  * [Normalized RNA-seq counts by DESeq2](data/normalized_counts.csv).
  * [FUS specific regulated genes](data/FusSpeRegulatedGenes.csv).
  * [DEG list of WT and KO cells](data/sigWTvsKO_DESeq2.csv).
  * [DEG list of WT and RE cells](data/sigWTvsRE_DESeq2.csv).
  * [DEG list of RE and KO cells](data/sigREvsKO_DESeq2.csv).

* __Fig. 8, 9 and Sup. Fig. 8:__ Analysis of FUS replication timing
    * Code for QC, mapping and log2(S/G1): [RT_FUS.sh](code/RT_FUS.sh).
    * RT signal Normalization and smoothing: [normalization_smoothing.R](code/normalization_smoothing.R).
    * [RT plot of Fig. 8 and Sup. Fig. 8A to D](/Fig_RT_Plots20190102.md).
      * Rmd file: [Fig_RT_Plots20190102.Rmd](code/Fig_RT_Plots20190102.Rmd).
      * Raw data of smoothing data: [all_sample_r1_r2_smooth.txt](data/all_sample_r1_r2_smooth.txt).
    * [RT domain plots for Fig. 9A to C](/Fig_RT_domain_plots.md).
      * Rmd file: [Fig_RT_domain_plots.Rmd](code/Fig_RT_domain_plots.Rmd).
      * Raw data of segway results:
        * [U2OS, FUS<sup>+/+</sup>](data/U2OS_segway.bed)
        * [Clone110, FUS<sup>-/-</sup>](data/Clone110_segway.bed)
        * [FUS/Clonee110, FUS<sup>-/-</sup>:FUS](data/FUSClone110_segway.bed)
