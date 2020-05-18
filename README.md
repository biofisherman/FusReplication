# Fused in sarcoma regulates DNA replication timing and progression
## Code for figures:
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
