setwd("/Users/weiyanjia/SynologyDrive/Randal S. Tibbetts/FUS project/Bioinfo_analysis/Replication_domain_ideatification/Resolution_10000/Anotation results/analysis results/GO analysis/RefSeq_3k/clusterProfiler")
getwd()
# GO ANALYSIS
# library
library(clusterProfiler)
library(org.Hs.eg.db)

# Load in data
ERD_MRD_LRD <- read.delim("ERD_MRD_LRD_ENTREID.txt", sep="\t", header = TRUE)
lapply(ERD_MRD_LRD, head)

#EnrichGO
#BP
gc_enrich_GO_BP <- compareCluster(geneCluster = ERD_MRD_LRD, 
                                  fun           = "enrichGO",
                                  OrgDb         = org.Hs.eg.db,
                                  keyType       = "ENTREZID",
                                  ont           = "BP",
                                  pAdjustMethod = "BH",
                                  pvalueCutoff  = 0.05,
                                  qvalueCutoff  = 0.05,
                                  readable      = TRUE)
#MF
gc_enrich_GO_MF <- compareCluster(geneCluster = ERD_MRD_LRD, 
                                  fun           = "enrichGO",
                                  OrgDb         = org.Hs.eg.db,
                                  keyType       = "ENTREZID",
                                  ont           = "MF",
                                  pAdjustMethod = "BH",
                                  pvalueCutoff  = 0.05,
                                  qvalueCutoff  = 0.05,
                                  readable      = TRUE)

#Plot

## BF
dotplot(gc_enrich_GO_BP,  showCategory = 5, color = "pvalue", title = "GO Analysis of FUS related RTD genes(BP) ")

## MF
dotplot(gc_enrich_GO_MF,  showCategory = 5, color = "pvalue", title = "GO Analysis of FUS related RTD genes(MF) ")

#write data
# write.table(gc_enrich_GO_BP, "gc_enrich_GO_BP_genes.txt", sep = "\t", col.names = T, row.names = F, quote =F)
# write.table(gc_enrich_GO_MF, "gc_enrich_GO_MF_genes.txt", sep = "\t", col.names = T, row.names = F, quote =F)