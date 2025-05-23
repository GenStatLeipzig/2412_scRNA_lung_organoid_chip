rm(list = ls())
require(toolboxH)
require(here)
require(ggplot2)
require(patchwork)
require(Seurat)
require(patchwork)

initializeSkript()

options(width = 222)

source(here("../../../../07_programme/github/scRNATexMex/R/scRNA_functions_25-05-09.R"))


# # load data ----
# Assuming your lung data is stored in an RDS file - update path as needed
all10x4 = readRDS(here("results/06_1_all10x4_celltypes_qced_integrated.rds"))
all10x4
anno_all10x4 = fread(here("results/06_1_anno_all10x4_celltypes_qced_integrated.txt.gz"))

# ## sankey and Dim plot ----
ggplotSankey(anno_all10x4[,.(sample, celltype_chatti4_1_v2_rpca)])
ggplotSankey(anno_all10x4[,.(celltype_chatti4_1_v2_rpca, sample )])

  ggplotSankey(anno_all10x4[,.(celltype_chatti4_1_purcons, celltype_chatti4_1_v2_rpca_pur )]) +
  ggplotSankey(anno_all10x4[grepl("Unknown",celltype_chatti4_1_v2_rpca_pur),.(celltype_chatti4_1_purcons, celltype_chatti4_1_v2_rpca )])


# save plot
# ggsave(here("results/s1380_lung_sankey_pp_vs_rs_250409.pdf"), width = 10, height = 20)
require(patchwork)
  DimPlot(all10x4, group.by = c('sample',"celltype_chatti4_1_purcons", "unintegrated_clusters"), label = T, reduction = "umap.unintegrated") + plot_annotation(title = "Unintegrated merged data")& NoLegend()  


p_dim =   DimPlot(all10x4, group.by = c('sample',"celltype_chatti4_1_v2_rpca_pur", "rpca_clusters"), label = T, reduction = "umap.rpca") + plot_annotation(title = "RPCA integrated merged data")& NoLegend()  
p_dim 
plotly::ggplotly(DimPlot(all10x4, group.by = c( "celltype_chatti4_1_v2_rpca_pur" ), label = T, reduction = "umap.rpca") + plot_annotation(title = "RPCA integrated merged data"))

p_feature =  FeaturePlot(all10x4, c('SFTPC',"AGER"), label = T, reduction = "umap.rpca") + plot_annotation(title = "RPCA integrated merged data") 
 p_feature
 
 plotly::ggplotly(p_feature)

 VlnPlot(all10x4, c('SFTPC',"AGER"), group.by = c('rpca_clusters'), alpha = 0.1 ) + plot_annotation(title = "RPCA integrated merged data")
 
 ggplotSankey(anno_all10x4[,.(RNA_snn_res.0.5, unintegrated_clusters )])
 ggplotSankey(anno_all10x4[,.(RNA_snn_res.0.5, celltype_chatti4_1_purcons )])

 ggplotSankey(anno_all10x4[,.(unintegrated_clusters, celltype_chatti4_1_purcons )]) 
 
 p_dim =   DimPlot(all10x4, group.by = c('sample',"celltype_chatti4_1_purcons","RNA_snn_res.0.5", "unintegrated_clusters"), label = T, reduction = "umap.rpca") + plot_annotation(title = "RPCA integrated merged data")& NoLegend()  
 p_dim 
 
 # Create marker table for mouse lung celltypes
 # epithel_markers <- rbind(
 #   data.table(celltype = "Liu 3", markers = c("Scgb3a2", "Cyp2f2", "Scgb1a1", "Cbr2")),  # Column 3 - yellow markers
 #   data.table(celltype = "Liu 6", markers = c("Gm26917", "mt-Atp8", "Abca3", "Hc")),  # Column 6 - yellow markers
 #   data.table(celltype = "Liu 2", markers = c(  "Lyz2", "Sftpc", "Cxcl15")),  # Column 2 - yellow markers
 #   data.table(celltype = "Liu 4", markers = c("Hc", "Lyz2", "Sftpc", "Cxcl15", "Chit1", "Il33", "Cd74", "H2-Aa")),  # Column 4 - yellow markers
 #   data.table(celltype = "Liu 5", markers = c( "Ly6a", "H2afz", "Ereg")),  # Column 5 - yellow markers
 #   data.table(celltype = "Liu 1", markers = c("Cldn4", "Clu", "Lgals3", "Krt8")),  # Column 1 - yellow markers
 #   data.table(celltype = "Liu 7", markers = c("Hopx", "Akap5", "Ager", "Rtkn2"))  # Column 7 - yellow markers
 # )
 # 
 
 
 # 
 # 
 # 
 # 
 # 
 # 
 # 
 # 
 # 
 # 
 liu_markers <- rbind(
   data.table(celltype = "3: BASC-like cell", markers = c("SCGB3A2", "CYP2F1", "SCGB1A1")),  # Column 3 - yellow markersn no "Cbr2" https://www.alliancegenome.org/gene/MGI:107200
   data.table(celltype = "6: Transient cell", markers = c( "MT-ATP8", "ABCA3", "C5")),  # Column 6 - yellow markers (Note: Gm26917 none, Hc -> C5)
   data.table(celltype = "2: Cxcl15 - high AT2 cell", markers = c("LYZ", "SFTPC" )),  # Column 2 - yellow markers no ortho for "Cxcl15"
   data.table(celltype = "4: Chit1 - high AT2 cell", markers = c("C5", "LYZ", "SFTPC",   "CHIT1", "IL33", "CD74", "HLA-DQA1", "HLA-DQA2")),  # Column 4 - yellow markers (H2-Aa -> HLA-DQA2 HLA-DQA1)
   data.table(celltype = "5: Ereg - high AT2 cell", markers = c("LY6H", "H2AZ1", "EREG")),  # Column 5 - yellow markers (Ly6a -> LY6H)
   data.table(celltype = "1: Cldn4 - high PreAT1 cell", markers = c("CLDN4", "CLU", "LGALS3", "KRT8")),  # Column 1 - yellow markers
   data.table(celltype = "7: AT1 cell", markers = c("HOPX", "AKAP5", "AGER", "RTKN2"))  # Column 7 - yellow markers
 )

 marker_plot_liu <- doMarkerDotPlot(seurat = all10x4, 
                                         marker_groups_peter = liu_markers, 
                                         grouping_factor = 'celltype_chatti4_1'
 ) 
 
 # Save plots
 ggsave(plot = marker_plot_liu, filename =  here("results/09b_unintegrated_dotplot_liu_celltype_chatti4_1.jpeg"), width = 16, height = 8)
 
 ggsave(plot = marker_plot_liu, filename =  here("results/09b_unintegrated_dotplot_liu_celltype_chatti4_1.pdf"), width = 16, height = 8)

 
 marker_plot_liu_feiner <- doMarkerDotPlot(seurat = all10x4, 
                                    marker_groups_peter = liu_markers, 
                                    grouping_factor = 'unintegrated_clusters'
 ) 
 
 
 # Save plots
 ggsave(plot = marker_plot_liu_feiner, filename =  here("results/09b_unintegrated_dotplot_liu_unintegrated_clusters.jpeg"), width = 16, height = 8)
 
 ggsave(plot = marker_plot_liu_feiner, filename =  here("results/09b_unintegrated_dotplot_liu_unintegrated_clusters.pdf"), width = 16, height = 8)
 
 
 
 finalizeSkript()
