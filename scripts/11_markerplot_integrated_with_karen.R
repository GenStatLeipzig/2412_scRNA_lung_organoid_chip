rm(list = ls())
require(toolboxH)
require(here)
require(ggplot2)
require(patchwork)
require(Seurat)
require(patchwork)

initializeSkript()

options(width = 222)

source(here("../../../../07_programme/github/scRNATexMex/R/scRNA_functions_25-05-15.R"))


# # load data ----
# Assuming your lung data is stored in an RDS file - update path as needed
all10x4 = readRDS(here("results/10_1_anno_all10x_celltypes_qced_integrated_withUnstretched.rds"))
all10x4
anno_all10x4 = fread(here("results/10_1_anno_all10x_celltypes_qced_integrated_withUnstretched.txt.gz"))
anno_all10x4[is.na(celltype_provis), celltype_provis := "NA"]

# ## sankey and Dim plot ----
ggplotSankey(anno_all10x4[,.(celltype_chatti4_1_v3_rpca, sample )])

ggplotSankey(anno_all10x4[,.(celltype_provis, celltype_chatti4_1_v3_rpca_pur )]) +
  ggplotSankey(anno_all10x4[grepl("Unknown",celltype_chatti4_1_v3_rpca_pur),.(celltype_provis, celltype_chatti4_1_v3_rpca )])


# save plot
# ggsave(here("results/s1380_lung_sankey_pp_vs_rs_250409.pdf"), width = 10, height = 20)
DimPlot(all10x4, group.by = c('sample',"celltype_provis", "unintegrated_clusters"), label = T, reduction = "umap.unintegrated") + plot_annotation(title = "Unintegrated merged data")& NoLegend()  


p_dim =   DimPlot(all10x4, group.by = c('sample',"celltype_chatti4_1_v3_rpca_pur", "rpca_clusters"), label = T, reduction = "umap.rpca") + plot_annotation(title = "RPCA integrated merged data")& NoLegend()  
p_dim 
plotly::ggplotly(DimPlot(all10x4, group.by = c( "celltype_chatti4_1_v3_rpca_pur" ), label = T, reduction = "umap.rpca") + plot_annotation(title = "RPCA integrated merged data"))

p_feature =  FeaturePlot(all10x4, c('SFTPC',"AGER"), label = T, reduction = "umap.rpca") + plot_annotation(title = "RPCA integrated merged data") 
 p_feature
 
 plotly::ggplotly(p_feature)

 VlnPlot(all10x4, c('SFTPC',"AGER"), group.by = c('rpca_clusters'), alpha = 0.1 ) + plot_annotation(title = "RPCA integrated merged data")
 
 epithel_markers <- rbind(
   data.table(celltype = "Basal epithelial", markers = c("KRT5", "TP63", "EGFR", "KRT14")),  # KRT5/TP63 are canonical markers
   data.table(celltype = "Type 2 alveolar", markers = c("SFTPC", "SFTPB", "SFTA2", "NAPSA")),  # SFTPC is canonical specific marker
   data.table(celltype = "Type 1 alveolar", markers = c("AGER", "HOPX", "CADM1", "RBMS3")),  # AGER is canonical specific marker
   data.table(celltype = "Club", markers = c("SCGB1A1", "CYP2B7P", "TMC5", "NEBL")),  # SCGB1A1 is canonical but not in data
   data.table(celltype = "Endothelial", markers = c("CDH5", "ERG", "CALCRL", "TIE1")),  # CDH5 is canonical specific marker
   data.table(celltype = "Fibroblasts", markers = c("VIM", "LGALS1", "COL1A1", "FABP5")),  # VIM/COL1A1 are canonical
   data.table(celltype = "Epithelial", markers = c("EPCAM", "KRT8", "KRT19", "CDH1"))  # EPCAM is canonical specific marker
 )

 marker_plot <- doMarkerDotPlot(seurat = all10x4, 
                                         marker_groups_peter = epithel_markers, 
                                         grouping_factor = 'celltype_chatti4_1_v3_rpca_pur'
 ) 
 
 # Save plots
 ggsave(plot = marker_plot, filename =  here("results/11_rpcaintegrated_dotplot_standard_celltypes_chatti4_1.jpeg"), width = 16, height = 8)
 
 ggsave(plot = marker_plot, filename =  here("results/11_rpcaintegrated_dotplot_standard_celltypes_chatti4_1.pdf"), width = 16, height = 8)
 
 
 
 marker_plot_feiner <- doMarkerDotPlot(seurat = all10x4, 
                                         marker_groups_peter = epithel_markers, 
                                         grouping_factor = 'celltype_chatti4_1_v3_rpca'
 ) 
 
 
 # Save plots
 ggsave(plot = marker_plot_feiner, filename =  here("results/11_rpcaintegrated_dotplot_standard_clusters_chatti4_1.jpeg"), width = 16, height = 8)
 
 ggsave(plot = marker_plot_feiner, filename =  here("results/11_rpcaintegrated_dotplot_standard_clusters_chatti4_1.pdf"), width = 16, height = 8)
 
 
 
 
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
                                         grouping_factor = 'celltype_chatti4_1_v3_rpca_pur'
 ) 
 
 # Save plots
 ggsave(plot = marker_plot_liu, filename =  here("results/11_rpcaintegrated_dotplot_liu_celltype_chatti4_1.jpeg"), width = 16, height = 8)
 
 ggsave(plot = marker_plot_liu, filename =  here("results/11_rpcaintegrated_dotplot_liu_celltype_chatti4_1.pdf"), width = 16, height = 8)

 
 marker_plot_liu_feiner <- doMarkerDotPlot(seurat = all10x4, 
                                    marker_groups_peter = liu_markers, 
                                    grouping_factor = 'celltype_chatti4_1_v3_rpca'
 ) 
 
 
 # Save plots
 ggsave(plot = marker_plot_liu_feiner, filename =  here("results/11_rpcaintegrated_dotplot_liu_clusters_chatti4_1.jpeg"), width = 16, height = 8)
 
 ggsave(plot = marker_plot_liu_feiner, filename =  here("results/11_rpcaintegrated_dotplot_liu_clusters_chatti4_1.pdf"), width = 16, height = 8)
 
 
 
 finalizeSkript()
