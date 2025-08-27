rm(list = ls())
require(toolboxH)
require(here)
require(ggplot2)
require(patchwork)
require(Seurat)
require(patchwork)

initializeSkript()

options(width = 222)

source(here("../../../../07_programme/github/scRNATexMex/R/scRNA_functions_25-07-04git.R"))


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
 
 # Create marker table for mouse lung celltypes
 epithel_markers <- rbind(
   data.table(celltype = "Basal epithelial", markers = c("KRT5", "TP63", "EGFR", "KRT14")),  # KRT5/TP63 are canonical markers
   data.table(celltype = "Type 2 alveolar", markers = c("SFTPC", "SFTPB", "SFTA2", "NAPSA")),  # SFTPC is canonical specific marker
   data.table(celltype = "Type 1 alveolar", markers = c("AGER", "HOPX", "CADM1", "RBMS3")),  # AGER is canonical specific marker
   data.table(celltype = "Club", markers = c("SCGB1A1", "CYP2B7P", "TMC5", "NEBL")),  # SCGB1A1 is canonical but not in data
   data.table(celltype = "Endothelial", markers = c("CDH5", "ERG", "CALCRL", "TIE1")),  # CDH5 is canonical specific marker
   data.table(celltype = "Fibroblasts", markers = c("VIM", "LGALS1", "COL1A1", "FABP5")),  # VIM/COL1A1 are canonical
   data.table(celltype = "Epithelial", markers = c("EPCAM", "KRT8", "KRT19", "CDH1"))  # EPCAM is canonical specific marker
 )

# Create species-specific plots - Basic version
 marker_plot_basic_v2 <- doMarkerDotPlot(seurat = all10x4, 
                                         marker_groups_peter = epithel_markers, 
                                         grouping_factor = 'celltype_chatti4_1_v2_rpca',
                                         custom_y_axis = 
                                           c("11: Endothelial cells", "12: Endothelial cells", "7: Endothelial cells", "10: Epithelial cells (Luminal, possible secretory/absorptive, i.e., glandular-like, MUC1+, CLDN7+, EPCAM+)", "14: Epithelial cells (Simple/luminal, KRT7+/KRT8+, CEACAM6+)", 
                                             
                                             "8: Epithelial cells (Luminal, likely glandular/absorptive, CLDN7+, EPCAM+, LSR+)",
                                             "3: Alveolar Type II cells", "1: Alveolar Type II cells", "2: Alveolar Type II cells", 
                                             "5: Alveolar Type II cells", "13: Unknown", "4: Basal Epithelial cells", 
                                             
                                             "9: Unknown", "0: Basal Epithelial cells", "6: Epithelial cells (Non-basal, cell-junction and cytoskeletal, likely intermediate/differentiating)"
                                           )
 ) 
 

# Save plots
 ggsave(plot = marker_plot_basic_v2, filename =  here("results/08_rpca_dotplot.jpeg"), width = 24, height = 12)
 
 ggsave(plot = marker_plot_basic_v2, filename =  here("results/08_rpca_dotplot.pdf"), width = 24, height = 12)

 
 finalizeSkript()
