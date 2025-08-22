library(tidyverse)
library(ggplot2)
library(dplyr)
library(Seurat)
library(patchwork)
library(ggrepel)
require(toolboxH)
library(xlsx)

OrganoidChip <- readRDS("E:/Seurat_Objects/KarenOrganoide/18_integrated_final_woClust27.rds")

OrganoidChip@active.ident = factor(OrganoidChip$celltype_250709b)
Cellorder = c('AT1-like', 'AT2', 'Basal', 'Endothelial', 'Basaloid', 'Epithelial', "Epithelial ", 'Neuroendocrine', 'Secretory', "Proliferating")
Colourpalette <-  c("#4D7698", "#AED2D3", "#FFCB47", "#94B67F", "#A0783B", "#E99DA4", "#E99DA4", "#D2D32E", "#9C777A", "#CA5071")
my_levels <- c(Cellorder)
OrganoidChip@active.ident <- factor(x = OrganoidChip@active.ident, levels = my_levels)
DimPlot(OrganoidChip, label = T, reduction = "umap.rpca", cols = c(Colourpalette)) + NoLegend()
#ggsave("Results/UMAP_OrganoidChipColourpaletteNEW.pdf", width = 7, height = 7)

Colourpalette10xrun <-  c("#2098CB", "#6A78B1", "#714BAC","#9FDEB0", "#216F46", "#E2CC7D", "#F8E503")
DimPlot(OrganoidChip, reduction = "umap.rpca", group.by = "run10x", shuffle = T, cols = c(Colourpalette10xrun))

#Run DGE
OrganoidChip$celltype_group <- paste(OrganoidChip$celltype_250709b, OrganoidChip$run10x, sep = "_")
Idents(OrganoidChip) <- "celltype_group"
diff.markers <- FindMarkers(OrganoidChip, ident.1 = "Endothelial_S4", ident.2 = "Endothelial_S3")   # repeat for AT1-like, AT2 and Basaloid cells
dplyr::filter(diff.markers, diff.markers$p_val<0.05)
diff.markers <- tibble::rownames_to_column(diff.markers, "Gene")
diff.markers$Bonferoni_adj <- "FALSE"
diff.markers$Bonferoni_adj[diff.markers$p_val_adj < 0.05] <- "TRUE"
diff.markers$diffexpressed <- "FALSE"
diff.markers$diffexpressed[diff.markers$avg_log2FC > 0 & diff.markers$p_val < 0.05] <- "UP"
diff.markers$diffexpressed[diff.markers$avg_log2FC < 0 & diff.markers$p_val < 0.05] <- "DOWN"
#write.xlsx(diff.markers, "diff.markers.Endothelial_S4S3.xlsx")

# reading in DGE data 
At1like_S4S3 <- read_excel2("Data/diff.markers.AT1-like_S4S3.xlsx")
Basaloid_S2S1 <- read_excel2("Data/diff.markers.Basaloid_S2S1.xlsx")
AT2_S4S3 <- read_excel2("Data/diff.markers.At2_S4S3.xlsx")
Endothel_S4S3 <- read_excel2("Data/diff.markers.Endothelial_S4S3.xlsx")


#Overlapping genes
At1like_S4S3_SigUP <- dplyr::filter(At1like_S4S3, c((At1like_S4S3$Bonferoni_adj=="TRUE") & (At1like_S4S3$diffexpressed=="UP"))) 
At1like_S4S3_SigDown <- dplyr::filter(At1like_S4S3, c((At1like_S4S3$Bonferoni_adj=="TRUE") & (At1like_S4S3$diffexpressed=="DOWN"))) 

Basaloid_S2S1_SigUP <- dplyr::filter(Basaloid_S2S1, c((Basaloid_S2S1$Bonferoni_adj=="TRUE") & (Basaloid_S2S1$diffexpressed=="UP"))) 
Basaloid_S2S1_SigDown <- dplyr::filter(Basaloid_S2S1, c((Basaloid_S2S1$Bonferoni_adj=="TRUE") & (Basaloid_S2S1$diffexpressed=="DOWN")))

mycolors5 =c("red", "salmon")
options(max.print=1000000)
venn2(At1like_S4S3_SigUP$Gene, Basaloid_S2S1_SigUP$Gene, fill_colors = alpha(mycolors5, 0.2), line_colors = mycolors5, text_colors = mycolors5)

# pdf("Results/Venn Diagram SigUP overlap_V2.pdf", 10,10)
# venn2(At1like_S4S3_SigUP$Gene, Basaloid_S2S1_SigUP$Gene, fill_colors = alpha(mycolors5, 0.2), line_colors = mycolors5, text_colors = mycolors5)
# dev.off()

mycolors6 =c("darkblue", "lightblue")
venn2(At1like_S4S3_SigDown$Gene, Basaloid_S2S1_SigDown$Gene, fill_colors = alpha(mycolors6, 0.2), line_colors = mycolors6, text_colors = mycolors6)

# pdf("Results/Venn Diagram SigDown overlap_V2.pdf", 10,10)
# venn2(At1like_S4S3_SigDown$Gene, Basaloid_S2S1_SigDown$Gene, fill_colors = alpha(mycolors6, 0.2), line_colors = mycolors6, text_colors = mycolors6)
# dev.off()



#Volcano Plots
#At1-like
custom_genes_At1like_S4S3 <- At1like_S4S3 %>% 
  filter(Gene %in% c("COL1A1", "COL1A2", "EMP2", "CLIC5", "SFTPA1", "SFTPB", "VIM", "ZEB1", "FTH1", "SCGB3A2"))

p <- ggplot(data=At1like_S4S3, aes(x=avg_log2FC, y=-log10(p_val),col=diffexpressed, label=Gene)) + geom_point(size = 4) + geom_point(data=custom_genes_At1like_S4S3 , aes(x=avg_log2FC,y=-log10(p_val)), col="pink", size = 4) + theme_classic() + geom_text_repel(size = 8, max.overlaps = 1) + geom_vline(xintercept=c(-1, 1), linetype = "dashed") + geom_hline(yintercept=-log10(0.05), linetype = "dashed") + scale_colour_manual(values=c("blue", "black", "red"))  + geom_text_repel(data = custom_genes_At1like_S4S3  %>%   filter(Gene %in% c("COL1A1", "COL1A2", "EMP2", "CLIC5", "SFTPA1", "SFTPB", "VIM", "ZEB1", "FTH1", "SCGB3A2")), size = 8) + NoLegend()

p
#ggsave("Results/Volcano_IntNoStretch_At1like_S4S3_Update2.pdf", width = 15, height = 15)


#AT2
custom_genes_At2_S4S3 <- AT2_S4S3 %>% 
  filter(Gene %in% c("IL6", "CXCL8"))

p1 <- ggplot(data=AT2_S4S3, aes(x=avg_log2FC, y=-log10(p_val),col=diffexpressed, label=Gene)) + geom_point(size = 4) + geom_point(data=custom_genes_At2_S4S3 , aes(x=avg_log2FC,y=-log10(p_val)), col="pink", size = 4) + theme_classic() + geom_text_repel(size = 8, max.overlaps = 1) + geom_vline(xintercept=c(-1, 1), linetype = "dashed") + geom_hline(yintercept=-log10(0.05), linetype = "dashed") + scale_colour_manual(values=c("blue", "black", "red"))  + geom_text_repel(data = custom_genes_At2_S4S3  %>%   filter(Gene %in% c("IL6", "CXCL8")), size = 8) + NoLegend()

p1
#ggsave("Results/Volcano_IntNoStretch_AT2_S4S3_Update2.pdf", width = 15, height = 15)


#Basaloid
custom_genes_Basaloid_S2S1 <- Basaloid_S2S1 %>% 
  filter(Gene %in% c("CDH1", "CDH5", "PECAM1", "VWF", "KRT81", "KRT19", "KRT18", "KRT7", "S100A2", "S100A14", "S100A6", "S100A11", "S100A10"))

p2 <- ggplot(data=Basaloid_S2S1, aes(x=avg_log2FC, y=-log10(p_val),col=diffexpressed, label=Gene)) + geom_point(size = 4) + geom_point(data=custom_genes_Basaloid_S2S1 , aes(x=avg_log2FC,y=-log10(p_val)), col="pink", size = 4) + theme_classic() + geom_text_repel(size = 8, max.overlaps = 1) + geom_vline(xintercept=c(-1, 1), linetype = "dashed") + geom_hline(yintercept=-log10(0.05), linetype = "dashed") + scale_colour_manual(values=c("blue", "black", "red"))  + geom_text_repel(data = custom_genes_Basaloid_S2S1  %>%   filter(Gene %in% c("CDH1", "CDH5", "PECAM1", "VWF", "KRT81", "KRT19", "KRT18", "KRT7", "S100A2", "S100A14", "S100A6", "S100A11", "S100A10")), size = 8) + NoLegend()

p2
#ggsave("Results/Volcano_IntNoStretch_Basaloid_S2S1_Update2.pdf", width = 15, height = 15)


#Endothelial
custom_genes_Endothel_S4S3 <- Endothel_S4S3 %>% 
  filter(Gene %in% c("CAV1", "IGFBP2", "WNT4", "CXCL10"))

p3 <- ggplot(data=Endothel_S4S3, aes(x=avg_log2FC, y=-log10(p_val),col=diffexpressed, label=Gene)) + geom_point(size = 4) + geom_point(data=custom_genes_Endothel_S4S3 , aes(x=avg_log2FC,y=-log10(p_val)), col="pink", size = 4) + theme_classic() + geom_text_repel(size = 8, max.overlaps = 1) + geom_vline(xintercept=c(-1, 1), linetype = "dashed") + geom_hline(yintercept=-log10(0.05), linetype = "dashed") + scale_colour_manual(values=c("blue", "black", "red"))  + geom_text_repel(data = custom_genes_Endothel_S4S3  %>%   filter(Gene %in% c("CAV1", "IGFBP2", "WNT4", "CXCL10")), size = 8) + NoLegend()

p3
#ggsave("Results/Volcano_IntNoStretch_Endothel_S4S3_Update2.pdf", width = 15, height = 15)


#Barplot

anno_all10x4 = fread(("Data/18_integrated_final_woClust27.csv.gz")) %>% .[seurat_clusters != 27]
sort(unique(anno_all10x4$celltype_250709b))
# all10x4$celltype_250709b %>% unique() %>% sort()

anno_all10x4[, stretched := fcase(
  sample == "S1", "5%",
  sample == "S2", "10%",
  sample == "S3", "5%",
  sample == "S4", "10%",default = '0%'
)]

anno_all10x4[, origin := ifelse(grepl("CPAEC", category), "CPAEC", 
                                ifelse(grepl("ODAEC", category), "ODAEC", category)
)]


anno_all10x4[, stretched := factor(stretched, levels = c("0%", "5%", "10%"))]
anno_all10x4[,.N, .(category, origin, stretched, sample)]
anno_all10x4[,.N, .(origin, stretched)]

anno_all10x4[, experiment := paste(origin, stretched)]
unique(anno_all10x4$experiment)
anno_all10x4[, experiment := factor(experiment, levels =c("ODAEC 0%","ODAEC 5%", "ODAEC 10%",  "CPAEC 5%", "CPAEC 10%") %>% rev())]

colors <- c(
  "ODAEC 0%" = "#2098CB",    # one color for all baseline Organoids
  "ODAEC 5%" = "#E2CC7D",    # S3= ODAEC 5% 
  "ODAEC 10%" = "#F8E503",   # S4=ODAEC 10% 
  "CPAEC 5%" = "#9FDEB0",    # S1=CPAEC 5%
  "CPAEC 10%" = "#216F46"    # S2=CPAEC 10% 
)


p_prop = ggplot(anno_all10x4, aes(x = celltype_250709b, fill = experiment)) +
  geom_bar(position = "fill", alpha= 0.9) +
  
  scale_y_continuous(labels = scales::percent, breaks = pretty_breaks(8)) +
  labs(x = "", y = "Proportion", fill = "") +
  theme_minimal(base_size = 20) +
  # scale_fill_tableau(drop = FALSE)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  scale_fill_manual(values = colors) 


p_prop

# ggsave(here("results/21_barplot_cells.jpeg"), p_prop, width = 7, height = 8, dpi = 300)
# ggsave(here("results/21_barplot_cells.pdf"), p_prop, width = 7, height = 8)



#FeaturePlots
FeaturePlot(OrganoidChip, features = "SFTPC", reduction = "umap.rpca")
#repeat for SFTPC; LPCAT1; NAPSA; NKX2-1; HOPX; EMP2; AGER (RAGE Receptor); CLIC5; KRT17; KRT8; TP63; KRT5; SCGB1A1; EREG; PECAM1; TOP2A
ggsave("Results/FeaturePlot_SFTPC.pdf", width = 7, height = 7)