require(toolboxH)
require(here)
require(Seurat)
require(ggplot2)
require(ggthemes)

source(here("../../../07_programme/github/scRNATexMex/R/scRNA_functions_25-07-04.R"))
# Get all analysis folders
folders <- list.dirs(here("../gitOrganoid_hk/demuxify"), recursive = TRUE)
folders <- folders[grepl("S[0-9]", basename(folders))]


# Get method subfolders
method_folders <- list.dirs(folders, recursive = FALSE, full.names = TRUE)

alldoublet = lapply(method_folders, function(method_path) {
  print(method_path)
  # method_path = method_folders[38]
  folder_name <- str_split(method_path, "demuxify/") %>% sapply("[", 2) %>% str_split("/") %>% sapply("[", 1)
  method_name <- basename(method_path)
  
  # Look for summary and results files
  file_summary <- list.files(method_path, pattern = ".*summary\\.tsv", full.names = TRUE)
  if(length(file_summary)==0) return(c())
  print(file_summary)
  # file_summary =   "J:/02_projekte/2412_scRNA_lung_organoid_chip/gitOrganoid_hkwin/../gitOrganoid_hk/demuxify/S1/scDblFinder/scDblFinder_doublet_summary.tsv"
  res_sum = fread(file_summary)
  # setnames(res_sum, c("Barcode", "doublet", "score"))
  res_sum$method = method_name
  res_sum$dataset = folder_name
  res_sum$fn = method_path
  
  
  res_sum
  
  file_results <- list.files(method_path, pattern = ".*results|.*singlets\\.tsv", full.names = TRUE)
  print(file_results)
  
  res_results = fread(file_results)
  if(method_name=="scds") setnames(res_results, c("Barcode","score", "doublet")) else setnames(res_results, c("Barcode", "doublet", "score"))
  res_results$method = method_name
  res_results$dataset = folder_name
  res_results$fn = file_results
  
  res_results
  
  res = c()  
  res$summary = res_sum
  res$results = res_results
  res
}
)

alldoublet_summary = lapply(alldoublet, function(x) x$summary) %>% rbindlist(fill = TRUE)
alldoublet_results = lapply(alldoublet, function(x) x$results) %>% rbindlist(fill = TRUE)

alldoublet_results[, Barcode2 := paste0("karenstretch_",dataset, "_", Barcode)]
stopifnot(nrow(alldoublet_results[allDuplicatedEntries(paste(Barcode2, method))])==0)



alldoublet_summary[, doublet_rate := `Droplet N` / sum(`Droplet N`), by = .(method, dataset)]
alldoublet_summary[, doublet_rate := fifelse(Classification == "doublet", doublet_rate, 0)]

ggplot(alldoublet_summary[Classification == "doublet"], aes(x = method, y = doublet_rate, fill = dataset)) +
  geom_col(position = "dodge") +
  labs(title = "Doublet Detection Rates by Method and Dataset",
       x = "Method", y = "Doublet Rate", fill = "Dataset") +
  theme_minimal() +
  scale_y_continuous(labels = label_percent(accuracy = 1))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


ggsave(here("results/18_doublet_detection_rates_per_sample.png"), width = 10, height = 6, dpi = 300)


require(UpSetR)
alldoublet_results[,.N, doublet]
upset_data <- alldoublet_results[doublet == "doublet", .(Barcode, method)]
upset_data[, method_dataset := paste0(method)]

binary_matrix <- dcast(upset_data, Barcode ~ method_dataset, value.var = "method_dataset", 
                       fun.aggregate = function(x) as.integer(length(x) > 0))

upset(binary_matrix, sets = colnames(binary_matrix)[-1], order.by = "freq")


pdf(here("results/18_doublet_upset_plot_per_sample.pdf"),width = 6, height = 5)
upset(binary_matrix, sets = colnames(binary_matrix)[-1], order.by = "freq")
dev.off()

# percluster ----
integrated = readRDS(here("results/13_1_anno_all10x_celltypes_qced_integrated_withUnstretched.rds"))
integrated

anno = integrated@meta.data %>% data.table(keep.rownames = "Barcode")
anno[,.N, .(seurat_clusters, celltype_opus_rpca)][order(seurat_clusters)]
integrated$seurat_clusters %>% as.character() %>% as.numeric() %>% unique() %>% sort()

# Von: Hoffmann, Karen <karen.hoffmann@charite.de> 
#   Gesendet: Mittwoch, 9. Juli 2025 11:20
# An: Holger Kirsten <holger.kirsten@imise.uni-leipzig.de>; Pennitz, Peter <peter.pennitz@charite.de>; Nouailles-Kursar, Geraldine <geraldine.nouailles@charite.de>
#   Betreff: [Extern] AW: AW: [ext] AW: Lung chip

karen_dt <- data.table(
  number = 0:31,
  celltype = c(
    "Basaloid cells",
    "Basaloid cells", 
    "AT1-like cells",
    "AT1-like cells",
    "AT1-like cells",
    "AT1-like cells",
    "Basaloid cells",
    "Basaloid cells",
    "AT1-like cells",
    "Basaloid cells",
    "AT1-like cells",
    "Epithelial cell",
    "Basal cells",
    "Basal cells",
    "AT1-like cells",
    "Endothelial cells",
    "Epithelial cells",
    "Secretory cells",
    "Epithelial cells",
    "Epithelial cells",
    "AT2 cells (from Lung chip)",
    "Epithelial cells",
    "AT1-like cells",
    "Epithelial cells",
    "AT2 cells (from Organoids)",
    "Endothelial cells",
    "Proliferating cells",
    "Proliferating cells",
    "Endothelial cells",
    "Epithelial cells",
    "Endothelial cells",
    "Neuroendocrine cells"
  )
)


integrated$celltype_250709 = karen_dt[match_hk(integrated$seurat_clusters, karen_dt$number),celltype]
anno$celltype_250709 = karen_dt[match_hk(anno$seurat_clusters, karen_dt$number),celltype]
ggplotSankey(anno[, .(celltype_250709,celltype_opus_rpca_pur,celltype_opus_rpca)])
DimPlot(integrated, group.by = c('sample',"seurat_clusters", "celltype_250709"), reduction = "umap.rpca") & NoLegend()  
customDimPlot(integrated, group_by = c("seurat_clusters"), reduction = "umap.rpca")

anno[seurat_clusters %in% c(20, 24), .N, .(celltype_opus_rpca,category)][order(celltype_opus_rpca,category)]

# add doublett info        
anno[, barcodeprefix := substr(Barcode, 1, nchar(Barcode) - 18)]
anno[, .N,.(category, barcodeprefix)]

qlist1 = venn2(alldoublet_results$Barcode2, anno[run10x %in% c("S1", "S2", "S3", "S4"),Barcode])
str(qlist1)

mytable(substr(qlist1$q1, 1, nchar(qlist1$q1) - 18))
mytable(substr(qlist1$q2, 1, nchar(qlist1$q2) - 18))
mytable(substr(qlist1$q3, 1, nchar(qlist1$q3) - 18))

# here q2 are QC filtered samplse and q3 are those, that were recovered by cellbender but not yet included in filtered 10x folder, which is used by pipeline. 

alldoublet_results[, qcfiltered:= Barcode2 %nin% anno$Barcode]
alldoublet_results[, uniqueN(Barcode2),qcfiltered]

alldoublet_results[, seurat_clusters := anno[match_hk(alldoublet_results$Barcode2, anno$Barcode), seurat_clusters]]

consensus_doublets <- alldoublet_results[doublet == "doublet", .N, by = Barcode2][N >= 3]
alldoublet_results[, consensus_doublet_2vote := Barcode2 %in% consensus_doublets$Barcode2]
alldoublet_results_unique  = alldoublet_results[,.(Barcode2, seurat_clusters,consensus_doublet_2vote)] %>% unique()

plot_data <- alldoublet_results_unique[!is.na(seurat_clusters), .(
  total = .N,
  doublet_count = sum(consensus_doublet_2vote)
), by = seurat_clusters]

plot_data[, doublet_prop := doublet_count / total]

ggplot(plot_data, aes(x = seurat_clusters, y = doublet_prop)) +
  geom_col(fill = "steelblue") +
  labs(title = "Consensus Doublet Rate by Seurat Cluster", 
       x = "Seurat Cluster", y = "Doublet Proportion") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# vice versa
alldoublet_results_wide = dcast.data.table(alldoublet_results, Barcode2 ~ method, value.var = c('doublet', "method"))
alldoublet_results_wide

anno2 = merge(anno, alldoublet_results_wide, by.x = "Barcode", by.y = "Barcode2", all.x = T, sort = FALSE)

consensus_doublets_2vote <- alldoublet_results[is.na(seurat_clusters)==F][doublet == "doublet", .N, by = Barcode2][N >= 2]
anno2[, consensus_doublet_2vote := Barcode %in% consensus_doublets_2vote$Barcode2]
anno2[, doublet_tested_again := Barcode %in% alldoublet_results$Barcode2]

plot_data_2vote <- anno2[ , .(
  total = .N,
  doublet_count = sum(consensus_doublet_2vote), 
  doublet_tested_again = sum(doublet_tested_again),
  unstretched = sum(category =='Unstretched ODAEC'), 
  seurat_clusters= unique(as.numeric(seurat_clusters))
  
), .(celltype = paste(seurat_clusters, celltype_250709))]

plot_data_2vote[, doublet_prop := doublet_count / total]
plot_data_2vote[, nodoublet_prop := (doublet_tested_again / total)-doublet_prop]
plot_data_2vote[, unstretched_prop := unstretched / total]


plot_data_2votem = melt(plot_data_2vote, measure.vars = c('doublet_prop', 'nodoublet_prop', "unstretched_prop"))
plot_data_2votem[, variable := factor(variable, levels = c('doublet_prop', 'nodoublet_prop', "unstretched_prop") %>% rev())]
p_2vote = ggplot(plot_data_2votem, aes(x = reorder(celltype,seurat_clusters), y = value, fill = variable)) +
  geom_col(position = "stack") +
  labs(title = "Consensus Doublet Rate by Seurat Cluster", 
       subtitle = "2-votes required among solo, scDblFinderm, scds, scrublet",
       x = "Seurat Cluster", y = "Doublet Proportion") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  scale_fill_manual(values = c("darkred", "lightblue", "grey") %>% rev()) + 
  geom_hline(yintercept = 0.5, col  = "darkred", lty = 2, alpha = 0.5)
p_2vote

consensus_doublets_3vote <- alldoublet_results[is.na(seurat_clusters)==F][doublet == "doublet", .N, by = Barcode2][N >= 3]
anno2[, consensus_doublet_3vote := Barcode %in% consensus_doublets_3vote$Barcode2]
anno2[, doublet_tested_again := Barcode %in% alldoublet_results$Barcode2]

plot_data_3vote <- anno2[ , .(
  total = .N,
  doublet_count = sum(consensus_doublet_3vote), 
  doublet_tested_again = sum(doublet_tested_again),
  unstretched = sum(category =='Unstretched ODAEC'), 
  seurat_clusters= unique(as.numeric(seurat_clusters))
  
), .(celltype = paste(seurat_clusters, celltype_250709))]

plot_data_3vote[, doublet_prop := doublet_count / total]
plot_data_3vote[, nodoublet_prop := (doublet_tested_again / total)-doublet_prop]
plot_data_3vote[, unstretched_prop := unstretched / total]


plot_data_3votem = melt(plot_data_3vote, measure.vars = c('doublet_prop', 'nodoublet_prop', "unstretched_prop"))
plot_data_3votem[, variable := factor(variable, levels = c('doublet_prop', 'nodoublet_prop', "unstretched_prop") %>% rev())]
p_3vote = ggplot(plot_data_3votem, aes(x = reorder(celltype,seurat_clusters), y = value, fill = variable)) +
  geom_col(position = "stack") +
  labs(title = "Consensus Doublet Rate by Seurat Cluster", 
       subtitle = "3-votes required among solo, scDblFinderm, scds, scrublet",
       x = "Seurat Cluster", y = "Doublet Proportion") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  scale_fill_manual(values = c("darkred", "lightblue", "grey") %>% rev()) + 
  geom_hline(yintercept = 0.5, col  = "darkred", lty = 2, alpha = 0.5)
p_3vote

combined_plot <- p_2vote + p_3vote
ggsave(here("results/18_doublet_consensus_comparison.png"), combined_plot, width = 20, height = 6, dpi = 300)

p_conse = ggplot(anno, aes(x = reorder(paste(seurat_clusters, celltype_250709),as.numeric(seurat_clusters)), fill = category)) +
  geom_bar(position = "fill")+
  theme_minimal() +
  facet_grid(.~ (seurat_clusters %in% c(20, 24)), scales = "free", space = "free")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  # scale_fill_manual(values = c("darkred", "lightblue", "grey") %>% rev()) + 
  scale_fill_viridis_d()+
  geom_hline(yintercept = 0.5, col  = "darkred", lty = 2, alpha = 0.5)
p_conse


ggsave(here("results/18_category_composition_by_cluster.jpeg"), p_conse, width = 12, height = 6, dpi = 300)


## exclude update ------
anno2[, QC_problem_doublett := seurat_clusters %in% 27]
anno2[, mytable(QC_problem_doublett)]
anno2[, celltype_250709b  := str_split(celltype_250709, " cells") %>% sapply("[", 1)] 


anno2[grep("cell", celltype_250709b),.N, celltype_250709b]
anno2[grep("cell", celltype_250709b), celltype_250709b := 'Epithelial ']
anno2[, .N, celltype_250709b]


integrated$celltype_250709b = anno2[match_hk(integrated$seurat_clusters, anno2$seurat_clusters, makeunique = T, importcol = anno2$celltype_250709b),celltype_250709b]
integrated$QC_problem_doublett = anno2[match_hk(integrated$seurat_clusters, anno2$seurat_clusters, makeunique = T, importcol = anno2$QC_problem_doublett),QC_problem_doublett]

integrated2 = integrated[, integrated$QC_problem_doublett==FALSE]
integrated
integrated2
NogpaletteReihe2 = 
c("#CB769E", "#DE639A", "#A85C85", "#0081AF",
  "#4F6D7A",
  "#7C6A0A", 
  "#7C6A0A",    "#F97E44", "#246A73", "#5CC1BC", "#62C370","#644536",
 "#F7C548",
  "#FB3640", "#B7245C",  "#3E2F5B", "#B2675E","#0D3B66"
)
dt = data.table(sort(integrated2$celltype_250709b%>% unique()) , NogpaletteReihe2[1:length(sort(integrated2$celltype_250709b%>% unique()))])
# Function to show hex codes in their actual colors
show_colored_hex <- function(dt) {
  require(crayon)
  options(crayon.enabled = TRUE)
  options(crayon.colors = 256)
  
  for(i in 1:nrow(dt)) {
    cell_type <- dt$V1[i]
    hex_color <- dt$V2[i]
    
    # Create colored text
    colored_text <- crayon::make_style(hex_color)(hex_color)
    cat(cell_type, colored_text, "\n")
  }
}


# Use the function
show_colored_hex(dt)


p_dimfin = customDimPlot(integrated2, group_by = c("celltype_250709b"), reduction = "umap.rpca", label.size = 4.5, colors = NogpaletteReihe2, label_darkening_alpha = 0.4)
p_dimfin
ggsave(here("results/18_dimplot_final_celltype_250709b.png"), p_dimfin, width = 12, height = 8, dpi = 300)
ggsave(here("results/18_dimplot_final_celltype_250709b.jpeg"), p_dimfin, width = 12, height = 8, dpi = 300)
ggsave(here("results/18_dimplot_final_celltype_250709b.pdf"), p_dimfin, width = 12, height = 8, dpi = 300)

fwrite(anno2, here("results/18_anno2_final.csv.gz"))
anno3 = integrated2@meta.data %>% data.table(keep.rownames = T)
saveRDS(integrated2, here("results/18_integrated_final_woClust27.rds"))
fwrite(anno3, here("results/18_integrated_final_woClust27.csv.gz"))



finalizeSkript()

