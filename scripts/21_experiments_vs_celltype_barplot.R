rm(list = ls())
require(toolboxH)
require(here)
require(ggplot2)
require(patchwork)
require(Seurat)
require(patchwork)

initializeSkript()

options(width = 222)

source(here("../../../../07_programme/github/scRNATexMex/R/scRNA_functions_25-07-04.R"))


# # load data ----
# Assuming your lung data is stored in an RDS file - update path as needed
# all10x4 = readRDS(here("results/18_integrated_final_woClust27.rds"))
# all10x4
# all10x4$celltype_250709b %>% unique() %>% sort()


anno_all10x4 = fread(here("results/18_integrated_final_woClust27.csv.gz")) %>% .[seurat_clusters != 27]
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
  "ODAEC 0%" = "#2E5C8A",    # Dark blue (baseline)
  "ODAEC 5%" = "#5B9BD5",    # Medium blue 
  "ODAEC 10%" = "#8BB6E8",   # Light blue
  "CPAEC 5%" = "#52A868",    # Medium green
  "CPAEC 10%" = "#7BC892"    # Light green
)

colors <- c(
  "ODAEC 0%" = "#2E5C8A",    # Dark blue (baseline)
  "ODAEC 5%" = "#5B9BD5",    # Medium blue 
  "ODAEC 10%" = "#8BB6E8",   # Light blue
  "CPAEC 5%" = "#8E44AD",    # Medium purple
  "CPAEC 10%" = "#BB8FCE"    # Light purple
)

colors <- c(
  "ODAEC 0%" = "#87CEEB",    # Light blue (baseline)
  "ODAEC 5%" = "#5B9BD5",    # Medium blue 
  "ODAEC 10%" = "#1F4E79",   # Dark blue
  "CPAEC 5%" = "#D2B4DE",    # Light purple
  "CPAEC 10%" = "#6A1B9A"    # Dark purple
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

ggsave(here("results/21_barplot_cells.jpeg"), p_prop, width = 7, height = 8, dpi = 300)
ggsave(here("results/21_barplot_cells.pdf"), p_prop, width = 7, height = 8)


finalizeSkript()
