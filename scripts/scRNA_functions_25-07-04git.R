

checkSCTslots = function(seuratobject) {
  lapply(X = SCTResults(object = seuratobject[["SCT"]],
                        slot = "cell.attributes"), FUN = function(x) head(x,1))

}



jpeg2 = function(fn, ...) {jpeg(fn, quality = 100, unit ="in", res = 150, ...)}

doMarkerDotPlot <- function(seurat, marker_groups_peter,grouping_factor =  "predicted.celltype.l1_1st", clustermethod ="complete", inverse_x_axis = FALSE, custom_x_axis = NULL, custom_y_axis = NULL) {
  # seurat = subcluster_endo_chosen
  # grouping_factor =  "predicted.celltype.l1.5"
  
  # seurat = seurat10x
  # marker_groups_peter = marker_groups_heart
  # grouping_factor =  "predicted.celltype.l2_subendo_major"

  # clustermethod ="complete"; inverse_x_axis = FALSE; custom_x_axis = NULL; custom_y_axis = NULL
  # seurat = blood
  # marker_groups_peter = marker_groups1
  # grouping_factor = 'pp_bloodmarker_221218v2'
  # 
  # 
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(tidyr)
  library(data.table)
  library(ggplot2)
  library(toolboxH)
  
  marker_groups2 = copy(marker_groups_peter)
  marker_groups2[, highmarker := markers]
  marker_groups2[, celltype2 := celltype]
  marker_groups = copy(marker_groups2)
  
  
  marker_groups2[, highmarker := str_trim(highmarker)]
  marker_groups2[,.N, celltype]
  marker_groups2[,.N, celltype2] %>% data.frame()
  
  qlist64 = toolboxH::venn2(rownames(seurat), marker_groups2$highmarker, plotte = F)
  
  # Create the base dotplot
  seurat[[grouping_factor]] %>% mytable()
  message("Not in seurat: ", marker_groups2[!highmarker %in% rownames(seurat), paste(unique(highmarker), collapse = "\n")])
  print(marker_groups2[!highmarker %in% rownames(seurat), ])
  # resi= c()
  # resi$seurat = seurat
  # resi$qlist64 = qlist64
  # return(resi)
  # DotPlot(seurat, features =  "Pax5", group.by = grouping_factor)
  
  DotPlot(seurat, features =  unique(qlist64$q1), group.by = grouping_factor)
  plot_data <- DotPlot(
    seurat,
    group.by = grouping_factor,
    # seurat_human2,
    #  group.by ="predicted.celltype.l2_1st",
    features =   unique(qlist64$q1),
    scale.by = "radius",
    scale = TRUE)
 
# print(plot_data$data)
  # Extract the plot data
  
# plot_df <-                           plot_data$data %>% as.data.table(keep.rownames = T)
plot_df <- data.table::as.data.table(plot_data$data, keep.rownames = TRUE)
plot_df
  plot_df[,.N, id]

  # Add grouping information
  marker_groups2[, .(highmarker, celltype2)] %>% unique() %>% .[ toolboxH::allDuplicatedEntries(highmarker)]
  plot_df2 = merge(plot_df, marker_groups2[, .(highmarker, group = celltype)] %>% unique(), by.x = 'features.plot', by.y = 'highmarker', all.x = T,allow.cartesian=TRUE )


  plot_df3 = plot_df2[is.na(features.plot)==F]
  plot_df3[is.na(avg.exp.scaled) &  avg.exp ==0, avg.exp.scaled := 0]

  # Create the modified plot

  #------------------------------------------------..
  # 2) Aggregate: e.g., take mean of 'avg.exp.scaled' by (features.plot, group)
  #    Then pivot to wide format (features in rows, groups in columns)
  #------------------------------------------------..

  # Compute mean of avg.exp.scaled by (group, id)
  agg_dt <- plot_df3[
    , .(mean_value = mean(avg.exp.scaled, na.rm = TRUE)),
    by = .(group, id)
  ]

  # Pivot: rows = group, columns = id, values = mean_value
  wide_dt <- dcast(
    agg_dt,
    formula = group ~ id,
    value.var = "mean_value",
    fill = 0  # fill any missing combos with 0 (or NA if you prefer)
  )

  #------------------------------------------------.
  # 3) Row clustering (cluster group)
  #------------------------------------------------.
  row_mat <- as.matrix(wide_dt[, -1])  # omit the 'group' column
  rownames(row_mat) =wide_dt$group
  hh(row_mat,11)

  row_dist <- dist(row_mat, method = "euclidean" )      # distance (e.g. "euclidean")
  

  row_hclust <- hclust(row_dist, method = clustermethod)  # linkage (e.g. "complete")
  plot(row_hclust)
  row_order <- row_hclust$order



  # Extract the new feature order from the dendrogram
  new_group_levels <- rownames(row_mat)[row_order]

  #------------------------------------------------.
  # 4) Column clustering (i.e., cluster groups)
  #------------------------------------------------.
  # Transpose the matrix so columns become rows
  col_mat <- t(row_mat)
  hh(col_mat,11)
  col_dist <- dist(col_mat, method = "euclidean")
  col_hclust <- hclust(col_dist, method = "complete")
  plot(col_hclust)
  col_order <- col_hclust$order

  new_id_levels <- rownames(col_mat)[col_order]


  #------------------------------------------------.
  # 5) For clustering genes within group Aggregate: e.g., take mean of 'avg.exp.scaled' by (features.plot, group)
  #------------------------------------------------.

  # Compute mean of avg.exp.scaled by (features.plot, id)
  agg_dt2 <- plot_df3[
    , .(mean_value = mean(avg.exp.scaled, na.rm = TRUE)),
    by = .(features.plot, id)
  ]

  # Pivot: rows = group, columns = id, values = mean_value
  wide_dt2 <- dcast(
    agg_dt2,
    formula = features.plot ~ id,
    value.var = "mean_value",
    fill = 0  # fill any missing combos with 0 (or NA if you prefer)
  )

  #------------------------------------------------.
  # 3) Row clustering (cluster group)
  #------------------------------------------------.
  row_mat2 <- as.matrix(wide_dt2[, -1])  # omit the 'group' column
  rownames(row_mat2) =wide_dt2$features.plot
  hh(row_mat2,11)

  row_dist2 <- dist(row_mat2, method = "euclidean")      # distance (e.g. "euclidean")
  row_hclust2 <- hclust(row_dist2, method = "complete")  # linkage (e.g. "complete")
  plot(row_hclust2)
  row_order2 <- row_hclust2$order


  # Extract the new feature order from the dendrogram
  new_features.plot_levels <- rownames(row_mat2)[row_order2]



  #------------------------------------------------.
  # 5) Reorder factor levels in the original data.table
  #------------------------------------------------.

  if(inverse_x_axis) {
    plot_df3[, group := factor(group, levels = (new_group_levels %>% rev()))]


  } else plot_df3[, group := factor(group, levels = new_group_levels )]

  if(is.null(custom_x_axis)==FALSE) {
    trier = try(stopifnot(all(plot_df3$group %in% custom_x_axis)))
    if(class(trier)=="try-error") {
      print(toolboxH::venn2(plot_df3$group , custom_x_axis))
      stop("no perfect overlap")
    }
    plot_df3[, group := factor(group, levels = custom_x_axis )]

  }


  plot_df3[, features.plot := factor(features.plot, levels = new_features.plot_levels)]

  plot_df3[, id := factor(id, levels = new_id_levels)]

  if(is.null(custom_y_axis)==FALSE) {
    triery = try(stopifnot(all(plot_df3$id %in% custom_y_axis)))
    if(class(triery)=="try-error") {
      print(toolboxH::venn2(plot_df3$id , custom_y_axis))
      stop("no perfect overlap")
    }
    plot_df3[, id := factor(id, levels = custom_y_axis)]
  }
  message("====================================\nCurrent order X axis:\n")
  print(dput(levels(plot_df3$group)))
  message("====================================\nCurrent order Y axis:\n")
  print(dput(levels(plot_df3$id)))
  #------------------------------------------------.
  # 6) Plot with ggplot2 in the new clustered order
  #------------------------------------------------.
  markerplot <- ggplot(plot_df3, aes(x = id, y = features.plot)) +
    geom_point(aes(size = pct.exp, color = avg.exp.scaled)) +
    scale_color_gradient(low = "lightgrey", high = "blue") +
    theme_classic() +
    theme(
      axis.text.x      = element_text(angle = 45, hjust = 1),
      strip.background = element_blank(),
      # strip.text       = element_text(face = "bold", hjust = 0),
      strip.text.x = element_text(angle = 90, hjust = 0,face = "bold", vjust = 0.5),
      panel.spacing    = unit(1, "lines")
    ) +
    facet_grid(. ~ group, switch = "x", scales = "free", space = "free") +
    labs(x = "Cell Type", y = "Genes",
         color = "Average Expression", size = "Percent Expressed") +
    coord_flip() +
    ggtitle(deparse(substitute(seurat)), subtitle = deparse(substitute(grouping_factor)))

  print(markerplot)

  
}


ggplotSankey <- function(input, sort_by_frequency =T, customcolors = NULL) {
  require(ggsankey)
  require(ggrepel)
  
  require(dplyr)
  sankey_gg_data = ggsankey::make_long(input, names(input))
  
  
  sankey_gg_data_dt = as.data.table(sankey_gg_data, keep.rownames = F)
  ordered = sankey_gg_data_dt[,.N, .(node, next_node)][order(is.na(next_node),-N, na.last = T)  ]

  
  if(sort_by_frequency ==T){
    sankey_gg_data$next_node = factor(sankey_gg_data$next_node, levels = unique(ordered$next_node) )
    sankey_gg_data$node = factor(sankey_gg_data$node, levels = unique(ordered$node))
  }
  p_sankey = ggplot(sankey_gg_data, aes(x = x, next_x = next_x,
                                        node = node, next_node = next_node,
                                        fill = factor(node), label = node)) +
    geom_sankey(flow.alpha = 0.6, node.color = alpha("grey55", 0.2), alpha = 0.4) +
    geom_sankey_text(size = 4, color = "black") +
    theme_sankey(base_size = 18) +
    labs(x = NULL) +
    theme(legend.position = "none",
          plot.title = element_text(hjust = .5),
          axis.text.x = element_text(color = "black"))
  
  if(is.null(customcolors)) p_sankey else p_sankey + scale_fill_manual(values = customcolors)  
  
  
}




plot_ggSankey_hk <- function(data ) {
stop("Function plot_ggSankey_hk() deprecated. Use function ggplotSankey() instead")
}

plotQC = function(seuratobject,plot_correlation=T,
                  splitBy = "run10x",
                  qccols = c('nCount_RNA', 'nFeature_RNA',  'pct.mito','pct.ribo'), plotlog = TRUE, dotplotcolor  = "pct.mito") {


  plotdata = seuratobject@meta.data %>% dplyr::select(all_of(c(splitBy,qccols[qccols %in% names(seuratobject@meta.data)])))
  setDT(plotdata)
  plotdata[, separator := get(splitBy) ]
  plotdatam=melt(plotdata, id.vars ='separator', measure.vars = qccols, variable.name = "metric" )

  p1 = ggplot( plotdatam,
               aes(x=separator,y=value,fill=separator)) +
    geom_violin() +
    facet_wrap(~metric,ncol=length(qccols),scales='free_y') +
    theme(axis.text.x=element_blank())
  if(plotlog) p1 = p1 + scale_y_log10(breaks = log_breaks(10))

  if(dotplotcolor %in% names(seuratobject@meta.data) & plot_correlation==T) {
    p4 = ggplot(plotdata, aes_string(x='nCount_RNA' ,y='nFeature_RNA', col = paste0("-", dotplotcolor))) +
      geom_point(alpha = 0.3) +
      facet_wrap(~separator,ncol=4,scales='free_y') +
      scale_color_gradient2_tableau()
    if(plotlog) p4 = p4 + scale_y_log10(breaks = log_breaks(10)) + scale_x_log10(breaks = log_breaks(10)) + theme(axis.text.x=element_text(angle = 45, hjust = 1))
  
    return((p1+p4 + patchwork::plot_layout(ncol = 1)) )
    # plot(p1+p4 + patchwork::plot_layout(ncol = 1)) 
  } else if( plot_correlation==T)
  {


    p4 = ggplot(plotdata, aes(x=nCount_RNA ,y=nFeature_RNA)) +
      geom_point(alpha = 0.3) +
      facet_wrap(~separator,ncol=4,scales='free_y')
    if(plotlog) p4 = p4 + scale_y_log10(breaks = log_breaks(10)) + scale_x_log10(breaks = log_breaks(10)) +theme(axis.text.x=element_text(angle = 45, hjust = 1))

    return((p1+p4 + patchwork::plot_layout(ncol = 1)) )
  } else {
    return(p1)
  }
}




doSCTransform = function(seuratobject, experimentIDcol = "run10x", variable.features.n = 4000, ...) {
  # seuratobject = mouse;experimentIDcol = "run10x"

  message(Sys.time(), "...Running Seurat SCT separately on following cells:")
  mytable(seuratobject[[experimentIDcol]])

  n_experiments = uniqueN(names(table(seuratobject[[experimentIDcol]])))

  if(n_experiments ==1) seuratobject_list = list(seuratobject) else seuratobject_list <- SplitObject(object = seuratobject, split.by = experimentIDcol)
  seuratobject_list


  for (i in 1:length(seuratobject_list)) {
    # i=1
    message(Sys.time(), "...SCT on ", i, "\n-------------------------------------------------------")
    seuratobject_list[[i]] <- SCTransform(seuratobject_list[[i]], verbose = T,variable.features.n = variable.features.n, ... ) # TODO include "G2M.Score", "S.Score" regression if necessary
  }
  print(table(names(warnings() ))) #
  message(Sys.time(), '...Notes about warning:\n"iteration limit reached"\nChristophH commented on 22 May 2019 - These warnings are showing that there are some genes for which it is hard to reliably estimate theta (presumably because of very few non-zero observations). Usually we donÂ´t worry about these warnings too much, since we regularize the parameters in a later step, thus averaging out uncertainty of individual gene parameters. https://github.com/ChristophH/sctransform/issues/25')


  list.features <- SelectIntegrationFeatures(object.list = seuratobject_list, nfeatures = 4000)

  if(length(seuratobject_list)>1) {
    seuratobject <- merge(seuratobject_list[[1]],
                          y = seuratobject_list[2:length(seuratobject_list)],
                          # project = "seuratobject",
                          merge.data = TRUE)

  } else seuratobject = seuratobject_list[[1]]

  VariableFeatures(seuratobject) <- list.features
  seuratobject
}

plot3clusterings = function(seuratobject, clustervar1 = 'SCT_snn_res.0.2', clustervar2 = 'SCT_snn_res.0.4', clustervar3 = 'SCT_snn_res.0.8') {
  # clustervar1 = 'SCT_snn_res.0.2'; clustervar2 = 'SCT_snn_res.0.4'; clustervar3 = 'SCT_snn_res.0.8'
  seuratobject_name = deparse(substitute(seuratobject))
  p1 = DimPlot(seuratobject, group.by = clustervar1, label = T)
  p2 = DimPlot(seuratobject, group.by = clustervar2, label = T)
  p3 = DimPlot(seuratobject, group.by = clustervar3, label = T)

  plot(p1+p2+p3 + patchwork::plot_annotation(title = seuratobject_name))
}







colplotLikeExcel = function (plotdat, mycolors = c("dodgerblue2", "white", "red"),
                             lowest_colorval = "minimum", middle_colorval = "median",
                             highest_colorval = "maximum", xlabel = "", ylabel = "",
                             x_axis_pos = "top", myround = 0, userdefined_labels = NULL,
                             row_names = NULL, sort_via_value = T)
{
  plotdat_m = reshape2::melt(plotdat) %>% data.table()
  var1name = names(plotdat_m)[1]
  var2name = names(plotdat_m)[2]

  if(sort_via_value==T) plotdat_m = plotdat_m[order(value, decreasing = T)]
  setnames(plotdat_m, c(var1name, var2name), c("Var1", "Var2"))

  if (sort_via_value==T) {
    plotdat_m$Var1 = factor(plotdat_m$Var1, levels = rev(unique(plotdat_m$Var1)))
  }  else plotdat_m$Var1 = factor(plotdat_m$Var1, levels = rev(unique(rownames(plotdat))))

  sort(unique(plotdat_m$Var1))


  if (sort_via_value==T) {
    plotdat_m$Var2 = factor(plotdat_m$Var2, levels = unique(plotdat_m$Var2))
  }  else plotdat_m$Var2 = factor(plotdat_m$Var2, levels = unique(colnames(plotdat)))

  if (is.numeric(plotdat_m$value) == F) {
    stop("Need numeric matrix as `plotdat` argument in order to know coloring according to provided numeric data! Stoping.\nConsider providing a userdefined matrix with names via parameter `userdefined_labels` in addition to providing numeric values via parameter `plotdat`.")}
  plotdat_m$value = round(plotdat_m$value, myround)
  if (lowest_colorval == "minimum") {
    lowest_colorval = min(plotdat_m$value, na.rm = T)
  }  else lowest_colorval = lowest_colorval
  if (middle_colorval == "median") {
    middle_colorval = stats::median(plotdat_m$value, na.rm = T)
  }  else middle_colorval = middle_colorval
  if (highest_colorval == "maximum") {
    highest_colorval = max(plotdat_m$value, na.rm = T)
  }  else highest_colorval = highest_colorval

  if (is.null(userdefined_labels)) {
    plot1 = ggplot2::ggplot(plotdat_m, ggplot2::aes(Var2,
                                                    Var1, label = value)) + ggplot2::geom_tile(ggplot2::aes(fill = value),
                                                                                               colour = "white") + ggplot2::scale_fill_gradientn(colours = mycolors,
                                                                                                                                                 values = scales::rescale(c(lowest_colorval, middle_colorval,
                                                                                                                                                                            highest_colorval)), guide = FALSE) + ggplot2::geom_text(show.legend = FALSE) +
      ggplot2::scale_x_discrete(position = x_axis_pos) +
      ggplot2::xlab(xlabel) + ggplot2::ylab(ylabel)
  }  else     {

    beschriftdat = as.matrix(userdefined_labels)
    beschriftdat = beschriftdat[rev(rownames(beschriftdat)),
    ]
    stopifnot(identical(dim(plotdat), dim(beschriftdat)))
    stopifnot(identical(rownames(plotdat), rownames(beschriftdat)))
    stopifnot(identical(colnames(plotdat), colnames(beschriftdat)))
    beschriftdat_m = data.table::melt(beschriftdat)
    plot1 = ggplot2::ggplot(plotdat_m, ggplot2::aes(Var2,
                                                    Var1)) + ggplot2::geom_tile(ggplot2::aes(fill = value),
                                                                                colour = "white") + ggplot2::scale_fill_gradientn(colours = mycolors,
                                                                                                                                  values = scales::rescale(c(lowest_colorval, middle_colorval,
                                                                                                                                                             highest_colorval)), guide = "none") + ggplot2::geom_text(label = beschriftdat_m$value,
                                                                                                                                                                                                                      show.legend = FALSE) + ggplot2::scale_x_discrete(position = x_axis_pos) +
      ggplot2::xlab(xlabel) + ggplot2::ylab(ylabel)
  }
  plot1 + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90,
                                                             hjust = 0)) + xlab(var2name) + ylab(var1name)
}



# alternative scalings
scale_y_log1p <- function(..., trans = log1p_trans()) {
  scale_y_continuous(trans = trans, ...)
}

log1p_trans <- function() {
  trans_new(
    "log1p",
    transform = function(x) log10(x + 1),
    inverse = function(x) 10^x - 1
  )
}

log1p_breaks <- function (n = 5, base = 10)
{
  scales:::force_all(n, base)
  n_default <- n
  function(x, n = n_default) {
    raw_rng <- suppressWarnings(range(x, na.rm = TRUE))
    if (any(!is.finite(raw_rng))) {
      return(numeric())
    }
    rng <- log(raw_rng+1, base = base)
    min <- floor(rng[1])
    max <- ceiling(rng[2])
    if (max == min) {
      return(base^min)
    }
    by <- floor((max - min)/n) + 1
    breaks <- base^seq(min, max, by = by)
    relevant_breaks <- base^rng[1] <= breaks & breaks <=
      base^rng[2]
    if (sum(relevant_breaks) >= (n - 2)) {
      return(breaks)
    }
    while (by > 1) {
      by <- by - 1
      breaks <- base^seq(min, max, by = by)
      relevant_breaks <- base^rng[1] <= breaks & breaks <=
        base^rng[2]
      if (sum(relevant_breaks) >= (n - 2)) {
        return(breaks)
      }
    }
    scales:::log_sub_breaks(rng, n = n, base = base)
  }
}


pct_counts_in_top_N_genes_SCT <- function(seurat_obj, top_N = 100, slot = "SCT") {
  # Check if the SCT slot exists


  # Get the SCT assay data matrix
  if(identical(slot, "SCT")) {
    data_matrix <- seurat_obj@assays$SCT@data
    if (!"SCT" %in% names(seurat_obj@assays)) {
      stop("The Seurat object does not have an SCT slot. Please make sure you have performed SCTransform on the object.")
    }
    } else
    if(identical(slot, "RNA")) {
      data_matrix <- seurat_obj@assays$RNA@data
      if (!"RNA" %in% names(seurat_obj@assays)) {
        stop("The Seurat object does not have an RNA slot.")
      }
      }else
      stop("Parameter slot must be either `SCT` or `RNA`")

  # hh(data_matrix)
  # Calculate the top N genes for each cell
  top_N_genes <- apply(data_matrix, 2, function(x) {
    top_genes <- sort(x, decreasing = TRUE)[1:top_N]
    sum_top_genes <- sum(top_genes)
    return(sum_top_genes)
  })

  # Calculate the total counts for each cell
  nCount_RNA <- colSums(data_matrix)

  # Calculate the percentage of counts in the top N genes
  pct_counts <- top_N_genes / nCount_RNA * 100

  return(pct_counts)
}



makeInputUpset = function(plotlist) {
  # better in toolboxH or genstathelpR
  library(data.table)
  if(length(names(plotlist))==0) {
    message("no names for list entries found -  providing standardnames  `list1...")
    names_list = paste0("list", seq(along = plotlist))
    names(plotlist) = names_list
  } else names_list = names(plotlist)

  plotlist2 = lapply(names_list, function(myname) {
    data.table(variable = myname,value = plotlist[[myname]])
  }
  )

  plotlist3 = rbindlist(plotlist2)

  plotlist4 = dcast.data.table(plotlist3, value ~ variable, fun.aggregate = function(x) as.numeric(length(x)>0 ))
  plotlist4

}

# NogpaletteReihe----
NogpaletteReihe <-  c("#CB769E", "#DE639A", "#A85C85", "#0081AF", "#4F6D7A", "#7C6A0A", "#368F8B", "#246A73", "#5CC1BC", "#62C370", "#F7C548", "#F97E44", "#FB3640", "#B7245C", "#0D3B66", "#3E2F5B", "#B2675E", "#644536")


get_full_gene_names <- function(gene_symbols, organism = "mouse") {
  
  input = data.table(gene_symbols)
  input[, initialorder := .I]
  
  
  if(organism == "mouse") {
    
    library(org.Mm.eg.db)
    results <- select(org.Mm.eg.db, keys = gene_symbols, columns = c("SYMBOL", "GENENAME"), keytype = "SYMBOL")%>% data.table(keep.rownames = F)
    results2 = results[,.(GENENAME=paste(GENENAME, collapse = " /// alternative_fullname:,")), by = SYMBOL]
    
    input[, GENENAME := results2[match_hk(input$gene_symbols, results2$SYMBOL), GENENAME]]
    
    return(input$GENENAME)
    
  } else if(organism == "human") {
    library(org.Hs.eg.db)
    
    results <- select(org.Hs.eg.db, keys = unique(gene_symbols), columns = c("SYMBOL", "GENENAME"), keytype = "SYMBOL") %>% data.table(keep.rownames = F)
    
    results2 = results[,.(GENENAME=paste(GENENAME, collapse = " /// alternative_fullname:,")), by = SYMBOL]
    
    input[, GENENAME := results2[match_hk(input$gene_symbols, results2$SYMBOL), GENENAME]]
    
    return(input$GENENAME)
  } else {
    stop("Unsupported organism. Please use 'mouse' or 'human'.")
  }
  
}


gptcelltype_hk <- function(input, tissuename=NULL, model='gpt-4', topgenenumber = 10, 
                           include_downreg = TRUE) {
  # adapted from https://github.com/Winnie09/GPTCelltype
  
  OPENAI_API_KEY <- Sys.getenv("OPENAI_API_KEY")
  if (OPENAI_API_KEY == "") {
    print("Note: OpenAI API key not found: returning the prompt itself.")
    API.flag <- 0
  } else {
    API.flag <- 1
  }
  
  if (class(input)=='list') {
    input <- sapply(input,paste,collapse=',')
  } else {
    # Create separate lists for upregulated and downregulated genes
    up_input <- input[input$avg_log2FC > 0,,drop=FALSE]
    up_genes <- tapply(up_input$gene, list(up_input$cluster), 
                       function(i) paste0(i[1:min(topgenenumber, length(i))], collapse=','))
    
    # Process downregulated genes if requested
    if (include_downreg) {
      down_input <- input[input$avg_log2FC < 0,,drop=FALSE]
      down_genes <- tapply(down_input$gene, list(down_input$cluster), 
                           function(i) paste0(i[1:min(topgenenumber, length(i))], collapse=','))
      
      # Create formatted input for each cluster
      formatted_input <- character(0)
      all_clusters <- unique(c(names(up_genes), names(down_genes)))
      
      for (cluster in all_clusters) {
        up_str <- if (cluster %in% names(up_genes)) up_genes[cluster] else ""
        down_str <- if (cluster %in% names(down_genes)) down_genes[cluster] else ""
        formatted_input[cluster] <- paste0("UP[", up_str, "] DOWN[", down_str, "]")
      }
    } else {
      # Use only upregulated genes if downregulated genes are not requested
      formatted_input <- paste0("UP[", up_genes, "] DOWN[]")
      names(formatted_input) <- names(up_genes)
    }
    
    input <- formatted_input
  }
  
  # Create the prompt header
  prompt_header <- paste0('Identify cell types of ', tissuename, ' cells using the following markers for each row representin a cluster of cells in a Seurat scRNA-seq analysis. ',
                          'Each row representing the cluster starts with a certain number as ID for this cluster and a colon, followed by genes differentially expressed in this cluster vs. the others, provided in two categories:\n',
                          '1. Upregulated genes are shown in brackets after "UP" (e.g., UP[Gene1,Gene2,Gene3])\n',
                          '2. Downregulated genes are shown in brackets after "DOWN" (e.g., DOWN[Gene4,Gene5,Gene6])\n\n',
                          'Only provide the cell type name. Do show cluster number and a colon before the name. If you do not find a celltype, return "Unknown". Return the same number of rows as the number of clusters provided as input\n ',
                          'Some clusters can be a mixture of multiple cell types.\n\n')
  
  # Generate full prompt with all genes included
  full_prompt <- paste0(prompt_header, paste0(names(input), ': ', unlist(input), collapse = '\n'))
  
  # Always print the full prompt for copying to online LLM
  print("Complete prompt (copy this to test with an online LLM):")
  cat(full_prompt)
  
  if (!API.flag){
    return(full_prompt)
  } else {
    print("Note: OpenAI API key found: returning the cell type annotations.")
    
    cutnum <- ceiling(length(input)/50)
    if (cutnum > 1) {
      cid <- as.numeric(cut(1:length(input),cutnum))	
      print(paste0("Note: Due to the number of clusters (", length(input), "), the API request will be split into ", 
                   cutnum, " batches."))
    } else {
      cid <- rep(1,length(input))
    }
    
    allres <- list()
    
    for (i in 1:cutnum) {
      # i=1
      id <- which(cid==i)
      
      # Generate batch-specific prompt
      batch_prompt <- paste0(prompt_header, paste0(names(input)[id], ': ', input[id], collapse = '\n'))
      
      # For each batch, print which clusters are included
      cat(paste0("Batch ", i, " includes clusters: ", paste(names(input)[id], collapse=", ")))
      
      # Use tryCatch to handle potential errors without looping
      result <- tryCatch({
        k <- openai::create_chat_completion(
          model = model,
          message = list(list("role" = "user", "content" = batch_prompt))
        )
        k_global <<- k
        res <- strsplit(k$choices[,'message.content'],'\n')[[1]]
        
        # Check if we have the right number of results
        if (length(res) != length(id)) {
          warning(paste("Expected", length(id), "results, but got", length(res), 
                        "- the model response may be incomplete or incorrectly formatted."), immediate. = TRUE )
        }
        
        names(res) <- names(input)[id]
        res
      }, error = function(e) {
        # Print the error message and return NULL instead of retrying
        print(paste("Error in API call for batch", i, ":", e$message))
        print("Suggestion: Try reducing the number of genes per cluster or clusters per batch.")
        return(res)
      })
      
      if (!is.null(result)) {
        allres[[i]] <- result
      } else {
        print(paste("Batch", i, "failed. Continuing with next batch if available."))
      }
    }
    
    # Combine results from successful batches
    if (length(allres) > 0) {
      final_results <- gsub(',$', '', unlist(allres))
      print('Note: It is always recommended to check the results returned by GPT-4 in case of\n AI hallucination, before going to down-stream analysis.')
      return(final_results)
    } else {
      print("All API calls failed. Please check error messages above.")
      return(NULL)
    }
  }
}

customDimPlot <- function(seurat_obj, group_by = "seurat_clusters", 
                          pt.size = 0.5, label.size = 4, pt.alpha = 0.3,
                          reduction = "umap", max.overlaps = Inf,
                          show_sticks = 1, label_darkening_alpha = 0.2, colors = NULL) {
  
  # Custom DimPlot with ggrepel labels
  
  # Usage examples:
  # Light background: customDimPlot(seurat_obj, label_darkening_alpha = 0.1)
  # Medium background: customDimPlot(seurat_obj, label_darkening_alpha = 0.3) 
  # Strong background: customDimPlot(seurat_obj, label_darkening_alpha = 0.5)
  
  # p_dim = customDimPlot(subseurat, group_by = "macro_subcluster", label.size = 9 , pt.size = 1) + customDimPlot(subseurat, group_by = "celltypes_250822", label.size = 6 , pt.size = 2)
  # p_dim
  
  
  # Extract coordinates and metadata
  coords <- Embeddings(seurat_obj, reduction = reduction)
  meta_data <- seurat_obj@meta.data[[group_by]]
  
  # Create data frame
  plot_df <- data.frame(
    x = coords[, 1],
    y = coords[, 2],
    cluster = meta_data
  )
  
  # Calculate label positions (centroids)
  label_pos <- plot_df %>%
    group_by(cluster) %>%
    summarise(
      x = median(x),
      y = median(y),
      .groups = 'drop'
    )
  
  # Get colors from Seurat's default palette
  if(length(colors)==0) colors <- scales::hue_pal()(length(unique(plot_df$cluster))) else colors = colors[1:length(unique(plot_df$cluster))]
  names(colors) <- levels(as.factor(plot_df$cluster)) 
  
  # Set seed for consistent positioning
  set.seed(42)
  
  # Create custom plot with ggrepel stick controls
  p <- ggplot(plot_df, aes(x = x, y = y, color = cluster)) +
    geom_point(size = pt.size, alpha = pt.alpha) +
    # First layer: black background text for better visibility
    geom_text_repel(data = label_pos, 
                    aes(x = x, y = y, label = cluster, color = cluster),
                    size = label.size, 
                    fontface = "bold",
                    max.overlaps = max.overlaps,
                    segment.color = "black",
                    segment.size = 0.3,
                    segment.alpha = 0.6,
                    segment.linetype = 1,
                    min.segment.length = if(show_sticks) 0 else show_sticks,
                    show.legend = FALSE,
                    seed = 42) +
    
    geom_text_repel(data = label_pos, 
                    aes(x = x, y = y, label = cluster),
                    size = label.size, 
                    fontface = "bold",
                    color = "black",
                    alpha = label_darkening_alpha,  # Controllable alpha parameter
                    max.overlaps = max.overlaps,
                    segment.color = "black",
                    segment.size = 0.3,
                    segment.alpha = 0.6,
                    segment.linetype = 1,
                    min.segment.length = if(show_sticks) 0 else show_sticks,
                    show.legend = FALSE,
                    seed = 42) +
    # Second layer: colored text on top
    
    scale_color_manual(values = colors) +
    theme_void() +
    theme(legend.position = "none") +
    labs(x = paste0(reduction, "_1"), y = paste0(reduction, "_2"))
  
  return(p)
}
