# title: "scRNA-seq data analysis with Seurat"
# author: "Tingting Zhao"
# email: tingting_zhao@dfci.harvard.edu
# date: "09/06/2022"
  
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(stringr)
library(SoupX)
library(reticulate)
library(glmpca)
library(SeuratWrappers)
library(scry)
library(reticulate)
library(monocle3)
library(FlexDotPlot)
library(cowplot)
library(googlesheets4)
library(tidyverse)
library(viridis)
library(scCustomize)
library(qs)
library(gridExtra)
library(plyr)
library(circlize)
library(ComplexHeatmap)

args <- commandArgs(trailingOnly = TRUE)

pwd = args[1]
indir = args[2]
outdir = args[3]
scrubletdir = args[4]
samples = unlist(strsplit(args[5], ','))
projectName = args[6]
marker_link = args[7]
marker_sheet = args[8]
flag = args[9]
mtPattern = args[10]
rbPattern = args[11]
qc_cutoff = as.numeric(args[12])
mito_cutoff = as.numeric(args[13])
sex = unlist(strsplit(args[14], ','))
genotypes = unlist(strsplit(args[15], ','))
refdir = args[16]
scriptdir = args[17]
geneN = as.numeric(args[18])

# message("Read in marker genes")
# gsurl=marker_link
# gs4_deauth()
# markers = read_sheet(gsurl, sheet = marker_sheet) %>%
#   select(Markers) %>% as.list()

message("Set working dir")
setwd(pwd)

message("If using Leiden algorithm in FindMarkers")
# use_condaenv("r_leiden", required=TRUE)
# py_config()

message("Load clustered dotplot function")
# scCustermize code updated by developer, https://github.com/samuel-marsh/scCustomize/issues/27
Clustered_DotPlot_relabel <- function(
  seurat_object,
  features,
  new_row_labels,
  colors_use_exp = viridis_plasma_dark_high,
  exp_color_min = -2,
  exp_color_middle = NULL,
  exp_color_max = 2,
  print_exp_quantiles = FALSE,
  colors_use_idents = NULL,
  x_lab_rotate = TRUE,
  k = 1,
  row_km_repeats = 1000,
  column_km_repeats = 1000,
  row_label_size = 8,
  raster = FALSE,
  plot_km_elbow = TRUE,
  elbow_kmax = NULL,
  assay = NULL,
  group.by = NULL,
  idents = NULL,
  show_parent_dend_line = TRUE,
  ggplot_default_colors = FALSE,
  seed = 123
) {
  # Check for packages
  ComplexHeatmap_check <- PackageCheck("ComplexHeatmap", error = FALSE)
  if (!ComplexHeatmap_check[1]) {
    stop(
      "Please install the ComplexHeatmap package to use Clustered_DotPlot",
      "\nThis can be accomplished with the following commands: ",
      "\n----------------------------------------",
      "\ninstall.packages('BiocManager')",
      "\nBiocManager::install('ComplexHeatmap')",
      "\n----------------------------------------",
      call. = FALSE
    )
  }
  
  # Check Seurat
  scCustomize:::Is_Seurat(seurat_object = seurat_object)
  
  # Check unique features
  features_unique <- unique(x = features)
  
  if (length(x = features_unique) != length(x = features)) {
    warning("Feature list contains duplicates, making unique.")
  }
  
  # Check exp min/max set correctly
  if (!exp_color_min < exp_color_max) {
    stop("The value for 'exp_color_min': ", exp_color_min, ", must be less than the value for 'exp_color_max': ", exp_color_max, ".")
  }
  
  # Get DotPlot data
  seurat_plot <- DotPlot(object = seurat_object, features = features_unique, assay = assay, group.by = group.by, scale = TRUE, idents = idents, col.min = NULL, col.max = NULL)
  
  data <- seurat_plot$data
  
  # Get expression data
  exp_mat <- data %>%
    select(-pct.exp, -avg.exp) %>%
    pivot_wider(names_from = id, values_from = avg.exp.scaled) %>%
    as.data.frame()
  
  row.names(x = exp_mat) <- exp_mat$features.plot
  
  # Check NAs if idents
  if (!is.null(x = idents)) {
    # Find NA features and print warning
    excluded_features <- exp_mat[rowSums(is.na(x = exp_mat)) > 0,] %>%
      rownames()
    warning("The following features were removed as there is no scaled expression present in subset (`idents`) of object provided: ", glue_collapse_scCustom(input_string = excluded_features, and = TRUE), ".")
    
    # Extract good features
    good_features <- rownames(exp_mat)
    
    # Remove rows with NAs
    exp_mat <- exp_mat %>%
      filter(features.plot %in% good_features)
  }
  
  exp_mat <- exp_mat[,-1] %>%
    as.matrix()
  
  # Get percent expressed data
  percent_mat <- data %>%
    select(-avg.exp, -avg.exp.scaled) %>%
    pivot_wider(names_from = id, values_from = pct.exp) %>%
    as.data.frame()
  
  row.names(x = percent_mat) <- percent_mat$features.plot
  
  # Subset dataframe for NAs if idents so that exp_mat and percent_mat match
  if (!is.null(x = idents)) {
    percent_mat <- percent_mat %>%
      filter(features.plot %in% good_features)
  }
  
  percent_mat <- percent_mat[,-1] %>%
    as.matrix()
  
  # print quantiles
  if (print_exp_quantiles) {
    message("Quantiles of gene expression data are:")
    print(quantile(exp_mat, c(0.1, 0.5, 0.9, 0.99)))
  }
  
  # set assay (if null set to active assay)
  assay <- assay %||% DefaultAssay(object = seurat_object)
  
  # Set default color palette based on number of levels being plotted
  if (is.null(x = group.by)) {
    group_by_length <- length(x = unique(x = seurat_object@active.ident))
  } else {
    group_by_length <- length(x = unique(x = seurat_object@meta.data[[group.by]]))
  }
  
  # Check colors use vs. ggplot2 color scale
  if (!is.null(x = colors_use_idents) && ggplot_default_colors) {
    stop("Cannot provide both custom palette to `colors_use` and specify `ggplot_default_colors = TRUE`.")
  }
  if (is.null(x = colors_use_idents)) {
    # set default plot colors
    if (is.null(x = colors_use_idents)) {
      colors_use_idents <- scCustomize_Palette(num_groups = group_by_length, ggplot_default_colors = ggplot_default_colors, color_seed = color_seed)
    }
  }
  
  # Pull Annotation and change colors to ComplexHeatmap compatible format
  Identity <- colnames(exp_mat)
  
  identity_colors <- DiscretePalette_scCustomize(num_colors = length(Identity), palette = "polychrome", shuffle_pal = F)
  names(identity_colors) <- Identity
  identity_colors_list <- list(Identity = identity_colors)
  
  # Create identity annotation
  column_ha <- ComplexHeatmap::HeatmapAnnotation(Identity = Identity,
                                                 col =  identity_colors_list,
                                                 na_col = "grey",
                                                 name = "Identity"
  )
  
  # Set middle of color scale if not specified
  if (is.null(x = exp_color_middle)) {
    exp_color_middle <- scCustomize:::Middle_Number(min = exp_color_min, max = exp_color_max)
  }
  
  palette_length <- length(colors_use_exp)
  palette_middle <- scCustomize:::Middle_Number(min = 0, max = palette_length)
  
  # Create palette
  col_fun = colorRamp2(c(exp_color_min, exp_color_middle, exp_color_max), colors_use_exp[c(1,palette_middle, palette_length)])
  
  # Calculate and plot Elbow
  if (plot_km_elbow) {
    # if elbow_kmax not NULL check it is usable
    if (!is.null(x = elbow_kmax) && elbow_kmax > (nrow(x = exp_mat) - 1)) {
      elbow_kmax <- nrow(x = exp_mat) - 1
      warning("The value provided for 'elbow_kmax' is too large.  Changing to (length(x = features)-1): ", elbow_kmax)
    }
    
    # if elbow_kmax is NULL set value based on input feature list
    if (is.null(x = elbow_kmax)) {
      # set to (length(x = features)-1) if less than 21 features OR to 20 if greater than 21 features
      if (nrow(x = exp_mat) > 21) {
        elbow_kmax <- 20
      } else {
        elbow_kmax <- nrow(x = exp_mat) - 1
      }
    }
    
    km_elbow_plot <- scCustomize:::kMeans_Elbow(data = exp_mat, k_max = elbow_kmax)
  }
  
  # prep heatmap
  if (raster) {
    layer_fun = function(j, i, x, y, w, h, fill) {
      grid.rect(x = x, y = y, width = w, height = h,
                gp = gpar(col = NA, fill = NA))
      grid.circle(x=x,y=y,r= sqrt(ComplexHeatmap::pindex(percent_mat, i, j)/100)  * unit(2, "mm"),
                  gp = gpar(fill = col_fun(ComplexHeatmap::pindex(exp_mat, i, j)), col = NA))
    }
  } else {
    cell_fun = function(j, i, x, y, w, h, fill) {
      grid.rect(x = x, y = y, width = w, height = h,
                gp = gpar(col = NA, fill = NA))
      grid.circle(x=x,y=y,r= sqrt(percent_mat[i, j]/100) * unit(2, "mm"),
                  gp = gpar(fill = col_fun(exp_mat[i, j]), col = NA))
    }
  }
  
  # Create legend for point size
  lgd_list = list(
    ComplexHeatmap::Legend(labels = c(0.25,0.5,0.75,1), title = "Percent Expressing",
                           graphics = list(
                             function(x, y, w, h) grid.circle(x = x, y = y, r = sqrt(0.25) * unit(2, "mm"),
                                                              gp = gpar(fill = "black")),
                             function(x, y, w, h) grid.circle(x = x, y = y, r = sqrt(0.5) * unit(2, "mm"),
                                                              gp = gpar(fill = "black")),
                             function(x, y, w, h) grid.circle(x = x, y = y, r = sqrt(0.75) * unit(2, "mm"),
                                                              gp = gpar(fill = "black")),
                             function(x, y, w, h) grid.circle(x = x, y = y, r = 1 * unit(2, "mm"),
                                                              gp = gpar(fill = "black")))
    )
  )
  
  # Set x label roration
  if (is.numeric(x = x_lab_rotate)) {
    x_lab_rotate <- x_lab_rotate
  } else if (isTRUE(x = x_lab_rotate)) {
    x_lab_rotate <- 45
  } else {
    x_lab_rotate <- 0
  }
  
  # Create Plot
  set.seed(seed = seed)
  if (raster) {
    cluster_dot_plot <- ComplexHeatmap::Heatmap(exp_mat,
                                                heatmap_legend_param=list(title="Expression"),
                                                col=col_fun,
                                                rect_gp = gpar(type = "none"),
                                                layer_fun = layer_fun,
                                                row_names_gp = gpar(fontsize = row_label_size),
                                                row_km = k,
                                                row_km_repeats = row_km_repeats,
                                                border = "black",
                                                top_annotation = column_ha,
                                                column_km_repeats = column_km_repeats,
                                                show_parent_dend_line = show_parent_dend_line,
                                                column_names_rot = x_lab_rotate)
  } else {
    cluster_dot_plot <- ComplexHeatmap::Heatmap(exp_mat,
                                                heatmap_legend_param=list(title="Expression"),
                                                col=col_fun,
                                                rect_gp = gpar(type = "none"),
                                                cell_fun = cell_fun,
                                                row_names_gp = gpar(fontsize = row_label_size),
                                                row_km = k,
                                                row_km_repeats = row_km_repeats,
                                                border = "black",
                                                top_annotation = column_ha,
                                                column_km_repeats = column_km_repeats,
                                                show_parent_dend_line = show_parent_dend_line,
                                                column_names_rot = x_lab_rotate)
  }
  
  # Add pt.size legend & return plots
  if (plot_km_elbow) {
    return(list(km_elbow_plot, ComplexHeatmap::draw(cluster_dot_plot, annotation_legend_list = lgd_list)))
  }
  return(ComplexHeatmap::draw(cluster_dot_plot + rowAnnotation(rn= anno_text(new_row_labels)), annotation_legend_list = lgd_list))
}

message("Read in gene name table")
geneTable <- read.csv(paste0(refdir, "geneAnnotationTable.csv"), header = T, row.names = 1)

message("Import seurat and FindAllMarkers object")
scrna.list <- readRDS(paste0(outdir, "individual/", "scrna.list.seurat.", projectName, ".rds"))
scrna.markers.list <- readRDS(paste0(outdir, "individual/", "scrna.markers.list.seurat.", projectName, ".rds"))
markers.topN.list <- readRDS(paste0(outdir, "individual/", "markers.topN.list.seurat.", projectName, ".rds"))

message("DE heatmap")
pdf(file = paste0(outdir, "individual/", "top10markers.heatmap.geneSymbol.pdf"))
for (sample in samples){
  plot1 <- DoHeatmap(scrna.list[[sample]], features = markers.topN.list[[sample]]$geneSymbol, size = 2, draw.lines = T, angle = 45, hjust = 0.2) + theme(axis.text.y = element_text(size = 5)) + NoLegend() + theme(aspect.ratio=10/10)
  print(plot1 + coord_fixed())
}
dev.off()

pdf(file = paste0(outdir, "individual/", "top10markers.heatmap.geneID.pdf"))
for (sample in samples){
  plot1 <- DoHeatmap(scrna.list[[sample]], features = markers.topN.list[[sample]]$geneSymbol, size = 2, draw.lines = T, angle = 45, hjust = 0.2) + theme(axis.text.y = element_text(size = 5)) + NoLegend() + scale_y_discrete(breaks=markers.topN.list[[sample]]$geneSymbol, labels=geneTable$geneID[match(markers.topN.list[[sample]]$geneSymbol, geneTable$geneSymbol)]) + theme(aspect.ratio=10/10)
  print(plot1 + coord_fixed())
}
dev.off()

message("Top identified genes, feature plot")
for (sample in samples){
  topN <- Extract_Top_Markers(scrna.markers.list[[sample]], num_genes = geneN, named_vector = FALSE, make_unique = TRUE, gene_column = "geneSymbol")
  pdf(paste0(outdir, "individual/", sample, "_", "top10markers.geneSymbol.pdf"))
  ggp = list()
  for (marker in topN){
    ggp[[marker]]=FeaturePlot(scrna.list[[sample]], features=marker)
    print(ggp[[marker]])
  }
  dev.off()
}

for (sample in samples){
  topN <- Extract_Top_Markers(scrna.markers.list[[sample]], num_genes = geneN, named_vector = FALSE, make_unique = TRUE, gene_column = "geneSymbol")
  pdf(paste0(outdir, "individual/", sample, "_", "top10markers.geneID.pdf"))
  ggp = list()
  for (marker in topN){
    ggp[[marker]]=FeaturePlot(scrna.list[[sample]], features=marker) + ggtitle(geneTable$geneID[match(marker, geneTable$geneSymbol)])
    print(ggp[[marker]])
  }
  dev.off()
}

message("Top identified genes, dot plot")
pdf(file = paste0(outdir, "individual/", "dotplot.DEtop10.geneSymbol.pdf"), width = 20, height = 10)
for (sample in samples){
  p1 <- DotPlot_scCustom(scrna.list[[sample]], features = topN, x_lab_rotate = TRUE, colors_use = "blue")
  print(p1)
}
dev.off()

pdf(file = paste0(outdir, "individual/", "dotplot.DEtop10.geneID.pdf"), width = 20, height = 10)
for (sample in samples){
  p1 <- DotPlot_scCustom(scrna.list[[sample]], features = topN, x_lab_rotate = TRUE, colors_use = "blue") + scale_x_discrete(breaks=topN, labels=geneTable$geneID[match(topN, geneTable$geneSymbol)])
  print(p1)
}
dev.off()

pdf(paste0(outdir, "individual/", "dotplot.DEtop10.clustered.geneSymbol.pdf"), width = 10, height = 15)
p1 <- Clustered_DotPlot_relabel(scrna.list[[sample]], features = topN, x_lab_rotate = F, plot_km_elbow = FALSE, new_row_labels = topN)
print(p1)
dev.off()

pdf(paste0(outdir, "individual/", "dotplot.DEtop10.clustered.geneID.pdf"), width = 10, height = 15)
p1 <- Clustered_DotPlot_relabel(scrna.list[[sample]], features = topN, plot_km_elbow = F, new_row_labels = geneTable$geneID[match(topN, geneTable$geneSymbol)])
print(p1)
dev.off()
