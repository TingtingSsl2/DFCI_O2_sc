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
library(monocle3)
library(FlexDotPlot)
library(cowplot)
library(googlesheets4)

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
mitoCutoff = as.numeric(args[12])
sex = unlist(strsplit(args[13], ','))
genotypes = unlist(strsplit(args[14], ','))
refdir = args[15]
scriptdir = args[16]

message("Samples are:")
print(samples)

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

message("Step 2: Pre-processing")
message("Remove ambient RNA by SoupX")
data.10x = list()
for (sample in samples) {
  data.10x[[sample]] = load10X(paste0(indir, sample, "/outs/"))
  data.10x[[sample]] <- autoEstCont(data.10x[[sample]], forceAccept=TRUE) #
  print(sample)
  data.10x[[sample]] <- adjustCounts(data.10x[[sample]])
}

message("Create Seurat object after SoupX")
scrna.list = list()
for (sample in samples) {
    scrna.list[[sample]] = CreateSeuratObject(counts = data.10x[[sample]], min.cells=3, project=sample)
}

message("Remove raw data to save memory")
rm(data.10x)

message("Add percent.mt and percent.rb to cell level metadata")
for (sample in samples) {
  scrna.list[[sample]][["percent.mt"]] <- PercentageFeatureSet(scrna.list[[sample]], pattern = mtPattern)
  scrna.list[[sample]][["percent.rb"]] <- PercentageFeatureSet(scrna.list[[sample]], pattern = rbPattern)
}

message("Rename nCount_RNA and nFeature_RNA")
for (sample in samples) {
  scrna.list[[sample]]$nUMI <- scrna.list[[sample]]$nCount_RNA
  scrna.list[[sample]]$nCount_RNA <- NULL
  scrna.list[[sample]]$nGene <- scrna.list[[sample]]$nFeature_RNA
  scrna.list[[sample]]$nFeature_RNA <- NULL
}

#message("Run doublet detection scripts")
#system2(command = "bash", args = c("run_scrublet_multi.sh"))

message("Read in doublet scores")
for (sample in samples){
  doublet_scores <- scan(paste0(scrubletdir, sample, "_srublet.score"))
  predicted_doublets <- scan(paste0(scrubletdir, sample, "_srublet.logic"))   
  ds <- as.data.frame(cbind(doublet_scores, predicted_doublets))
  ds$predicted_doublets <- as.logical(ds$predicted_doublets)
  rownames(ds) <- rownames(scrna.list[[sample]]@meta.data) 
  scrna.list[[sample]] <- AddMetaData(scrna.list[[sample]], ds)
  scrna.list[[sample]] <- subset(scrna.list[[sample]], subset=predicted_doublets == FALSE)
}


message("Feature plot before QC")
pdf(file = paste0(outdir, "individual/", "QC.plot.before.pdf"), width = 8, height = 8)
for (sample in samples){
  plot1 <- VlnPlot(scrna.list[[sample]], features = c("nGene", "nUMI", "percent.mt", "percent.rb"), ncol = 4, pt.size = 0.1) & theme(plot.title = element_text(size=10)) + theme(aspect.ratio=40/10)
  print(plot1 + coord_fixed())
  plot2 <- FeatureScatter(scrna.list[[sample]], feature1 = "nUMI", feature2 = "percent.mt") + theme(aspect.ratio=10/10)
  print(plot2 + coord_fixed())
  plot3 <- FeatureScatter(scrna.list[[sample]], feature1 = "nUMI", feature2 = "nGene") + theme(aspect.ratio=10/10)
  print(plot3 + coord_fixed())
}
dev.off()

message("Filtered cells with 3SD of mean nCount and nFeature, percent of mito")
qc_cutoff = 3
mito_cutoff = mitoCutoff
for (sample in samples){
  mean.nCount <- mean(scrna.list[[sample]]@meta.data$nUMI)
  sd.nCount <- sd(scrna.list[[sample]]@meta.data$nUMI)
  mean.nFeature <- mean(scrna.list[[sample]]@meta.data$nGene)
  sd.nFeature <- sd(scrna.list[[sample]]@meta.data$nGene)
  scrna.list[[sample]] <- subset(scrna.list[[sample]], subset = nUMI > mean.nCount - qc_cutoff*sd.nCount & nUMI < mean.nCount + qc_cutoff*sd.nCount & nGene > mean.nFeature - qc_cutoff*sd.nFeature & nGene < mean.nFeature + qc_cutoff*sd.nFeature & percent.mt < mito_cutoff)
}

message("Feature plot after QC")
pdf(file = paste0(outdir, "individual/", "QC.plot.after.pdf"))
for (sample in samples){
  plot1 <- VlnPlot(scrna.list[[sample]], features = c("nGene", "nUMI", "percent.mt", "percent.rb"), ncol = 4, pt.size = 0.1) & theme(plot.title = element_text(size=10)) + theme(aspect.ratio=40/10)
  print(plot1 + coord_fixed())
  plot2 <- FeatureScatter(scrna.list[[sample]], feature1 = "nUMI", feature2 = "percent.mt") + theme(aspect.ratio=10/10)
  print(plot2 + coord_fixed())
  plot3 <- FeatureScatter(scrna.list[[sample]], feature1 = "nUMI", feature2 = "nGene") + theme(aspect.ratio=10/10)
  print(plot3 + coord_fixed())
}
dev.off()


message("Step 4: Sample processing, Normalization, Find variable features, Data scaling")
for (sample in samples){
  scrna.list[[sample]] <- NormalizeData(scrna.list[[sample]], normalization.method = "LogNormalize", scale.factor = 10000)
  scrna.list[[sample]] <- FindVariableFeatures(scrna.list[[sample]], selection.method = "vst", nfeatures = 2000)
  all.genes <- rownames(scrna.list[[sample]])
  scrna.list[[sample]] <- ScaleData(scrna.list[[sample]], features = all.genes, verbose = FALSE)
}


message("Step 5: Determine the ‘dimensionality’ of the dataset")
pdf(file = paste0(outdir, "individual/", "elbow.plot.pdf"), width = 8, height = 8)
for (sample in samples){
  scrna.list[[sample]] <- RunPCA(scrna.list[[sample]], verbos = FALSE)
  plot1 <- ElbowPlot(scrna.list[[sample]]) + ggtitle(sample) + theme(aspect.ratio=5/10) + theme(plot.margin = unit(c(3, 3, 3, 3), "cm"))
  print(plot1 + coord_fixed())
}
dev.off()


message("Step 6: Cell clustering")
scrna.default.list = list()
for (sample in samples){
  scrna.default.list[[sample]] <- scrna.list[[sample]]
}

message("Cell clustering using default settings: PCA, Louvain. CHANGE dims according to elbow plot !!!")
if (flag==1 | flag==0){
  pdf(file = paste0(outdir, "individual/", "cluster.umap.louvain.pdf"), width = 6, height = 6)
for (sample in samples){
  scrna.list[[sample]] <- RunPCA(scrna.list[[sample]], npcs = 15, verbose = FALSE) 
  scrna.list[[sample]] <- FindNeighbors(scrna.list[[sample]], reduction = "pca", dims = 1:15) # picked 15
  scrna.list[[sample]] <- FindClusters(scrna.list[[sample]], resolution = 0.5)
  scrna.list[[sample]] <- RunUMAP(scrna.list[[sample]], reduction = "pca", dims = 1:15) 
  plot1 <- DimPlot(scrna.list[[sample]], reduction = "umap", label = TRUE, pt.size = 1) + ggtitle(sample)
  print(plot1 + coord_fixed())
}
dev.off()
}

message("Cell clustering using GLMPCA")
if (flag==2 | flag==0){
  pdf(file = paste0(outdir, "individual/", "cluster.umap.glmpca.pdf"), width = 6, height = 6)
scrna.glmpca.list = list()
for (sample in samples){
  scrna.glmpca.list[[sample]] <- RunGLMPCA(scrna.default.list[[sample]], L = 15)
  scrna.glmpca.list[[sample]] <- FindNeighbors(scrna.glmpca.list[[sample]], reduction = "glmpca", dims = 1:15)
  scrna.glmpca.list[[sample]] <- FindClusters(scrna.glmpca.list[[sample]], resolution = 0.5)
  scrna.glmpca.list[[sample]] <- RunUMAP(scrna.glmpca.list[[sample]], reduction = "glmpca", dims = 1:15)
  plot1 <- DimPlot(scrna.glmpca.list[[sample]], reduction = "umap", label = TRUE, pt.size = 1) + ggtitle(sample)
  print(plot1 + coord_fixed())
}
dev.off()
}

message("Cell clustering using Leiden")
if (flag==3 | flag==0){
pdf(file = paste0(outdir, "individual/", "cluster.umap.leiden.pdf"), width = 6, height = 6)
scrna.leiden.list = list()
for (sample in samples){
  scrna.leiden.list[[sample]] <- RunPCA(scrna.default.list[[sample]], npcs = 15, verbose = FALSE) 
  scrna.leiden.list[[sample]] <- FindNeighbors(scrna.leiden.list[[sample]], reduction = "pca", dims = 1:15) # picked 10
  scrna.leiden.list[[sample]] <- FindClusters(scrna.leiden.list[[sample]], resolution = 0.5, algorithm = 4)
  scrna.leiden.list[[sample]] <- RunUMAP(scrna.leiden.list[[sample]], reduction = "pca", dims = 1:15)
  plot1 <- DimPlot(scrna.leiden.list[[sample]], reduction = "umap", label = TRUE, pt.size = 1) + ggtitle(sample)
  print(plot1 + coord_fixed())
}
dev.off()  
}
