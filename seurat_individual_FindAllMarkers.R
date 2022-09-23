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

message("Read in gene name table")
geneTable <- read.csv(paste0(refdir, "geneAnnotationTable.csv"), header = T, row.names = 1)

message("Import seurat object")
scrna.list <- readRDS(paste0(outdir, "individual/", "scrna.list.seurat.", projectName, ".rds"))

message("Finding differentially expressed markers")
scrna.markers.list <- list()
markers.topN.list <- list()
for (sample in samples){
  scrna.markers <- FindAllMarkers(scrna.list[[sample]], only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
  names(scrna.markers)[names(scrna.markers) == "gene"] <- "geneSymbol"
  scrna.markers <- cbind(scrna.markers, geneID=geneTable$geneID[match(scrna.markers$geneSymbol, geneTable$geneSymbol)])
  scrna.markers.list[[sample]] <- scrna.markers
  topN <- scrna.markers %>% group_by(cluster) %>% top_n(n = geneN, wt = avg_log2FC)
  markers.topN.list[[sample]] <- topN
  rm(scrna.markers)
  rm(topN)
}
openxlsx::write.xlsx(scrna.markers.list, paste0(outdir, "individual/", "markers.xls"))
openxlsx::write.xlsx(markers.topN.list, paste0(outdir, "individual/", "markers.topN.xls"))

message("Save the FindAllMarkers data")
saveRDS(scrna.markers.list, paste0(outdir, "individual/", "scrna.markers.list.seurat.", projectName, ".rds"))
saveRDS(markers.topN.list, paste0(outdir, "individual/", "markers.topN.list.seurat.", projectName, ".rds"))
