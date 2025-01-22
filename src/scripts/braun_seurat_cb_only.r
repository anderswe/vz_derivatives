# braun convert h5ad to seurat---
# module load R/4.2.1
# R
.libPaths(c("/hpf/largeprojects/mdtaylor/aerickson/data/clones/G4MB/src/envs/symphony/renv/library/R-4.2/x86_64-pc-linux-gnu", .libPaths())) 

library(SeuratDisk)
library(magrittr)
library(Seurat)
library(ggplot2)
source("src/scripts/utils.R")
options(future.globals.maxSize = 50 * 1024 ^ 3) # for 50G globals

# h5ad directory
h5_dir <- "/hpf/largeprojects/mdtaylor/aerickson/data/scRNAseq/linnarsson_loom/first_trimester/subset_cerebellum.h5ad"


# # convert
Convert(h5_dir, dest = "h5seurat", overwrite = TRUE)
so <- LoadH5Seurat((h5_dir %>% stringr::str_replace(".h5ad", ".h5seurat")), meta.data = FALSE, misc = FALSE) %>% Seurat::DietSeurat()
so <- so@assays$RNA@counts %>% Seurat::as.sparse() %>% SeuratObject::CreateSeuratObject()

# rescue metadata
md <- data.table::fread("/hpf/largeprojects/mdtaylor/aerickson/data/scRNAseq/linnarsson_loom/first_trimester/metadata_subset_cerebellum.csv") %>% 
  tibble::column_to_rownames("V1")
so %<>% Seurat::AddMetaData(md)
so %>% qs::qsave("/hpf/largeprojects/mdtaylor/aerickson/data/scRNAseq/linnarsson_loom/first_trimester/cb.qs")

if(file.exists("/hpf/largeprojects/mdtaylor/aerickson/data/scRNAseq/linnarsson_loom/first_trimester/cb.qs")){
  file.remove("/hpf/largeprojects/mdtaylor/aerickson/data/scRNAseq/linnarsson_loom/first_trimester/subset_cerebellum.h5seurat")
}
# NOT USED:
# # ID HOX-neg cells
# hox_genes <- grep("^HOX", rownames(so@assays$RNA), value = T)
# hox_neg_cells <- so@assays$RNA@counts[hox_genes,] %>%
#   as.matrix() %>% 
#   t() %>% 
#   Matrix::rowSums() %>% 
#   .[. == 0] %>% 
#   names()
# ss <- subset(so, cells = hox_neg_cells) %>% 

# # preprocessing

ss <- so %>% 
  # subset(CellClass %in% c("Glioblast", "Neural crest", "Neuroblast", "Neuron", "Neuronal IPC", "Oligo", "Radial glia")) %>% 
  NormalizeData(normalization.method = "LogNormalize", scale.factor = 10000) %>% 
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>% 
  ScaleData() %>% 
  RunPCA() %>% 
  harmony::RunHarmony(group.by.vars = "Donor", # Chemistry
             reduction = "pca",
             assay.use = "RNA",
             theta = 2) %>%
  FindNeighbors(dims = 1:30, reduction = "harmony") %>%
  FindClusters(resolution = 0.2) %>% 
  RunUMAP(dims = 1:30, reduction = "harmony")

# PICKUP!!!!
# 17:19:43 Using Annoy for neighbor search, n_neighbors = 30
# 17:19:43 Building Annoy index with metric = cosine, n_trees = 50
# 0%   10   20   30   40   50   60   70   80   90   100%
# [----|----|----|----|----|----|----|----|----|----|
# **************************************************|
# 17:20:56 Writing NN index file to temp file /localhd/8052315/RtmpdnCm1r/file25a937217a9b2c
# 17:20:56 Searching Annoy index using 1 thread, search_k = 3000
# 17:27:09 Annoy recall = 100%
# 17:27:09 Commencing smooth kNN distance calibration using 1 thread with target n_neighbors = 30
# 17:27:15 Initializing from normalized Laplacian + noise (using irlba)
# Error in irlba::irlba(L, nv = n, nu = 0, maxit = iters) : 
#   function 'as_cholmod_sparse' not provided by package 'Matrix'
# In addition: Warning messages:
# 1: Quick-TRANSfer stage steps exceeded maximum (= 7118100) 
# 2: Quick-TRANSfer stage steps exceeded maximum (= 7118100) 
# 3: Quick-TRANSfer stage steps exceeded maximum (= 7118100) 
# 4: Quick-TRANSfer stage steps exceeded maximum (= 7118100) 
# > ss %>% qs::qsave("/hpf/largeprojects/mdtaylor/aerickson/data/scRNAseq/linnarsson_loom/first_trimester/cb_clean.qs")
# Error in qs::qsave(., "/hpf/largeprojects/mdtaylor/aerickson/data/scRNAseq/linnarsson_loom/first_trimester/cb_clean.qs") : 
#   object 'ss' not found

# save
ss %>% qs::qsave("/hpf/largeprojects/mdtaylor/aerickson/data/scRNAseq/linnarsson_loom/first_trimester/cb_clean.qs")

# # plots
# ss %<>% CellCycleScoring(s.features = Seurat::cc.genes.updated.2019$s.genes, g2m.features = Seurat::cc.genes.updated.2019$g2m.genes, set.ident = FALSE)
# d1 <- DimPlot(ss, group.by = "Donor", raster = T, cols = sample(colorRampPalette(RColorBrewer::brewer.pal(9, "Spectral"))(length(unique(ss$Donor))))) + NoAxes() + NoLegend() + ggplot2::ggtitle("Donor")
# d2 <- DimPlot(ss, group.by = "Age", raster = T, cols = colorRampPalette(RColorBrewer::brewer.pal(9, "Spectral"))(length(unique(ss$Age)))) + NoAxes() + ggplot2::ggtitle("Age")
# d3 <- DimPlot(ss, group.by = "CellClass", raster = T, cols = tableau20) + NoAxes() + ggplot2::ggtitle("CellClass")
# d4 <- DimPlot(ss, group.by = "Chemistry", raster = T, cols = tableau10) + NoAxes() + ggplot2::ggtitle("Chemistry")
# d5 <- DimPlot(ss, group.by = "Tissue", raster = T, cols = tableau10) + NoAxes() + ggplot2::ggtitle("Tissue")
# d6 <- DimPlot(ss, group.by = "Subregion", raster = T, cols = tableau10) + NoAxes() + ggplot2::ggtitle("Subregion")
# d7 <- DimPlot(ss, group.by = "Subdivision", raster = T, cols = tableau10) + NoAxes() + ggplot2::ggtitle("Subdivision")
# d8 <- DimPlot(ss, group.by = "Phase", raster = T, cols = phase_pal) + NoAxes() + ggplot2::ggtitle("Phase")
# d9 <- DimPlot(ss, group.by = "seurat_clusters", raster = T, cols = tableau20) + NoAxes() + ggplot2::ggtitle("Cluster")

# pdf("tmp/dimplots_cb_hb_clean.pdf", w = 15, h = 15)
# cowplot::plot_grid(d1, d2, d3, d4, d5, d6, d7, d8, d9, ncol = 3)
# dev.off()

# fp(ss, "cb_hb_clean", dir = "tmp", reduction = "umap")


# reload----
lin_dir <- "/hpf/largeprojects/mdtaylor/aerickson/data/scRNAseq/linnarsson_loom/first_trimester"
# so <- qs::qread(glue::glue("{lin_dir}/cb_hb_clean.qs"))

# cc diff
# so$cc_difference <- so$S.Score - so$G2M.Score

# reprocess
# future::plan("multisession", workers = 8)

# ss <- so %>% subset(CellClass %in% c("Glioblast", "Neural crest", "Neuroblast", "Neuron", "Neuronal IPC", "Oligo", "Radial glia")) %>% 
#   NormalizeData(normalization.method = "LogNormalize", scale.factor = 10000) %>% 
#   FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>% 
#   ScaleData(vars.to.regress = "cc_difference") %>% 
#   RunPCA() %>% 
#   harmony::RunHarmony(group.by.vars = "Donor",
#              reduction = "pca",
#              assay.use = "RNA",
#              theta = 2) %>%
#   FindNeighbors(dims = 1:30, reduction = "harmony") %>%
#   FindClusters(resolution = 0.2) %>% 
#   RunUMAP(dims = 1:30, reduction = "harmony", min.dist = 0.2, spread = 1.2)
# qs::qsave(ss, glue::glue("{lin_dir}/hb_cb_final.qs"))


# recluster again---
# options(future.globals.maxSize = 10 * 1024^3)
# future::plan("multisession", workers = 2)

# import
# ss <- qs::qread(glue::glue("{lin_dir}/hb_cb_final.qs"))
# message("Object read.")

# # qc threshes
# ss$percent_mt <- PercentageFeatureSet(ss, pattern = "^MT-")
# mito_thresh <- 5 #median(ss$percent_mt) + 5*mad(ss$percent_mt)
# min_gn_thresh <- 500
# max_gn_thresh <- median(ss$nFeature_RNA) + 5*mad(ss$nFeature_RNA)
  
# # reprocess again
# sss <- ss %>% 
#   subset(subset = nFeature_RNA > min_gn_thresh & nFeature_RNA < max_gn_thresh & percent_mt < mito_thresh) %>% 
#   NormalizeData(normalization.method = "LogNormalize", scale.factor = 10000) %>% 
#   FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>% 
#   ScaleData(vars.to.regress = "cc_difference") %>% 
#   RunPCA() %>% 
#   harmony::RunHarmony(group.by.vars = "Donor",
#              reduction = "pca",
#              assay.use = "RNA",
#              theta = 2) %>%
#   FindNeighbors(dims = 1:30, reduction = "harmony") %>%
#   FindClusters(resolution = 0.4) %>% 
#   RunUMAP(dims = 1:30, reduction = "harmony", min.dist = 0.2, spread = 1.2)
# qs::qsave(sss, glue::glue("{lin_dir}/hb_cb_final_v2.qs"))
# message("Object reprocessed and saved.")

# # dimplots
# d1 <- DimPlot(sss, group.by = "Donor", raster = T, cols = sample(colorRampPalette(RColorBrewer::brewer.pal(9, "Spectral"))(length(unique(ss$Donor))))) + NoAxes() + NoLegend() + ggplot2::ggtitle("Donor")
# d2 <- DimPlot(sss, group.by = "Age", raster = T, cols = colorRampPalette(RColorBrewer::brewer.pal(9, "Spectral"))(length(unique(ss$Age)))) + NoAxes() + ggplot2::ggtitle("Age")
# d3 <- DimPlot(sss, group.by = "CellClass", raster = T, cols = tableau10) + NoAxes() + ggplot2::ggtitle("CellClass")
# d4 <- DimPlot(sss, group.by = "Chemistry", raster = T, cols = tableau10) + NoAxes() + ggplot2::ggtitle("Chemistry")
# d5 <- DimPlot(sss, group.by = "Tissue", raster = T, cols = tableau10) + NoAxes() + ggplot2::ggtitle("Tissue")
# d6 <- DimPlot(sss, group.by = "Phase", raster = T, cols = phase_pal) + NoAxes() + ggplot2::ggtitle("Phase")
d7 <- DimPlot(sss, group.by = "seurat_clusters", raster = T, cols = tableau20, label = T, repel = T) + NoAxes() + ggplot2::ggtitle("Cluster")

# pdf(glue::glue("{lin_dir}/dimplots_cb_hb_final_v2.pdf"), w = 15, h = 15)
# cowplot::plot_grid(d1, d2, d3, d4, d5, d6, d7, ncol = 3)
# dev.off()

# # featureplots
# fp_genes <- genes %>% .[. %in% rownames(sss)]
# pdf(glue::glue("{lin_dir}/featureplots_cb_hb_final_v2.pdf"), w = 15, h = (length(fp_genes)/2))
# FeaturePlot(sss, fp_genes, order = T, raster = T, ncol = 6) & NoAxes() & NoLegend()
# dev.off()

# # qc
# f1 <- FeaturePlot(sss, c("nFeature_RNA", "nCount_RNA", "percent_mt"), order = F, raster = T, ncol = 3) & NoAxes()
# f2 <- VlnPlot(sss, c("nFeature_RNA", "nCount_RNA", "percent_mt"), ncol = 3, pt.size = 0, cols = tableau20)

# pdf(glue::glue("{lin_dir}/qc_cb_hb_final_v2.pdf"), w = 18, h = 10)
# f2 / (d7 + patchwork::plot_spacer() + patchwork::plot_spacer()) + patchwork::plot_layout(guides = "collect")
# dev.off()

# markers
# mks_11 <- FindMarkers(sss, ident.1 = "11", only.pos = T) %>% dplyr::mutate(unique = pct.1-pct.2)
# genes_11 <- mks_11 %>% dplyr::arrange(p_val_adj, desc(avg_log2FC)) %>% head(20) %>% rownames()
# pdf(glue::glue("{lin_dir}/cl11_mks.pdf"), w = 20, h = 24)
# FeaturePlot(sss, genes_11, raster = T, ncol = 5) & NoAxes()
# dev.off()

# extra qc
extra_qc <- data.table::fread(glue::glue("{lin_dir}/extra_qc.csv")) %>% 
  dplyr::select(-V1) %>% 
  dplyr::mutate(cell_id = stringr::str_replace(cell_id, "b\\'", "") %>% stringr::str_replace("\\'", "")) %>% 
  dplyr::filter(cell_id %in% rownames(sss@meta.data)) %>%
  tibble::column_to_rownames("cell_id")
sss %<>% AddMetaData(extra_qc)
f3 <- VlnPlot(sss, c("doublet_score", "unspliced_fraction"), ncol = 2, pt.size = 0, cols = tableau20)
d8 <- DimPlot(sss, group.by = "doublet_flag", raster = T, cols = tableau10) + NoAxes() + ggplot2::ggtitle("Doublet flag")
pdf(glue::glue("{lin_dir}/extra_qc_cb_hb_final_v2.pdf"), w = 10, h = 10)
f3 / (d7 + d8) + patchwork::plot_layout(guides = "collect")
dev.off()



# reprocess...again again.... ---
sss <- qs::qread(glue::glue("{lin_dir}/hb_cb_final_v2.qs"))

# parallelize
# options(future.globals.maxSize = 10 * 1024^3)
# future::plan("multisession", workers = 2)


# reprocess
ssss <- sss %>% 
  subset(subset = seurat_clusters != "11") %>% # expresses LMX1A and SKOR2...presumed doublets
  NormalizeData(normalization.method = "LogNormalize", scale.factor = 10000) %>% 
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>% 
  ScaleData(vars.to.regress = "cc_difference") %>% 
  RunPCA() %>% 
  harmony::RunHarmony(group.by.vars = "Donor",
             reduction = "pca",
             assay.use = "RNA",
             theta = 2) %>%
  FindNeighbors(dims = 1:30, reduction = "harmony") %>%
  FindClusters(resolution = 0.4) %>% 
  RunUMAP(dims = 1:30, reduction = "harmony")
qs::qsave(ssss, glue::glue("{lin_dir}/hb_cb_final_v3.qs"))
# ssss <- qs::qread(glue::glue("{lin_dir}/hb_cb_final_v3.qs"))

d1 <- DimPlot(ssss, group.by = "Donor", raster = T, cols = sample(colorRampPalette(RColorBrewer::brewer.pal(9, "Spectral"))(length(unique(ssss$Donor))))) + NoAxes() + NoLegend() + ggplot2::ggtitle("Donor")
d2 <- DimPlot(ssss, group.by = "Age", raster = T, cols = colorRampPalette(RColorBrewer::brewer.pal(9, "Spectral"))(length(unique(ssss$Age)))) + NoAxes() + ggplot2::ggtitle("Age")
d3 <- DimPlot(ssss, group.by = "CellClass", raster = T, cols = tableau10) + NoAxes() + ggplot2::ggtitle("CellClass")
d4 <- DimPlot(ssss, group.by = "Chemistry", raster = T, cols = tableau10) + NoAxes() + ggplot2::ggtitle("Chemistry")
d5 <- DimPlot(ssss, group.by = "Tissue", raster = T, cols = tableau10) + NoAxes() + ggplot2::ggtitle("Tissue")
d6 <- DimPlot(ssss, group.by = "Phase", raster = T, cols = phase_pal) + NoAxes() + ggplot2::ggtitle("Phase")
d7 <- DimPlot(ssss, group.by = "seurat_clusters", raster = T, cols = colorRampPalette(tableau20)(length(unique(ssss$seurat_clusters))), label = T, repel = T) + NoAxes() + ggplot2::ggtitle("Cluster")

pdf(glue::glue("{lin_dir}/dimplots_cb_hb_final_v3.pdf"), w = 15, h = 15)
cowplot::plot_grid(d1, d2, d3, d4, d5, d6, d7, ncol = 3)
dev.off()

# featureplots
fp_genes <- genes %>% .[. %in% rownames(ssss)]
pdf(glue::glue("{lin_dir}/featureplots_cb_hb_final_v3.pdf"), w = 15, h = (length(fp_genes)/2))
FeaturePlot(ssss, fp_genes, order = T, raster = T, ncol = 6) & NoAxes() & NoLegend()
dev.off()

# qc
f1 <- FeaturePlot(ssss, c("nFeature_RNA", "nCount_RNA", "percent_mt"), order = F, raster = T, ncol = 3) & NoAxes()
f2 <- VlnPlot(ssss, c("nFeature_RNA", "nCount_RNA", "percent_mt"), ncol = 3, pt.size = 0, cols = colorRampPalette(tableau20)(length(unique(ssss$seurat_clusters))))

pdf(glue::glue("{lin_dir}/qc_cb_hb_final_v3.pdf"), w = 18, h = 10)
f2 / (d7 + patchwork::plot_spacer() + patchwork::plot_spacer()) + patchwork::plot_layout(guides = "collect")
dev.off()



# pdf(glue::glue("{lin_dir}/featureplots_genes_12.pdf"), w = 20, h = 16)
# FeaturePlot(ss, genes_12, order = T, raster = T, ncol = 6) & NoAxes() & NoLegend()
# dev.off()






# rerun CB only + microenv -----------------------------------------------------
# 2023-10-26
# module load R/4.2.1
# R

.libPaths(c("/hpf/largeprojects/mdtaylor/aerickson/data/clones/G4MB/src/envs/symphony/renv/library/R-4.2/x86_64-pc-linux-gnu", .libPaths())) 

library(SeuratDisk)
library(magrittr)
library(Seurat)
library(ggplot2)
source("src/scripts/utils.R")

# h5ad directory
# h5_dir <- "/hpf/largeprojects/mdtaylor/aerickson/data/scRNAseq/linnarsson_loom/first_trimester/subset_cerebellum.h5ad"

# convert
# Convert(h5_dir, dest = "h5seurat", overwrite = TRUE)
# so <- LoadH5Seurat((h5_dir %>% stringr::str_replace(".h5ad", ".h5seurat")), meta.data = FALSE, misc = FALSE) #%>% Seurat::DietSeurat()
# so <- so@assays$RNA@counts %>% Seurat::as.sparse() %>% SeuratObject::CreateSeuratObject()

# rescue metadata
# md <- data.table::fread("/hpf/largeprojects/mdtaylor/aerickson/data/scRNAseq/linnarsson_loom/first_trimester/metadata_subset_cerebellum.csv") %>% 
#   tibble::column_to_rownames("V1")
# so %<>% Seurat::AddMetaData(md)
# so %>% qs::qsave("out/braun/braun_cb_microenv.qs")
so <- qs::qread("out/braun/braun_cb_microenv.qs")

# qc
dir.create("out/braun/qc", showWarnings = FALSE, recursive = TRUE)
so[["percent_mt"]] <- PercentageFeatureSet(so, pattern = "^MT-")

pdf("out/braun/qc/qc_initial.pdf", w = 12, h = 4)
VlnPlot(so, features = c("nFeature_RNA", "nCount_RNA", "percent_mt"), ncol = 3, pt.size = 0) & NoLegend() & theme(axis.title.x = element_blank())
dev.off()

# preprocess
min_gn_thresh <- 500
max_gn_thresh <- median(so$nFeature_RNA) + 5*mad(so$nFeature_RNA)
so %<>% 
  subset(subset = nFeature_RNA > min_gn_thresh & nFeature_RNA < max_gn_thresh & percent_mt < 5) %>% 
  NormalizeData(normalization.method = "LogNormalize", scale.factor = 10000) %>% 
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>% 
  ScaleData() %>% 
  RunPCA() %>% 
  harmony::RunHarmony(group.by.vars = "Donor", # Chemistry
              reduction = "pca",
              assay.use = "RNA",
              theta = 2) %>%
  FindNeighbors(dims = 1:20, reduction = "harmony") %>%
  FindClusters(resolution = 0.2) %>% 
  RunUMAP(dims = 1:20, reduction = "harmony")

pdf("out/braun/qc/dim_prelim.pdf", w = 6, h = 5)
DimPlot(so, group.by = "CellClass", raster = T, cols = c(custom_colors, tableau20)) & NoAxes() & ggplot2::ggtitle("")
dev.off()

pdf("out/braun/qc/feats.pdf", w = 25, h = 90)
FeaturePlot(so, genes, raster = T, ncol = 5) & NoAxes() & NoLegend()
dev.off()

so %>% qs::qsave("out/braun/braun_cb_microenv.qs")
so <- qs::qread("out/braun/braun_cb_microenv.qs")