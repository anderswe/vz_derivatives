# Info --------------------------------------------------------------------
# Anders E.
# Oct 23, 2023
# Taylor lab

# Notes -------------------------------------------------------------------
# 
# idea is to project older purkinje cells onto CS21-23 object
# for a closer, semi-supervised look at the representation of PC subtypes at later timepoints

# RProfile ----------------------------------------------------------------
# conda activate r_general
# R

# Libraries ---------------------------------------------------------------
library(qs)
library(Seurat)
library(Signac)
library(magrittr)
library(tibble)
library(ggplot2)
library(harmony)
source("src/scripts/utils.R")
yml <- yaml::read_yaml("config.yaml")

# Inputs ------------------------------------------------------------------

# outdir
outdir <- glue::glue("out/spca")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# objects
ref <- qs::qread("out/cs20s_annot.qs")

# query <- readRDS(glue::glue("{yml$ald_path}/cbl_integrated_cleanCC_210111.rds")) %>% 
    # subset(subset = fig_cell_type == "H-PC") %>% 
    # SCTransform(method = 'glmGamPoi') %>% 
    # RunPCA() %>% 
    # FindNeighbors(dims = 1:20) %>%
    # FindClusters(resolution = 0.2) %>% 
    # RunUMAP(dims = 1:20)
# qs::qsave(query, glue::glue("{yml$g4mb_clone_path}/sandbox/aldinger_pc_subset.qs"))
query <- qs::qread(glue::glue("{yml$g4mb_clone_path}/sandbox/aldinger_pc_subset.qs"))


# sPCA ----------------------------------------------------------------------

# default assays
DefaultAssay(ref) <- "SCT"
DefaultAssay(query) <- "SCT"

# rerun umap to return model
ref %<>% RunUMAP(
    assay = "RNA",
    verbose = TRUE,
    return.model = TRUE,
    dims = 1:20)

# run sPCA
ref %<>% RunSPCA(assay = 'SCT', graph = 'SCT_snn')

# find sPCA neighbors
ref %<>% FindNeighbors(
    reduction = "spca",
    dims = 1:50,
    graph.name = "spca.annoy.neighbors", 
    k.param = 50,
    cache.index = TRUE,
    return.neighbor = TRUE,
    l2.norm = TRUE)

# find anchors
anchors <- FindTransferAnchors(
        reference = ref,
        query = query,
        # k.filter = NA,
        reference.reduction = "spca", 
        reference.neighbors = "spca.annoy.neighbors", 
        dims = 1:50)

# map
query  <- MapQuery(anchorset = anchors, 
        query = query,
        reference = ref, 
        refdata = list( 
          celltype_fine = "celltype_fine",
          celltype_broad = "celltype_broad",
          cytotrace = "CytoTRACE"),
        reference.reduction = "spca",
        reduction.model = "umap")

# pdf("sandbox/cytotrace.pdf", w = 4, h = 4)
# FeaturePlot(ref, "CytoTRACE", cols = viridis::viridis(100), raster = T)+NoAxes()+NoLegend()+ggtitle("")
# dev.off()

set.seed(667)
ct_pal <- sample(tableau20, length(unique(ref$celltype_fine))) %>% magrittr::set_names(levels(ref$celltype_fine) %>% rev()) # #79706E, #8CD17D, #D4A6C8, #499894, #FABFD2, #B07AA1, #86BCB6, #D37295, #A0CBE8, #D7B5A6

# plots - celltype
p1 <- DimPlot(ref, reduction = 'umap', group.by = 'celltype_fine', cols = ct_pal)+NoAxes()+ggtitle("CS21-23")
p2 <- DimPlot(query, reduction = 'ref.umap', group.by = 'predicted.celltype_fine', cols = ct_pal)+NoAxes()+NoLegend()+ggtitle("Aldinger PCs (Proj.)")
p3 <- DimPlot(query, reduction = 'umap', group.by = 'seurat_clusters', raster = T, cols = tableau20)+NoAxes()+NoLegend()+ggtitle("Aldinger PCs")

# plots - pseudotime
p4 <- FeaturePlot(ref, "CytoTRACE", reduction = 'umap', cols = viridis::viridis(100))+NoAxes()+ggtitle("CytoTRACE")
p5 <- FeaturePlot(query, "predicted.cytotrace", reduction = 'ref.umap', cols = viridis::viridis(100))+NoAxes()+NoLegend()#+viridis::scale_color_viridis(option = "viridis", limits = c(min(ref$CytoTRACE, na.rm = TRUE), max(ref$CytoTRACE, na.rm = TRUE)))
p6 <- FeaturePlot(query, "predicted.cytotrace", reduction = 'umap', cols = viridis::viridis(100))+NoAxes()+NoLegend()#+viridis::scale_color_viridis(option = "viridis", limits = c(min(ref$CytoTRACE, na.rm = TRUE), max(ref$CytoTRACE, na.rm = TRUE)))

# plots - phase
p7 <- DimPlot(ref, reduction = 'umap', group.by = 'Phase', cols = phase_pal)+NoAxes()
p8 <- DimPlot(query, reduction = 'ref.umap', group.by = 'Phase', cols = phase_pal)+NoAxes()+ggtitle("")
p9 <- DimPlot(query, reduction = 'umap', group.by = 'Phase', raster = T, cols = phase_pal)+NoAxes()+ggtitle("")

# print plots
pdf(glue::glue("{outdir}/cs_aldinger_mapped.pdf"), w = 15, h = 15)
p1 + p2 + p3 + p4 + p5 + p6 + p7 + p8 + p9 + patchwork::plot_layout(ncol = 3, guides = "collect")
dev.off()

# query md
qmd <- query@meta.data %>% 
  dplyr::select(seurat_clusters, predicted.celltype_fine.score:predicted.cytotrace) %>% 
  tibble::rownames_to_column("barcode") %>% 
  dplyr::left_join(Embeddings(query, "ref.umap") %>% tibble::as_tibble(rownames = "barcode"), by = "barcode")

# cleaner plot
g <- ggplot()+
    geom_point(data = as.data.frame(Embeddings(ref, "umap")), aes(UMAP_1, UMAP_2), color = "lightgrey", size = 0.5)+
    geom_point(data = qmd, aes(refUMAP_1, refUMAP_2, color = predicted.celltype_fine), size = 0.5)+
    scale_color_manual(values = ct_pal, name = "")+
    theme_classic()+
    labs(x = "", y = "")+
    guides(colour = guide_legend(override.aes = list(size=3)))+
    theme(axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank())

    
pdf(glue::glue("{outdir}/proj_clean_cs_ald.pdf"), h = 5, w = 5)
print(g)
dev.off()

# final plots ---

# dimplots on original Purkinje embedding
d <- DimPlot(query, group.by = "predicted.celltype_fine", cols = ct_pal, reduction = "umap") + NoAxes() + ggplot2::ggtitle("")
d2 <- DimPlot(query, group.by = "seurat_clusters", cols = tableau10, reduction = "umap") + NoAxes() + NoLegend() + ggplot2::ggtitle("")

# dimplots on ref embedding
e <- ggplot()+
    geom_point(data = as.data.frame(Embeddings(ref, "umap")), aes(UMAP_1, UMAP_2), color = "lightgrey", size = 0.5)+
    geom_point(data = qmd, aes(refUMAP_1, refUMAP_2, color = predicted.celltype_fine), size = 0.5)+
    scale_color_manual(values = ct_pal, name = "")+
    theme_classic()+
    labs(x = "", y = "")+
    guides(colour = guide_legend(override.aes = list(size=3)))+
    theme(axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank())
e2 <- ggplot()+
    geom_point(data = as.data.frame(Embeddings(ref, "umap")), aes(UMAP_1, UMAP_2), color = "lightgrey", size = 0.5)+
    geom_point(data = qmd, aes(refUMAP_1, refUMAP_2, color = seurat_clusters), size = 0.5)+
    scale_color_manual(values = tableau10, name = "")+
    theme_classic()+
    labs(x = "", y = "")+
    guides(colour = guide_legend(override.aes = list(size=3)))+
    theme(axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank())

# plot
pdf(glue::glue("{outdir}/proj_assignments.pdf"), h = 10, w = 12)
print(d + d2 + e + e2 + patchwork::plot_layout(guides = "collect", ncol = 2))
dev.off()


# Save ---------------------------------------------------------------------

# query metadata
qmd %>% data.table::fwrite(glue::glue("{outdir}/spca_metadata.csv"))
