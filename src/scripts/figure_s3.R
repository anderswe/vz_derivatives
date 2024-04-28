# Info --------------------------------------------------------------------

# Anders E.
# Taylor lab
# April 27, 2024


# Notes -------------------------------------------------------------------

# Supplemental Figure 3A: Featureplots for interneuron markers



# Libraries ---------------------------------------------------------------
# conda activate ../G4MB/src/envs/r_general
# R
library(Seurat)
library(magrittr)
library(ggplot2)
library(patchwork)
library(yaml)
options(future.globals.maxSize = 1.2 * 1024 ^ 3)
yml <- yaml::read_yaml("config.yaml")
source("src/scripts/utils.R")


# Import data -------------------------------------------------------------


# import dev human data
ald <- readRDS(glue::glue("{yml$ald_path}/cbl_integrated_cleanCC_210111.rds"))
br <- qs::qread(glue::glue("{yml$lin_dir}/hb_cb_final_v3.qs"))
cs20 <- qs::qread("out/cs20s_annot.qs")

# import adult human data
lake <- qs::qread("src/data/lake_cerebellum.qs")


# Add celltype annotation for just PCs and PIPs ---------------------------------------------

ald@meta.data %<>% dplyr::mutate(celltype_aggr = dplyr::case_when(
  figure_clusters == "01-PC" ~ "PC",
  figure_clusters == "07-PIP" ~ "PIP",
  TRUE ~ "Other") %>% factor(levels = c("PC", "PIP", "Other")),
  age_num = stringr::str_replace(age, " PCW", "") %>% as.numeric())

br@meta.data %<>% dplyr::mutate(celltype_aggr = dplyr::case_when(
  seurat_clusters %in% c(1,0,3,15) ~ "PC",
  seurat_clusters %in% c(14) ~ "PIP",
  TRUE ~ "Other") %>% factor(levels = c("PC", "PIP", "Other")))

cs20@meta.data %<>% dplyr::mutate(celltype_pc_only = ifelse(grepl("PC", celltype_broad), "PC", "Other") %>% factor(levels = c("PC", "Other")))

lake@meta.data %<>% dplyr::mutate(celltype_aggr = dplyr::case_when(
  grepl("Ast", celltype) ~ "Astr.",
  celltype == "End" ~ "Endo.",
  celltype == "Gran" ~ "Gran.",
  celltype == "Mic" ~ "Micro.",
  celltype == "Oli" ~ "Oligo.",
  grepl("OPC", celltype) ~ "OPC",
  celltype == "Per" ~ "Peri.",
  grepl("Purk", celltype) ~ "PC"))
lake@meta.data %<>% dplyr::mutate(celltype_pc_only = ifelse(celltype_aggr == "PC", "PC", "Other") %>% factor(levels = c("PC", "Other")))



# S3A ----------------------------------------------------


# dimplots
d_ald <- DimPlot(ald, group.by = "celltype_aggr", cols = c(tableau20[1:2], "lightgrey"), pt.size = 0.5, raster = T) & NoAxes() & NoLegend() & ggtitle("Aldinger et al.")
d_br <- DimPlot(br, group.by = "celltype_aggr", cols = c(tableau20[1:2], "lightgrey"), pt.size = 0.5, raster = T) & NoAxes() & NoLegend() & ggtitle("Braun et al.")
d_cs20 <- DimPlot(cs20, group.by = "celltype_pc_only", cols = c(tableau20[1], "lightgrey"), pt.size = 0.5, raster = T) & NoAxes() & NoLegend() & ggtitle("This study")
d_lake <- DimPlot(lake, group.by = "celltype_pc_only", cols = c(tableau20[1], "lightgrey"), pt.size = 0.5, raster = T) & NoAxes() & NoLegend() & ggtitle("Lake et al.")

pc_dimplots <- d_cs20 | d_ald | d_br | d_lake

# featureplots
interneuron_genes <- c('KCNA2', 'GRIK3', 'NXPH1', 'SLC6A5', 'ALDH1A3', 'LGI2','GJD2', 'SORCS3') # cf. Yu 2021 Nat NSc

DefaultAssay(ald) <- "RNA"
f_ald <- FeaturePlot(ald, interneuron_genes, raster = T, order = T, ncol = 1, pt.size = 2.5) & NoAxes() & NoLegend()
f_br <- FeaturePlot(br, interneuron_genes, raster = T, order = T, ncol = 1, pt.size = 2.5) & NoAxes() & NoLegend()
f_lake <- FeaturePlot(lake, interneuron_genes, raster = T, order = T, ncol = 1, pt.size = 2.5) & NoAxes() & NoLegend()
f_cs20 <- FeaturePlot(cs20, interneuron_genes, raster = T, order = T, ncol = 1, pt.size = 2.5) & NoAxes() & NoLegend()


pdf("out/feats_interneurons.pdf", w = 16, h = 36)
(pc_dimplots) / (f_cs20 | f_ald | f_br | f_lake) + patchwork::plot_layout(heights = c(1, length(interneuron_genes)))
dev.off()