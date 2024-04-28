# Info --------------------------------------------------------------------

# Anders E.
# Taylor lab
# April 21, 2023


# Notes -------------------------------------------------------------------

# 
# original plots by Liam H.
# 


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

# Params ------------------------------------------------------------------

dims <- 20
res <- 0.2
pip_colour <- "#62b10a"
mypal <- RColorBrewer::brewer.pal(9, "RdBu")[c(1,5,9)]

# Aldinger ------------------------------------------------------------------

# from Aldinger et al. 2021 Nat NSc.
# so <- readRDS("src/data/cbl_integrated_cleanCC_210111.rds")
# on hpf:
so <- readRDS(glue::glue("{yml$ald_path}/cbl_integrated_cleanCC_210111.rds")) %>% Seurat::UpdateSeuratObject()

# update metadata
so@meta.data %<>% dplyr::mutate(celltype_aggr = dplyr::case_when(
  figure_clusters == "01-PC" ~ "PC",
  figure_clusters == "07-PIP" ~ "PIP",
  TRUE ~ "Other") %>% factor(levels = c("PC", "PIP", "Other")),
  age_num = stringr::str_replace(age, " PCW", "") %>% as.numeric()
)


# featureplots
# DefaultAssay(so) <- "SCT"
# ald_genes <- genes %>% .[. %in% rownames(so)]
# pdf("out/fp_ald.pdf", w = 25, h = length(ald_genes) / 1.5)
# FeaturePlot(so, ald_genes, raster = T, order = F, ncol = 5) & NoAxes() & NoLegend()
# dev.off()

# pip markers
# future::plan("multisession", workers = 6)
# Idents(so) <- so$figure_clusters
# pip_mks <- FindMarkers(so, ident.1 = "07-PIP", only.pos = T) %>% dplyr::mutate(unique = pct.1 - pct.2)
# qs::qsave(pip_mks, "src/out/pip_markers.qs")
# pip_mks <- qs::qread("src/out/pip_markers.qs")

# pip heatmap
# pip_keep <- pip_mks %>% dplyr::slice_min(p_val_adj, n = 60) %>% rownames()
# pip_mks %>% dplyr::arrange((pct.2))

# subset to pips
# ss_cells <- so@meta.data %>% dplyr::filter(figure_clusters == "07-PIP") %>% rownames()

# pdf("src/out/pip_highlight_dimplot.pdf", h = 5, w = 5)
# DimPlot(so, cells.highlight = ss_cells, cols.highlight = pip_colour, sizes.highlight = 0.1, na.value = "lightgrey", pt.size = 0.1, cols = "lightgrey")+theme_void()+NoLegend()
# dev.off()

# ss <- so %>% 
#   subset(cells = ss_cells) %>% 
#   SCTransform(assay = "RNA") %>%
#   RunPCA() %>%
#   FindNeighbors(dims = 1:dims) %>% 
#   FindClusters(resolution = res) %>%
#   RunUMAP(dims = 1:dims)
# qs::qsave(ss, "src/data/ald_pip_subset.qs")
# ss <- qs::qread("src/data/ald_pip_subset.qs")


# pip plots
# DefaultAssay(ss) <- "RNA"
# dp <- DimPlot(ss, group.by = "age", cols = RColorBrewer::brewer.pal(10, "Spectral"))+NoAxes()+labs(title="", color = "Age")
# fp <- FeaturePlot(ss, c("PVALB", "SST", "PAX2", "SKOR2"), order = T, cols = c("lightgrey", "darkblue"), ncol = 4) & NoAxes()
# dp + fp + patchwork::plot_layout(ncol = 2, widths = c(1,4))

# subset to PCs
# pc_cells <- so@meta.data %>% dplyr::filter(figure_clusters == "01-PC") %>% rownames()

# pc <- so %>% 
#   subset(cells = pc_cells) %>% 
#   NormalizeData(assay = "RNA") %>% 
#   ScaleData() %>% 
#   FindVariableFeatures() %>% 
#   RunPCA() %>%
#   FindNeighbors(dims = 1:dims) %>% 
#   FindClusters(resolution = res) %>%
#   RunUMAP(dims = 1:dims)
# qs::qsave(pc, "src/data/ald_pc_subset.qs")
# pc <- qs::qread("src/data/ald_pc_subset.qs")

# quick plots
# DefaultAssay(pc) <- "RNA"
# dp <- DimPlot(pc, group.by = "age", cols = RColorBrewer::brewer.pal(10, "Spectral"))+NoAxes()+labs(title="", color = "Age")
# fp <- FeaturePlot(pc, c("PVALB", "SST", "PAX2", "SKOR2", "PCP4"), order = T, cols = c("lightgrey", "darkblue"), ncol = 5) & NoAxes()
# dp + fp + patchwork::plot_layout(ncol = 2, widths = c(1,5))


# Linnarsson --------------------------------------------------------------

# from Braun 2023 Science
lin <- qs::qread(glue::glue("{yml$lin_dir}/hb_cb_final_v3.qs"))
lin@meta.data %<>% dplyr::mutate(celltype_aggr = dplyr::case_when(
  seurat_clusters %in% c(1,0,3,15) ~ "PC",
  seurat_clusters %in% c(14) ~ "PIP",
  TRUE ~ "Other") %>% factor(levels = c("PC", "PIP", "Other")))

# featurescatter
# pdf("out/feature_scatter_lin.pdf", h = 5, w = 10)
# FeatureScatter(lin, "PAX2", "SKOR2", group.by = "seurat_clusters", pt.size = 0.1)+NoLegend()+geom_jitter()+ggtitle("")
# dev.off()

# lin pip markers
# Idents(lin) <- lin$celltype
# lin_pip_mks <- FindMarkers(lin, ident.1 = "PIP (PAX2+)", only.pos = T) %>% dplyr::mutate(unique = pct.1 - pct.2, gene = rownames(.))
# qs::qsave(lin_pip_mks, "out/markers_linnarsson_pips.qs")

# lin pc markers
# lin_pc_mks <- FindMarkers(lin, ident.1 = "Purkinje (SKOR2+ PCP4+)", only.pos = T) %>% dplyr::mutate(unique = pct.1 - pct.2, gene = rownames(.))
# qs::qsave(lin_pc_mks, "out/markers_linnarsson_pcs.qs")

# Lake -------------------------------------------------------

# from Lake et al. 2017 Nat Biotech. (GSE97930)
# lake <- data.table::fread("src/data/GSE97930_CerebellarHem_snDrop-seq_UMI_Count_Matrix_08-01-2017.txt") %>% 
#   tibble::column_to_rownames("V1") %>% 
#   Seurat::as.sparse() %>% 
#   CreateSeuratObject() %>% 
#   SCTransform() %>%
#   RunPCA() %>%
#   FindNeighbors(dims = 1:dims) %>%
#   FindClusters(resolution = res) %>%
#   RunUMAP(dims = 1:dims)

# metadata from supp table 2
# lake_md <- readxl::read_xlsx("src/data/41587_2018_BFnbt4038_MOESM11_ESM.xlsx", sheet = 2, skip = 5) %>%
#   .[, 1:6] %>%
#   janitor::clean_names() %>%
#   dplyr::rename(sample_umi = 1, library = 2, experiment = 3, patient = 4, brain_region = 5, celltype = 6) %>%
#   dplyr::mutate(barcode = glue::glue("{celltype}_{sample_umi}")) %>%
#   dplyr::filter(barcode %in% colnames(lake)) %>%
#   tidyr::separate(sample_umi, into = c("sample", "umi"), remove = FALSE) %>% 
#   tibble::column_to_rownames("barcode")

# clin_md <- readxl::read_xlsx("src/data/41587_2018_BFnbt4038_MOESM11_ESM.xlsx", sheet = 1, skip = 3) %>% 
#   .[-1] %>%
#   janitor::clean_names() %>% 
#   dplyr::select(patient = patient_umb_number, age, sex) %>% 
#   dplyr::filter(patient %in% unique(lake_md$patient)) %>% 
#   dplyr::distinct() %>% 
#   dplyr::mutate(age = age %>% stringr::str_replace(" years", "") %>% as.numeric())
  

# add metadata
# lake %<>% AddMetaData(lake_md)
# lake@meta.data %<>% dplyr::mutate(age = plyr::mapvalues(patient, clin_md$patient, clin_md$age),
#                                   sex = plyr::mapvalues(patient, clin_md$patient, clin_md$sex))

# featureplots
# lake_genes <- genes %>% .[. %in% rownames(lake)]
# pdf("out/fp_lake.pdf", w = 25, h = length(fp_genes) / 1.5)
# FeaturePlot(lake, lake_genes, raster = T, order = T, ncol = 5) & NoAxes() & NoLegend()
# dev.off()

# save
# qs::qsave(lake, "src/data/lake_cerebellum.qs")
lake <- qs::qread("src/data/lake_cerebellum.qs")

# clean celltypes
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

# dimplots
# sample_plot <- DimPlot(lake, group.by = "sample", cols = RColorBrewer::brewer.pal(9, "Spectral"))+NoAxes()+labs(color = "Library", title = "")
# pt_plot <- DimPlot(lake, group.by = "patient", cols = RColorBrewer::brewer.pal(4, "Set2"))+NoAxes()+labs(color = "Patient", title = "")
# ct_plot <- DimPlot(lake, group.by = "celltype_aggr", cols = RColorBrewer::brewer.pal(11, "Set3"))+NoAxes()+labs(color = "Lake et al. 2017", title = "")+theme(legend.title = element_text(size = 12))
# age_plot <- DimPlot(lake, group.by = "age", cols = RColorBrewer::brewer.pal(9, "Spectral")[c(6,8,9)])+NoAxes()+labs(color = "", title = "Age (Years)")+theme(plot.title = element_text(size = 12, face = "bold")) #RColorBrewer::brewer.pal(9, "YlOrRd")[c(2,5,7)]
# sex_plot <- DimPlot(lake, group.by = "sex", cols = c(Male = "#73b3b2", Female = "#dfb5c6"))+NoAxes()+labs(color = "Sex", title = "")
# pdf("out/lake_dimplots.pdf", w = 20, h = 4)
# cowplot::plot_grid(sample_plot, pt_plot, age_plot, sex_plot, ct_plot, ncol = 5)
# dev.off()

# featureplots
# Idents(lake) <- lake$celltype
# lake_mks <- FindAllMarkers(lake, only.pos = T) %>% dplyr::mutate(unique = pct.1 - pct.2)




# DimPlots (panel N) ----------------------------------------------------------------

# features
feats <- c("PAX2", "SST", "PVALB")

# clean ages
ages <- unique(so$age) %>% stringr::str_replace(" PCW", "") %>% as.numeric() %>% c(., unique(lin$Age)) %>% unique() %>% sort()

age_bins <- c("5-6", "7-8", "9-10", "11-12", "13-14", "15-16", "17-18", "19-20", "21-22")
age_dict <- data.frame(old = ages) %>% dplyr::mutate(new = cut(old, breaks = c(5,7,9,11,13,15,17,19,21,23), right = FALSE, labels = age_bins) %>% as.character())
so@meta.data %<>% dplyr::mutate(age_clean = age %>% stringr::str_replace(" PCW", ""),
                                age_binned = plyr::mapvalues(age_clean, age_dict$old, age_dict$new, warn_missing = FALSE) %>% factor(levels = age_bins))
lin@meta.data %<>% dplyr::mutate(age_binned = plyr::mapvalues(Age, age_dict$old, age_dict$new, warn_missing = FALSE) %>% factor(levels = age_bins))
age_pal_new <- colorRampPalette(RColorBrewer::brewer.pal(9, "Spectral"))(length(age_bins)) %>% magrittr::set_names(as.character(age_bins))

# params
title_size <- 15
legend_title_size <- 15

# lake
DefaultAssay(lake) <- "RNA"
lake_ct <- DimPlot(lake, group.by = "celltype_pc_only", raster = T, cols = c(tableau20[1], "lightgrey"))+NoAxes()+labs(color = "", title = "Lake et al. 2017")+theme(plot.title = element_text(size = title_size, face = "bold"))
lake_age <- DimPlot(lake, group.by = "age", raster = T, cols = RColorBrewer::brewer.pal(11, "BrBG")[c(8,4,2)])+NoAxes()+labs(color = "Years", title = "Age")+theme(plot.title = element_text(size = title_size, face = "bold")) #RColorBrewer::brewer.pal(9, "YlOrRd")[c(2,5,7)]
lake_fps <- FeaturePlot(lake, feats, order = T, raster = T, cols = c("cadetblue1", "deeppink3"), ncol = 1) & NoAxes() & theme( legend.key.size = unit(8, "pt"))
lake_dot <- DotPlot(lake, features = feats, assay = "RNA", group.by = "celltype_pc_only", cols = c("cadetblue1", "deeppink3")) +
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
  guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white")))+
  labs(x="", y="")+
  # theme(legend.position = "bottom", legend.box = "vertical")+
  theme(legend.position = "none")+
  coord_flip()

# aldinger
DefaultAssay(so) <- "RNA"
ald_ct <- DimPlot(so, group.by = "celltype_aggr", raster = T, cols = c(tableau20[c(1:2)], "lightgrey"))+NoAxes()+labs(color = "", title = "Aldinger et al. 2021")+theme(plot.title = element_text(size = title_size, face = "bold"))
ald_age <- DimPlot(so, group.by = "age_binned", raster = T, cols = age_pal_new)+NoAxes()+labs(color = "PCW", title = "Age")+theme(plot.title = element_text(size = title_size, face = "bold"))+guides(color=guide_legend(ncol=1, override.aes = list(size=2.5)))
ald_fps <- FeaturePlot(so, feats, order = T, raster = T, cols = c("cadetblue1", "deeppink3"), ncol = 1) & NoAxes() & theme(legend.key.size = unit(8, "pt"))
ald_dot <- DotPlot(so, features = feats, assay = "RNA", group.by = "celltype_aggr", cols = c("cadetblue1", "deeppink3")) +
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
  guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white")))+
  labs(x="", y="")+
  # theme(legend.position = "bottom", legend.box = "vertical")+
  theme(legend.position = "none")+
  coord_flip()

# linnarsson
DefaultAssay(lin) <- "RNA"
lin_ct <- DimPlot(lin, group.by = "celltype_aggr", raster = T, cols = c(tableau20[c(1:2)], "lightgrey"))+NoAxes()+labs(color = "", title = "Braun et al. 2022")+theme(plot.title = element_text(size = title_size, face = "bold"))
lin_age <- DimPlot(lin, group.by = "age_binned", raster = T, cols = age_pal_new)+NoAxes()+labs(color = "PCW", title = "Age")+theme(plot.title = element_text(size = title_size, face = "bold"))
lin_fps <- FeaturePlot(lin, c(feats), order = T, raster = T, cols = c("cadetblue1", "deeppink3"), ncol = 1) & NoAxes() & theme(legend.key.size = unit(8, "pt"))
lin_dot <- DotPlot(lin, features = feats, assay = "RNA", group.by = "celltype_aggr", cols = c("cadetblue1", "deeppink3")) +
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
  guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white")))+
  labs(x="", y="")+
  # theme(legend.position = "bottom", legend.box = "vertical")+
  theme(legend.position = "none")+
  coord_flip()

# plot all
pdf(glue::glue("out/fig3_feature_and_dotplots_{paste(feats, collapse = '_')}.pdf"), h = 10, w = 10)
(((lin_ct / lin_age) | (ald_ct / ald_age) | (lake_ct / lake_age)) /
  (lin_dot2 | ald_dot2 | lake_dot2))+patchwork::plot_layout(heights = c(1, 0.5))
dev.off()

# dotplots redone ---
feats <- c("PAX2", "SST", "PVALB")
DefaultAssay(lake) <- "SCT"
DefaultAssay(so) <- "SCT"
DefaultAssay(lin) <- "RNA" # SCT not available for this obj...too many cells.

lake_dot2 <- DotPlot(lake, features = feats, group.by = "celltype_pc_only", cols = c("cadetblue1", "deeppink3")) +
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
  guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white")))+
  labs(x="", y="")+
  theme(legend.position = "bottom", legend.box = "vertical")+
  coord_flip()
ald_dot2 <- DotPlot(so, features = feats,  group.by = "celltype_aggr", cols = c("cadetblue1", "deeppink3")) +
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
  guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white")))+
  labs(x="", y="")+
  theme(legend.position = "bottom", legend.box = "vertical")+
  coord_flip()
lin_dot2 <- DotPlot(lin, features = feats, group.by = "celltype_aggr", cols = c("cadetblue1", "deeppink3")) +
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
  guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white")))+
  labs(x="", y="")+
  theme(legend.position = "bottom", legend.box = "vertical")+
  coord_flip()


# plot for legends
pdf(glue::glue("out/fig3_feature_and_dotplots_{paste(feats, collapse = '_')}_legend_sct.pdf"), h = 3.33, w = 15)
lin_dot2 | ald_dot2 | lake_dot2
dev.off()
 
# plot all (OLD)
# pdf(glue::glue("out/fig2_featureplots_{paste(feats, collapse = '_')}.pdf"), h = 15, w = 10)
# (lin_ct / lin_age / lin_fps + patchwork::plot_layout(heights = c(1,1,3))) |
#   (ald_ct / ald_age / ald_fps + patchwork::plot_layout(heights = c(1,1,3))) |
#   (lake_ct / lake_age / lake_fps + patchwork::plot_layout(heights = c(1,1,3)))
# dev.off()


# extra featureplots --------------------------------------------------------

ald <- readRDS(glue::glue("{yml$ald_path}/cbl_integrated_cleanCC_210111.rds"))
lake <- qs::qread("src/data/lake_cerebellum.qs")
br <- qs::qread(glue::glue("{yml$lin_dir}/hb_cb_final_v3.qs"))

fp_fig_2 <- function(so, feats){ 
  tmp <- feats %>% .[. %in% rownames(so)]
  return(FeaturePlot(so, tmp, raster = T, ncol = length(feats)) & NoAxes() & NoLegend())
}

feats <- c("SORCS3", "NXPH1", "GAD2", "GRM8") %>% sort()
pdf(glue::glue("out/fig2_extra_features.pdf"), w = 4 * length(feats), h = 12)
print(fp_fig_2(ald, feats) / fp_fig_2(lake, feats) / fp_fig_2(br, feats))
dev.off()

# FeatureScatter loop --------------------------------------------------------


fs_out <- "out/featurescatters"
dir.create(fs_out, showWarnings = F, recursive = T)
DefaultAssay(so) <- DefaultAssay(lin) <- DefaultAssay(lake) <- "RNA"
iters <- c("PAX2", "SKOR2", "MKI67", "SOX2") %>% combn(m = 2) %>% as.data.frame() %>% as.list()

# ald -->
so$null_grp <- "null"
for(i in iters){
  
  pdf(glue::glue("{fs_out}/ald_{paste(i, collapse = '_')}.pdf"), h = 5, w = 5)
  print(FeatureScatter(so, i[1], i[2], group.by = "null_grp", raster = T, slot = "counts", pt.size = 0) + NoLegend() + ggtitle("") + geom_jitter())
  dev.off()

}

# lin -->
lin$null_grp <- "null"
for(i in iters){
  
  pdf(glue::glue("{fs_out}/lin_{paste(i, collapse = '_')}.pdf"), h = 5, w = 5)
  print(FeatureScatter(lin, i[1], i[2], group.by = "null_grp", raster = T, slot = "counts", pt.size = 0) + NoLegend() + ggtitle("") + geom_jitter())
  dev.off()

}

# lake -->
lake$null_grp <- "null"
for(i in iters){
  
  tryCatch({

    pdf(glue::glue("{fs_out}/lake_{paste(i, collapse = '_')}.pdf"), h = 5, w = 5)
    print(FeatureScatter(lake, i[1], i[2], group.by = "null_grp", raster = T, slot = "counts", pt.size = 0) + NoLegend() + ggtitle("") + geom_jitter())
    dev.off()

  }, error = function(e) { NA })
  

}

# FeatureScatters SKOR2 PAX2 (cleaned...panel O) ----------------------------

# import dev mouse data
mv_dir <- paste0(yml$mv_dir, "/original_manuscript_mouse.rds")
mv <- readRDS(mv_dir)
mm_lin <- qs::qread("out/linnarsson_hindbrain_mouse.qs")

# assign groups
so@meta.data %<>% dplyr::mutate(skor2_pax2 = ifelse(rownames(.) %in% WhichCells(so, expression = SKOR2 > 0 & PAX2 > 0), "double_pos", "not"))
lin@meta.data %<>% dplyr::mutate(skor2_pax2 = ifelse(rownames(.) %in% WhichCells(lin, expression = SKOR2 > 0 & PAX2 > 0), "double_pos", "not"))
mv@meta.data %<>% dplyr::mutate(skor2_pax2 = ifelse(rownames(.) %in% WhichCells(mv, expression = Skor2 > 0 & Pax2 > 0), "double_pos", "not"))
mm_lin@meta.data %<>% dplyr::mutate(skor2_pax2 = ifelse(rownames(.) %in% WhichCells(mm_lin, expression = Skor2 > 0 & Pax2 > 0), "double_pos", "not"))

# plots
jit <- TRUE
pt_size <- 0.5
grp <- "skor2_pax2"
pal <- c("deeppink3", "cadetblue1")
fs1 <- FeatureScatter(lin, "SKOR2", "PAX2", group.by = grp, raster = T, slot = "counts", pt.size = pt_size, jitter = jit, cols = pal) + NoLegend() + ggtitle("Braun et al. 2022")
fs2 <- FeatureScatter(so, "SKOR2", "PAX2", group.by = grp, raster = T, slot = "counts", pt.size = pt_size, jitter = jit, cols = pal) + NoLegend() + ggtitle("Aldinger et al. 2021")
fs3 <- FeatureScatter(mv, "Skor2", "Pax2", group.by = grp, raster = T, slot = "counts", pt.size = pt_size, jitter = jit, cols = pal) + NoLegend() + ggtitle("Vladoiu et al. 2019")
fs4 <- FeatureScatter(mm_lin, "Skor2", "Pax2", group.by = grp, raster = T, slot = "counts", pt.size = pt_size, jitter = jit, cols = pal) + NoLegend() + ggtitle("La Manno et al. 2021")

# print
pdf(glue::glue("{fs_out}/all_skor2_pax2.pdf"), h = 5, w = 20)
print(fs1 | fs2 | fs3 | fs4)
dev.off()


# stacked barplot quantification
ald_skor2 <- WhichCells(so, expression = SKOR2 > 0 & PAX2 == 0) %>% length()
ald_pax2 <- WhichCells(so, expression = SKOR2 == 0 & PAX2 > 0) %>% length()
ald_skor2_pax2 <- WhichCells(so, expression = SKOR2 > 0 & PAX2 > 0) %>% length()
braun_skor2 <- WhichCells(lin, expression = SKOR2 > 0 & PAX2 == 0) %>% length()
braun_pax2 <- WhichCells(lin, expression = SKOR2 == 0 & PAX2 > 0) %>% length()
braun_skor2_pax2 <- WhichCells(lin, expression = SKOR2 > 0 & PAX2 > 0) %>% length()
vlad_skor2 <- WhichCells(mv, expression = Skor2 > 0 & Pax2 == 0) %>% length()
vlad_pax2 <- WhichCells(mv, expression = Skor2 == 0 & Pax2 > 0) %>% length()
vlad_skor2_pax2 <- WhichCells(mv, expression = Skor2 > 0 & Pax2 > 0) %>% length()
lam_skor2 <- WhichCells(mm_lin, expression = Skor2 > 0 & Pax2 == 0) %>% length()
lam_pax2 <- WhichCells(mm_lin, expression = Skor2 == 0 & Pax2 > 0) %>% length()
lam_skor2_pax2 <- WhichCells(mm_lin, expression = Skor2 > 0 & Pax2 > 0) %>% length()

bpdf <- data.frame("ald" = c(ald_skor2, ald_pax2, ald_skor2_pax2),
                    "braun" = c(braun_skor2, braun_pax2, braun_skor2_pax2),
                    "vlad" = c(vlad_skor2, vlad_pax2, vlad_skor2_pax2),
                    "lam" = c(lam_skor2, lam_pax2, lam_skor2_pax2)) %>% 
  magrittr::set_rownames(c("skor2", "pax2", "skor2_pax2")) %>% 
  t() %>% 
  tibble::as_tibble(rownames = "study") %>% 
  tidyr::pivot_longer(-study, names_to = "marker", values_to = "count") %>% 
  dplyr::mutate(study = factor(study, levels = c("braun", "ald", "lam", "vlad")))

ggbp <- ggplot(bpdf, aes(study, count, fill = marker))+
  geom_bar(position = "fill", stat = "identity")+ #position = "stack"
  scale_fill_manual(values = c("cadetblue1", "cadetblue4", "deeppink3"), labels = c("SKOR2+", "PAX2+", "SKOR2+ PAX2+"))+
  scale_x_discrete(labels=c("Braun", "Aldinger", "La Manno", "Vladoiu"))+
  theme_classic()+
  labs(y = "Cells (%)", x = "", fill = "")+
  theme(axis.text.x=element_text(colour="black"),
    axis.text.y=element_text(colour="black"))

pdf(glue::glue("{fs_out}/all_skor2_pax2_barplot.pdf"), h = 4, w = 5)
print(ggbp)
dev.off()

pdf(glue::glue("{fs_out}/all_skor2_pax2_barplot_scatters.pdf"), h = 4, w = 27)
print((fs1 | fs2 | fs4 | fs3 | ggbp) + patchwork::plot_layout(widths = c(1,1,1,1,1.2)))
dev.off()


# DotPlots (OLD) ----------------------------------------------------------------

# gene_list <- c("PAX2", "PVALB", "SST", "ITPR1")
# 
# # lake
# lake_dp <- DotPlot(lake, 
#         group.by = "celltype",
#         assay = "SCT", 
#         features = gene_list)+
#   scale_colour_gradient2(low = mypal[3], mid = mypal[2], high = mypal[1], breaks = c(-3,0,3), limits = c(-3,3))+
#   # NoLegend()+
#   scale_size(breaks = c(0,30,60), limits = c(0,60))+
#   labs(x="", y="")+
#   guides(color = guide_colorbar(title = "Expr. (SCT)", title.position = "top"),
#          size = guide_legend(title = "Pct. Expr.", title.position = "top"))
#   
# # ald
# so@meta.data %<>% dplyr::mutate(figure_clusters = factor(figure_clusters, levels = rev(sort(unique(so$figure_clusters)))))
# ald_dp <- DotPlot(so, 
#         group.by = "figure_clusters",
#         assay = "SCT", 
#         features = gene_list)+
#   scale_colour_gradient2(low = mypal[3], mid = mypal[2], high = mypal[1], breaks = c(-3,0,3), limits = c(-3,3))+
#   labs(x="", y="")+
#   guides(color = guide_colorbar(title = "Expr. (SCT)", title.position = "top"),
#          size = guide_legend(title = "Pct. Expr.", title.position = "top"))
# 
# # print
# ald_dp / lake_dp + patchwork::plot_layout(guides = "collect") &
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
#         legend.position = "bottom",
#         legend.title = element_text(size = 12, hjust = 0.5))


# PAX2 and SKOR2 fractions for SOX2 and MKI67 ----------------------------------
# pax2_ald <- so@assays$RNA@counts["PAX2", ] %>% .[. != 0] %>% names()