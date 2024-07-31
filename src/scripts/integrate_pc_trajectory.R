# Info --------------------------------------------------------------------
# Anders E.
# May 25, 2024
# Taylor lab

# Notes ---------------------------------------------------------------
# 
# 
# purpose of this script is to use RPCA to integrate all available
# normal PC datasets.
# 
# thoughts on cao dataset:
# downsample because the cao object is so relatively big (365k cells of 483k total)
# 
# 
# thoughts on luo dataset:
# tread carefully,
# because the luo object has multiple clusters which do not integrate at all even internally
# and because it is lacking entries (not just 0 expression, I mean the whole row is missing) for many key genes e.g. MKI67, MEIS2, MECOM, NEUROD1???
# 



# Libraries ---------------------------------------------------------------
# conda activate ../G4MB/src/envs/r_general
# R

# libraries
library(qs)
library(Seurat)
library(Signac)
library(magrittr)
library(ggplot2)
source("src/scripts/utils.R")
yml <- yaml::read_yaml("config.yaml")

# outdir
outdir <- "out/integration/pc"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# method
meth <- "norm"


# Import CB objects ---------------------------------------------------------------

mannens <- qs::qread("../G4MB/out/integration/mannens_multiome/integrated_10X280_1_ABCD_1___10X346_1_ABCD_1___10X365_2_ABCDE_2___10X406_1_ABCD_2.qs")
aldinger <- readRDS(glue::glue("{yml$ald_path}/cbl_integrated_cleanCC_210111.rds"))
braun <- qs::qread("../G4MB/out/braun/braun_cb_microenv.qs")
cao <- qs::qread(glue::glue("{yml$public_scrnaseq_dir}/2021_Cao_Science/cao_cerebellum.qs"))
sepp <- qs::qread("../G4MB/out/sepp/sepp_hgnc.qs")
# luo <- readRDS(glue::glue("{yml$public_scrna}/luo_2022_nature/GSM5952337_Fetal_cerebellum_final_cluster_1.rds")) # omit due to ... concerns
siletti <- qs::qread(glue::glue("{yml$public_scrna}/siletti_2022_biorxiv/sil_cerebellum.qs"))
ts <- qs::qread("out/cs20s_annot.qs") # ts  = "this study"


# Subset to GABA lineages ---------------------------------------------------------------

# initial look function
il <- function(so, label, red){

    if(label == "sepp"){ DefaultAssay(so) <- "hgnc" } else { DefaultAssay(so) <- "RNA" }
    dim <- DimPlot(so, raster = T, label = T, repel = T, reduction = red) & NoAxes() & NoLegend() & ggtitle("")
    feat <- FeaturePlot(so, c("SOX2", "WLS", "PTF1A", "PRDM13", "PAX2", "SKOR2"), raster = T, order = T, reduction = red, ncol = 6) & NoLegend() & NoAxes()
    pdf(glue::glue("{outdir}/prelim_{label}.pdf"), w = 35, h = 5)
    print((dim | feat) + patchwork::plot_layout(widths = c(1,6)))
    dev.off()

}

# recluster function
recluster <- function(so, label){
    so %>% 
        NormalizeData(normalization.method = "LogNormalize", scale.factor = 10000) %>% 
        FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>% 
        ScaleData() %>% 
        RunPCA() %>% 
        FindNeighbors(dims = 1:20) %>%
        FindClusters(resolution = 0.2) %>% 
        RunUMAP(dims = 1:20) %>% 
        AddMetaData(label, col.name = "sample_id") %>% 
        return()
}


# mannens
il(mannens, "mannens", "wnn.umap.recip")
DefaultAssay(mannens) <- "RNA"

mannens_ss <- mannens %>%
    subset(subset = seurat_clusters %in% c(2,8,9,4,6,16,3,0)) %>% 
    recluster("mannens")
qs::qsave(mannens_ss, glue::glue("{outdir}/mannens_ss.qs"))

# aldinger
il(aldinger, "aldinger", "umap")
DefaultAssay(aldinger) <- "RNA"

aldinger_ss <- aldinger %>%
    subset(subset = figure_clusters %in% c("07-PIP", "01-PC")) %>% 
    recluster("aldinger")
qs::qsave(aldinger_ss, glue::glue("{outdir}/aldinger_ss.qs"))

# braun
il(braun, "braun", "umap")
DefaultAssay(braun) <- "RNA"

braun_ss <- braun %>%
    subset(subset = seurat_clusters %in% c(0,1,6,4,9,12,8,2)) %>% 
    recluster("braun")
qs::qsave(braun_ss, glue::glue("{outdir}/braun_ss.qs"))

# cao
il(cao, "cao", "umap")
DefaultAssay(cao) <- "RNA"

cao_ss <- cao %>%
    subset(subset = Main_cluster_name %in% c("Purkinje neurons", "Inhibitory interneurons")) %>% 
    recluster("cao")
qs::qsave(cao_ss, glue::glue("{outdir}/cao_ss.qs"))

# sepp
sepp %<>% recluster("sepp")
il(sepp, "sepp", "umap")
DefaultAssay(sepp) <- "hgnc"

pdf(glue::glue("{outdir}/sepp_cell_types.pdf"), w = 12, h = 10)
DimPlot(sepp, group.by = "author_cell_type", raster = T, label = T) & NoAxes() & ggtitle("")
dev.off()

sepp_ss <- sepp %>% 
    subset(subset = author_cell_type %in% c("Purkinje", "VZ_neuroblast", "interneuron")) %>% 
    recluster("sepp")
qs::qsave(sepp_ss, glue::glue("{outdir}/sepp_ss.qs"))

# luo
# il(luo, "luo", "umap")
# luo_ss <- luo %>% 
#     subset(subset = celltype %in% c("01xNSC", "03xDev.purkinje", "05xGABA_interneuron", "04xPurkinje")) %>% 
#     recluster("luo")
# qs::qsave(luo_ss, glue::glue("{outdir}/luo_ss.qs"))

# siletti
il(siletti, "siletti", "umap")
siletti_ss <- siletti %>% 
    subset(subset = celltype %in% c("GABA (NXPH1+)", "GABA DCN (SOX14+)", "PC (PCP4+)", "PIP (PAX2+)")) %>% 
    recluster("siletti")
qs::qsave(siletti_ss, glue::glue("{outdir}/siletti_ss.qs"))


# this study
il(ts, "ts", "umap")
ts_ss <- ts %>% 
    subset(subset = celltype_broad %in% c("NE & VZ", "Early PC", "Late PC")) %>% 
    recluster("ts")
qs::qsave(ts_ss, glue::glue("{outdir}/ts_ss.qs"))


# Reread subsets & downsample for balanced integration ---------------------------------------------------------------

# loop reread
id_list <- c("mannens", "aldinger", "braun", "cao", "sepp", "siletti", "ts") # "luo"

rna_list <- purrr::map(id_list, function(x){
  glue::glue("{outdir}/{x}_ss.qs") %>% 
    qs::qread() %>% 
    return()
}) %>% magrittr::set_names(id_list)

# downsample to same number as in aldinger...where applicable
set.seed(666)
n_aldinger <- rna_list[["aldinger"]]@meta.data %>% nrow()
cao_cells <- rna_list[["cao"]]@meta.data %>% dplyr::slice_sample(n = n_aldinger) %>% rownames()
sepp_cells <- rna_list[["sepp"]]@meta.data %>% dplyr::slice_sample(n = n_aldinger) %>% rownames()
braun_cells <- rna_list[["braun"]]@meta.data %>% dplyr::slice_sample(n = n_aldinger) %>% rownames()

rna_list[["cao"]] <- rna_list[["cao"]] %>% subset(cells = cao_cells)
rna_list[["sepp"]] <- rna_list[["sepp"]] %>% subset(cells = sepp_cells)
rna_list[["braun"]] <- rna_list[["braun"]] %>% subset(cells = braun_cells)

# more cleaning
rna_list[["sepp"]][["RNA"]] <- rna_list[["sepp"]][["hgnc"]]

# Integrate ---------------------------------------------------------------


# re-normalize
rna_list %<>% purrr::map(\(x) { 
    DefaultAssay(x) <- "RNA"
    x %>% 
    DietSeurat(counts = TRUE, data = TRUE, assays = "RNA") %>% 
    NormalizeData() %>% 
    FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>% 
    return()
}) %>% magrittr::set_names(id_list)

# get integration features
int_features <- SelectIntegrationFeatures(object.list = rna_list)

# pca using integration features
rna_list %<>% purrr::map(\(x) { x %>% ScaleData(features = int_features) %>% RunPCA(features = int_features) %>% return() }) %>% 
    magrittr::set_names(id_list)

# find anchors
int_anchors <- FindIntegrationAnchors(object.list = rna_list, anchor.features = int_features, reduction = "rpca")

# integrate
int_rna <- IntegrateData(anchorset = int_anchors)

DefaultAssay(int_rna) <- "integrated"

# process
int_rna %<>% 
    ScaleData() %>%
    RunPCA() %>% 
    RunUMAP(dims = 1:20) %>% 
    FindNeighbors(dims = 1:20) %>%
    FindClusters(resolution = 0.2)

# save
qs::qsave(int_rna, glue::glue("{outdir}/pc_integrated.qs"))


# Clean integrated object ---------------------------------------------------------------

# reload
pc <- qs::qread(glue::glue("{outdir}/pc_integrated.qs"))

# cc
DefaultAssay(pc) <- "RNA"
pc %<>% CellCycleScoring(s.features = cc.genes.updated.2019$s.genes, g2m.features = cc.genes.updated.2019$g2m.genes, set.ident = FALSE)
pc$cc_diff <- pc$S.Score - pc$G2M.Score

# plot
d <- DimPlot(pc, group.by = c("seurat_clusters", "sample_id", "Phase"), raster = T, label = T, ncol = 3) & NoAxes() & ggtitle("")
f <- FeaturePlot(pc, c("SOX2", "WLS", "PTF1A", "PRDM13", "PAX2", "SKOR2"), raster = T, ncol = 3) & NoLegend() & NoAxes()
pdf(glue::glue("{outdir}/pc.pdf"), width = 17, height = 15)
(d / f) + patchwork::plot_layout(heights = c(1,2))
dev.off()
pdf(glue::glue("{outdir}/pc_feats.pdf"), width = 25, height = 85)
FeaturePlot(pc, genes, raster = T, ncol = 5) & NoAxes() & NoLegend()
dev.off()

# subset
DefaultAssay(pc) <- "integrated"
ss <- pc %>% 
    subset(subset = seurat_clusters %in% setdiff(unique(as.numeric(pc$seurat_clusters)), c(3,7,11,12))) %>% # rhombic lip and glial celltypes or sample-specific clusters
    FindVariableFeatures(selection.method = "vst") %>% 
    ScaleData() %>% 
    RunPCA() %>% 
    FindNeighbors(dims = 1:20) %>%
    FindClusters(resolution = 0.2) %>% 
    RunUMAP(dims = 1:20)

# replot
dss <- DimPlot(ss, group.by = c("seurat_clusters", "sample_id", "Phase"), raster = T, label = T, ncol = 3) & NoAxes() & ggtitle("")
fss <- FeaturePlot(ss, c("SOX2", "WLS", "PTF1A", "PRDM13", "PAX2", "SKOR2", "SOX14"), raster = T, ncol = 3) & NoLegend() & NoAxes()
pdf(glue::glue("{outdir}/ss.pdf"), width = 17, height = 20)
(dss / fss) + patchwork::plot_layout(heights = c(1,3))
dev.off()
DefaultAssay(ss) <- "RNA"
pdf(glue::glue("{outdir}/ss_feats.pdf"), width = 25, height = 90)
FeaturePlot(ss, sort(c(genes, "KCNA2", 'GRIK3', 'NXPH1', 'SLC6A5', 'ALDH1A3', 'LGI2', 'GJD2', 'SORCS3')), raster = T, ncol = 5) & NoAxes() & NoLegend()
dev.off()
DefaultAssay(ss) <- "integrated"
qs::qsave(ss, glue::glue("{outdir}/ss.qs"))


mks <- FindAllMarkers(ss, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.75) %>% dplyr::mutate(unique = pct.1 - pct.2)
mks %>% dplyr::group_by(cluster) %>% dplyr::slice_max(avg_log2FC, n = 5) %>% as.data.frame()
writexl::write_xlsx(mks, glue::glue("{outdir}/pc_subset_markers.xlsx"))



# quick save
qs::qsave(ss, glue::glue("{outdir}/ss.qs"))
ss <- qs::qread(glue::glue("{outdir}/ss.qs"))


# remove cl 10 --> low quality cells, uninformative markers
DefaultAssay(ss) <- "integrated"
sss <- ss %>% 
    subset(subset = seurat_clusters != "10") %>%
    FindVariableFeatures(selection.method = "vst") %>% 
    ScaleData() %>% 
    RunPCA() %>% 
    FindNeighbors(dims = 1:20) %>%
    FindClusters(resolution = 0.2) %>% 
    RunUMAP(dims = 1:20)

# replot
dsss <- DimPlot(sss, group.by = c("seurat_clusters", "sample_id", "Phase"), raster = T, label = T, ncol = 3) & NoAxes() & ggtitle("")
fsss <- FeaturePlot(sss, c("SOX2", "WLS", "PTF1A", "PRDM13", "PAX2", "SKOR2", "SOX14"), raster = T, ncol = 3) & NoLegend() & NoAxes()
pdf(glue::glue("{outdir}/sss.pdf"), width = 17, height = 20)
(dsss / fsss) + patchwork::plot_layout(heights = c(1,3))
dev.off()
DefaultAssay(sss) <- "RNA"
pdf(glue::glue("{outdir}/sss_feats.pdf"), width = 25, height = 90)
FeaturePlot(sss, sort(c(genes, "KCNA2", 'GRIK3', 'NXPH1', 'SLC6A5', 'ALDH1A3', 'LGI2', 'GJD2', 'SORCS3')), raster = T, ncol = 5) & NoAxes() & NoLegend()
dev.off()

# quick save
qs::qsave(sss, glue::glue("{outdir}/sss.qs"))
sss <- qs::qread(glue::glue("{outdir}/sss.qs"))



mannens_age_dict <- sss@meta.data %>% 
    tibble::rownames_to_column("barcode") %>% 
    dplyr::filter(sample_id == "mannens") %>% 
    tidyr::separate(barcode, into = c("key1", "key2", "key3", "key4", "key5"), sep = "_") %>% 
    dplyr::select(key1, Age) %>% 
    dplyr::distinct() %>% 
    dplyr::filter(!is.na(Age))

# annotate & clean more metadata
sss@meta.data %<>% dplyr::mutate(celltype = dplyr::case_when(
    seurat_clusters %in% c(5) ~ "NE",
    seurat_clusters %in% c(1) ~ "PC 1",
    seurat_clusters %in% c(2) ~ "PC 2",
    seurat_clusters %in% c(0) ~ "PC 3",
    seurat_clusters %in% c(7) ~ "PC 4",
    seurat_clusters %in% c(4) ~ "PIP 1",
    seurat_clusters %in% c(3,10) ~ "PIP 2",
    seurat_clusters %in% c(9) ~ "IN 1",
    seurat_clusters %in% c(8) ~ "IN 2",
    seurat_clusters %in% c(6) ~ "IN 3"),
    mannens_key = stringr::str_extract(rownames(.), "10X[0-9]+"),
    Age = ifelse(sample_id == "mannens" & is.na(Age), plyr::mapvalues(mannens_key, mannens_age_dict$key1, mannens_age_dict$Age, warn_missing = FALSE), Age),
    age_aggr_raw = dplyr::case_when( 
        grepl("PCW", batch) ~ batch,
        !is.na(age) ~ age,
        !is.na(Age) ~ as.character(Age),
        !is.na(Development_day) ~ as.character(round(as.numeric(Development_day) / 7, 1)),
        !is.na(development_stage) ~ development_stage,
        grepl("CS21_|CS22_|CS23_", rownames(.)) ~ "8",
        TRUE ~ "NA"),
    age_aggr_numeral = age_aggr_raw %>% 
        stringr::str_replace(" PCW", "") %>% 
        stringr::str_replace("th week post-fertilization human stage", "") %>% 
        stringr::str_replace("^PCW", "") %>% 
        stringr::str_replace("Carnegie stage 18|Carnegie stage 19", "6") %>%  # we will be rounding all DOWN to the nearest week
        stringr::str_replace("Carnegie stage 22", "7") %>% 
        stringr::str_replace("newborn human stage", "38") %>% 
        stringr::str_replace("-year-old human stage|-month-old human stage", "") %>% 
        as.numeric() %>% 
        floor(),
    age_aggr = dplyr::case_when(
        grepl("year|42|29|50", age_aggr_raw) ~ paste0(as.character(age_aggr_numeral), " yrs"),
        grepl("month", age_aggr_raw) ~ paste0(as.character(age_aggr_numeral), " mths"),
        TRUE ~ paste0(as.character(age_aggr_numeral), " PCW")
    ))

lvls <- c("6 PCW", "7 PCW", "8 PCW", "9 PCW", "10 PCW", "11 PCW", "12 PCW", "13 PCW", "14 PCW", "15 PCW", "16 PCW", "17 PCW", "18 PCW", "20 PCW", "21 PCW", "38 PCW", "6 mths", "7 mths", "9 mths", "3 yrs", "29 yrs", "42 yrs", "44 yrs", "46 yrs", "50 yrs", "52 yrs")
sss$age_factor <- factor(sss$age_aggr, levels = lvls)

celltype_lvls <- c("NE", "PC 1", "PC 2", "PC 3", "PC 4", "PIP 1", "PIP 2", "IN 1", "IN 2", "IN 3")
sss$celltype_factor <- factor(sss$celltype, levels = celltype_lvls)


# clean study and sample IDs
sss$study_id <- sss$sample_id %>% stringr::str_replace('ts', 'this study') %>% snakecase::to_sentence_case()
sss$sample_id_final <- dplyr::coalesce(sss$donor_id, sss$batch, sss$Fetus_id, sss$experiment, sss$mannens_key, sss$orig.ident)

# plot
age_pal <- colorRampPalette(RColorBrewer::brewer.pal(11, 'Spectral'))(length(unique(lvls)))
celltype_pal <- tableau10 %>% magrittr::set_names(celltype_lvls)

sssp8 <- DimPlot(sss, group.by = c("celltype_factor"), reduction = "umap", raster = TRUE, cols = celltype_pal) & NoAxes() & ggtitle("")
sssp9 <- DimPlot(sss, group.by = "age_factor", reduction = "umap", raster = TRUE, cols = age_pal) & NoAxes() & ggtitle("")
sssp10 <- DimPlot(sss, group.by = "Phase", reduction = "umap", raster = TRUE, cols = phase_pal) & NoAxes() & ggtitle("")
sssp11 <- DimPlot(sss, group.by = "study_id", reduction = "umap", raster = TRUE, cols = tableau10) & NoAxes() & ggtitle("")
sssp12 <- DimPlot(sss, group.by = "sample_id_final", reduction = "umap", raster = TRUE) & NoAxes() & NoLegend() & ggtitle("")
pdf(glue::glue("{outdir}/final_object/umaps.pdf"), width = 27, height = 5)
sssp8 | sssp9 | sssp10 | sssp11 | sssp12
dev.off()

# save
dir.create(glue::glue("{outdir}/final_object"), showWarnings = FALSE, recursive = TRUE)
qs::qsave(sss, glue::glue("{outdir}/final_object/integrated.qs"))
sss <- qs::qread(glue::glue("{outdir}/final_object/integrated.qs"))

# lemma
# study- and sample-
head(sss@meta.data)
names(sss@meta.data)
sss$study_id <- sss$sample_id %>% stringr::str_replace('ts', 'this study') %>% snakecase::to_sentence_case()
table(sss$study_id, useNA = 'always')
system('clear')

table(sss$donor_id, useNA = 'always')
table(sss$batch, useNA = 'always')
table(sss$Fetus_id, useNA = 'always')
table(sss$experiment, useNA = 'always')
table(sss$sample, useNA = 'always')

table(sss$sample_id_final, useNA = 'always')


# end lemma

# features
DefaultAssay(sss) <- "RNA"
bonus_genes <- c("HES5", "KIRREL2", "MKI67", "NOTCH1", "PIEZO2", "SHROOM3", "SOX2", "TCF7L1", "TCF7L2", "VIM", "LHX1", "PDRM13", "PTF1A", "CALB1", "EBF1", "EBF2", "ESRRB", "FOXP1", "FOXP2", "PCP4" , "SKOR2", "GRIK3", "KCNA2", "NXPH1", "PAX2", "PVALB", "SLC6A5", "SORCS3", 'GRIK3', 'NXPH1', 'SLC6A5', 'ALDH1A3', 'LGI2', 'GJD2', 'SORCS3', "KCNA2") %>% unique()
to_plot <- sort(unique(c(genes, bonus_genes)))
pdf(glue::glue("{outdir}/final_object/featureplots.pdf"), width = 25, height = 5 * ceiling(length(to_plot) / 5))
FeaturePlot(sss, to_plot, raster = T, ncol = 5) & NoAxes() & NoLegend()
dev.off()



# Slingshot ---------------------------------------------------------------
# reload
sss <- qs::qread(glue::glue("{outdir}/final_object/integrated.qs"))
pto <- slingshot::slingshot(data = Embeddings(sss, reduction = "umap"),
    clusterLabels = sss$celltype, 
    start.clus = "NE",
    end.clus = c("PC 4", "PIP 2", "IN 1"))
sds <- slingshot::as.SlingshotDataSet(pto)

# see number of curves
slingshot::slingCurves(sds) %>% length()

# build df for ggplot
ggdf <- Embeddings(sss, reduction = "umap") %>%
  as.data.frame() %>% 
  dplyr::mutate(celltype = sss$celltype) %>% 
  cbind(pto@assays@data@listData$pseudotime) %>% 
  tibble::rownames_to_column("barcode") %>% 
  dplyr::mutate(sling_aggr = dplyr::coalesce(Lineage1, Lineage2, Lineage3))#, sling_gc_postnatal)) # setting sling_gc first since longer overall trajectory --> sets RL-VZ as earlier in both cases.

# get curve coords
crv1 <- slingshot::slingCurves(sds)[[1]]$s %>% as.data.frame() 
crv2 <- slingshot::slingCurves(sds)[[2]]$s %>% as.data.frame() 
crv3 <- slingshot::slingCurves(sds)[[3]]$s %>% as.data.frame() 
crv4 <- slingshot::slingCurves(sds)[[4]]$s %>% as.data.frame()
crv5 <- slingshot::slingCurves(sds)[[5]]$s %>% as.data.frame()

# plot
p1 <- ggplot(ggdf, aes(UMAP_1, UMAP_2, color = celltype))+
        geom_point()+
        theme_void()+
        theme(legend.text=element_text(size=10),
              legend.title=element_blank())+
        xlab("UMAP 1")+
        ylab("UMAP 2")+
        scale_colour_manual(values = colors_dutch)+
        geom_path(data = crv1, aes(UMAP_1, UMAP_2), size = 0.75, color = "black")+
        geom_path(data = crv2, aes(UMAP_1, UMAP_2), size = 0.75, color = "black")+
        geom_path(data = crv3, aes(UMAP_1, UMAP_2), size = 0.75, color = "black")+
        geom_path(data = crv4, aes(UMAP_1, UMAP_2), size = 0.75, color = "black")+
        geom_path(data = crv5, aes(UMAP_1, UMAP_2), size = 0.75, color = "black")+
        guides(color = guide_legend(override.aes = list(size = 3)))

p2 <- ggplot(ggdf, aes(UMAP_1, UMAP_2, color = sling_aggr))+
        geom_point()+
        theme_void()+
        theme(legend.text=element_text(size=10),
              legend.title=element_blank())+
        xlab("UMAP 1")+
        ylab("UMAP 2")+
        viridis::scale_color_viridis(option = "rocket")+
        guides(color = guide_legend(override.aes = list(size = 3)))

pdf(glue::glue("{outdir}/final_object/slingshot.pdf"), width = 12, height = 5)
p1 | p2
dev.off()

# save
# DO NOT USE
# ggdf %>% data.table::fwrite(glue::glue("{outdir}/final_object/sling_coords.csv"))

# notes:
# interneuron trajectories p much useless
# one trajectory for PC is sensible, though


# PC pseuodtime trajectory ---------------------------------------------------------------

DefaultAssay(sss) <- "integrated"
pc <- sss %>% 
    subset(subset = celltype %in% c("NE", "PC 1", "PC 2", "PC 3", "PC 4")) %>%
    FindVariableFeatures(selection.method = "vst") %>% 
    ScaleData() %>% 
    RunPCA() %>% 
    FindNeighbors(dims = 1:20) %>%
    FindClusters(resolution = 0.3) %>% 
    RunUMAP(dims = 1:20)

# quick save
qs::qsave(pc, glue::glue("{outdir}/pc.qs"))
pc <- qs::qread(glue::glue("{outdir}/pc.qs"))

# palettes
lvls <- c("6 PCW", "7 PCW", "8 PCW", "9 PCW", "10 PCW", "11 PCW", "12 PCW", "13 PCW", "14 PCW", "15 PCW", "16 PCW", "17 PCW", "18 PCW", "20 PCW", "21 PCW", "38 PCW", "6 mths", "7 mths", "9 mths", "3 yrs", "29 yrs", "42 yrs", "44 yrs", "46 yrs", "50 yrs", "52 yrs")
celltype_lvls <- c("NE", "PC 1", "PC 2", "PC 3", "PC 4", "PIP 1", "PIP 2", "IN 1", "IN 2", "IN 3")
age_pal <- colorRampPalette(RColorBrewer::brewer.pal(11, 'Spectral'))(length(unique(lvls)))
celltype_pal <- tableau10 %>% magrittr::set_names(celltype_lvls)

# plot
pc_plot1 <- DimPlot(pc, group.by = "celltype", reduction = "umap", raster = TRUE, cols = celltype_pal) & NoAxes() & ggtitle("")
pc_plot2 <- DimPlot(pc, group.by = "age_factor", reduction = "umap", raster = TRUE, cols = age_pal) & NoAxes() & ggtitle("")
pc_plot3 <- DimPlot(pc, group.by = "Phase", reduction = "umap", raster = TRUE, cols = phase_pal) & NoAxes() & ggtitle("")
pdf(glue::glue("{outdir}/pc_umaps.pdf"), width = 16, height = 5)
pc_plot1 | pc_plot2 | pc_plot3 
dev.off()

# slingshot
pto <- slingshot::slingshot(data = Embeddings(pc, reduction = "umap"),
    clusterLabels = pc$celltype, 
    start.clus = "NE",
    end.clus = c("PC 4"))
sds <- slingshot::as.SlingshotDataSet(pto)

# see number of curves
slingshot::slingCurves(sds) %>% length()

# build df for ggplot
ggdf <- Embeddings(pc, reduction = "umap") %>%
  as.data.frame() %>% 
  dplyr::mutate(celltype = pc$celltype) %>% 
  cbind(pto@assays@data@listData$pseudotime) %>% 
  tibble::rownames_to_column("barcode")

# get curve coords
crv1 <- slingshot::slingCurves(sds)[[1]]$s %>% as.data.frame() 

# plot
p1 <- ggplot(ggdf, aes(UMAP_1, UMAP_2, color = celltype))+
        geom_point()+
        theme_void()+
        theme(legend.text=element_text(size=10),
              legend.title=element_blank())+
        xlab("UMAP 1")+
        ylab("UMAP 2")+
        scale_colour_manual(values = celltype_pal)+
        geom_path(data = crv1, aes(UMAP_1, UMAP_2), size = 0.75, color = "black")+
        guides(color = guide_legend(override.aes = list(size = 3)))

p2 <- ggplot(ggdf, aes(UMAP_1, UMAP_2, color = Lineage1))+
        geom_point()+
        theme_void()+
        theme(legend.text=element_text(size=10),
              legend.title=element_blank())+
        xlab("UMAP 1")+
        ylab("UMAP 2")+
        viridis::scale_color_viridis(option = "rocket")+
        guides(color = guide_legend(override.aes = list(size = 3)))

pdf(glue::glue("{outdir}/slingshot_pc.pdf"), width = 12, height = 5)
p1 | p2
dev.off()



# markers
Idents(pc) <- pc$celltype
pc_mks <- FindAllMarkers(pc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.5) %>% dplyr::mutate(unique = pct.1 - pct.2)
pc_mks %>% readr::write_csv(glue::glue("{outdir}/pc_markers.csv"))
pc_mks <- readr::read_csv(glue::glue("{outdir}/pc_markers.csv"))

# subset to corresponding marker genes
pc_genes <- pc_mks %>% 
  dplyr::group_by(cluster) %>% 
  dplyr::slice_min(p_val_adj, n = 20) %>% 
  dplyr::slice_max(avg_log2FC, n = 20) %>% 
  dplyr::arrange(factor(cluster, levels = c("NE", "PC 1", "PC 2", "PC 3", "PC 4"))) %>%
  dplyr::pull(gene) %>% 
  unique() %>% c(., "SOX2")
  
DefaultAssay(pc) <- "RNA"
pdf(glue::glue("{outdir}/pc_traj_featureplots.pdf"), width = 25, height = 75)
FeaturePlot(pc, pc_genes, raster = T, ncol = 5) & NoAxes() & NoLegend()
dev.off()

bonus_genes <- c("HES5", "KIRREL2", "MKI67", "NOTCH1", "PIEZO2", "SHROOM3", "SOX2", "TCF7L1", "TCF7L2", "VIM", "LHX1", "PDRM13", "PTF1A", "CALB1", "EBF1", "EBF2", "ESRRB", "FOXP1", "FOXP2", "PCP4" , "SKOR2", "GRIK3", "KCNA2", "NXPH1", "PAX2", "PVALB", "SLC6A5", "SORCS3", 'GRIK3', 'NXPH1', 'SLC6A5', 'ALDH1A3', 'LGI2', 'GJD2', 'SORCS3', "KCNA2",
"SRGAP2C", "SYNGAP1", "GATA3", "SLC1A6", "ALDOC", "PLCB4", "PLCB3", "TRPC3", "GRM1", "ITPR1", "PRKCD", "SCN1A", "SCN2A", "SCN3A", "SCN8A", "CACNA1S", "CACNA1C", "CACNA1D", "CACNA1B", "CACNA1E") %>% unique()

to_plot <- sort(unique(c(genes, bonus_genes)))
pdf(glue::glue("{outdir}/pc_featureplots.pdf"), width = 25, height = 5 * ceiling(length(to_plot) / 5))
FeaturePlot(pc, to_plot, raster = T, ncol = 5) & NoAxes() & NoLegend()
dev.off()

# add pt
pt_to_add <- ggdf %>% 
    tibble::column_to_rownames("barcode") %>% 
    dplyr::select(pc_pseudotime = Lineage1)
pc %<>% AddMetaData(pt_to_add)

# cut pt into bins
pc@meta.data %<>% dplyr::mutate(pt_bin = cut(pc_pseudotime, 50, labels = paste0("bin", 1:50)))

# get weighted mean pseudotime to assign gene order in heatmap
wmp <- apply(pc@assays$RNA@data[pc_genes,], 1, function(x){
    x <- minmax_scale(x) # so all weights are >= 0
    non0 <- which(x!=0)
    wex <- weighted.mean(pc$pc_pseudotime[non0], x[non0])
    return(wex)}) %>% sort()

# get mean expression by bin
xp_df <- pc@assays$RNA@data[pc_genes, ] %>% 
    as.matrix() %>% 
    t() %>% 
    tibble::as_tibble(rownames = "sample_barcode") %>% 
    dplyr::mutate(pt_bin = plyr::mapvalues(sample_barcode, rownames(pc@meta.data), pc@meta.data$pt_bin)) %>% 
    tidyr::pivot_longer(c(-pt_bin, -sample_barcode), names_to = "gene", values_to = "xp") %>% 
    dplyr::group_by(pt_bin, gene) %>% 
    dplyr::summarise(mean_xp = mean(xp)) %>% 
    dplyr::group_by(gene) %>% 
    dplyr::mutate(scaled_mean_xp = minmax_scale(mean_xp),
    gene = factor(gene, levels = rev(names(wmp)))) #


# density plot
dens <- ggplot(pc@meta.data, aes(x = pc_pseudotime, fill = celltype_factor), colour = "black")+
  geom_density(alpha = 0.8, adjust = 3)+
  scale_fill_manual(values = celltype_pal)+
  theme_void()+
  theme(legend.position = "none")+
  labs(x = "", y = "")

# heatmap
h <- ggplot(xp_df, aes(x = as.numeric(pt_bin), y = gene, fill = scaled_mean_xp))+
  geom_tile()+
  scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(9, "RdYlBu")), breaks = c(0,1))+
  theme_void()+
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(face = "italic", size = 2),
        legend.position = "bottom",
        legend.key.size = unit(5, "pt"))+
  labs(x = "", y = "", fill = "Scaled\nmean\nNorm. Counts")

pdf(glue::glue("{outdir}/pc_heatmap.pdf"), h = 4, w = 3)
dens / h + patchwork::plot_layout(heights = c(1,4))
dev.off()


# IN subclustering ---------------------------------------------------------------
in_genes <- sort(c(genes, "KCNA2", 'GRIK3', 'NXPH1', 'SLC6A5', 'ALDH1A3', 'LGI2', 'GJD2', 'SORCS3'))

DefaultAssay(sss) <- "integrated"
inn <- sss %>% 
    subset(subset = celltype %in% c("NE", "PIP 1", "PIP 2", "IN 1", "IN 2", "IN 3")) %>% 
    FindVariableFeatures(selection.method = "vst") %>% 
    ScaleData() %>% 
    RunPCA() %>% 
    FindNeighbors(dims = 1:20) %>%
    FindClusters(resolution = 0.3) %>% 
    RunUMAP(dims = 1:20)


# quick save
qs::qsave(inn, glue::glue("{outdir}/inn.qs"))
inn <- qs::qread(glue::glue("{outdir}/inn.qs"))

# palettes
lvls <- c("6 PCW", "7 PCW", "8 PCW", "9 PCW", "10 PCW", "11 PCW", "12 PCW", "13 PCW", "14 PCW", "15 PCW", "16 PCW", "17 PCW", "18 PCW", "20 PCW", "21 PCW", "38 PCW", "6 mths", "7 mths", "9 mths", "3 yrs", "29 yrs", "42 yrs", "44 yrs", "46 yrs", "50 yrs", "52 yrs")
celltype_lvls <- c("NE", "PC 1", "PC 2", "PC 3", "PC 4", "PIP 1", "PIP 2", "IN 1", "IN 2", "IN 3")
age_pal <- colorRampPalette(RColorBrewer::brewer.pal(11, 'Spectral'))(length(unique(lvls)))
celltype_pal <- tableau10 %>% magrittr::set_names(celltype_lvls)

# quick plots
inn_plot1 <- DimPlot(inn, group.by = "celltype_factor", reduction = "umap", raster = TRUE, cols = celltype_pal) & NoAxes() & ggtitle("")
inn_plot2 <- DimPlot(inn, group.by = "Phase", reduction = "umap", raster = TRUE, cols = phase_pal) & NoAxes() & ggtitle("")
inn_plot3 <- DimPlot(inn, group.by = "age_factor", reduction = "umap", raster = TRUE, cols = age_pal) & NoAxes() & ggtitle("")
pdf(glue::glue("{outdir}/interneuron_umaps.pdf"), width = 18, height = 5)
inn_plot1 | inn_plot2 | inn_plot3
dev.off()

DefaultAssay(inn) <- "RNA"
pdf(glue::glue("{outdir}/interneuron_featureplots.pdf"), width = 25, height = 95)
FeaturePlot(inn, in_genes, raster = T, ncol = 5) & NoAxes() & NoLegend()
dev.off()

in_select_feats <- c("KCNA2", 'GRIK3', 'NXPH1', 'SLC6A5', 'ALDH1A3', 'LGI2', 'GJD2', 'SORCS3')
pdf(glue::glue("{outdir}/interneuron_featureplots_select.pdf"), width = 10, height = 20)
FeaturePlot(inn, in_select_feats, raster = T, order = F, ncol = 2) & NoAxes() & NoLegend()
dev.off()


# slingshot
pto <- slingshot::slingshot(data = Embeddings(inn, reduction = "umap"),
    clusterLabels = inn$celltype, 
    start.clus = "NE",
    end.clus = c("IN 1", "PIP 2"))
sds <- slingshot::as.SlingshotDataSet(pto)

# see number of curves
slingshot::slingCurves(sds) %>% length()

# build df for ggplot
ggdf <- Embeddings(inn, reduction = "umap") %>%
  as.data.frame() %>% 
  dplyr::mutate(celltype = inn$celltype) %>% 
  cbind(pto@assays@data@listData$pseudotime) %>% 
  tibble::rownames_to_column("barcode")

# get curve coords
crv1 <- slingshot::slingCurves(sds)[[1]]$s %>% as.data.frame() 
crv2 <- slingshot::slingCurves(sds)[[2]]$s %>% as.data.frame()
crv3 <- slingshot::slingCurves(sds)[[3]]$s %>% as.data.frame()

# plot
p1 <- ggplot(ggdf, aes(UMAP_1, UMAP_2, color = celltype))+
        geom_point()+
        theme_void()+
        theme(legend.text=element_text(size=10),
              legend.title=element_blank())+
        xlab("UMAP 1")+
        ylab("UMAP 2")+
        scale_colour_manual(values = celltype_pal)+
        geom_path(data = crv1, aes(UMAP_1, UMAP_2), size = 0.75, color = "black")+
        geom_path(data = crv2, aes(UMAP_1, UMAP_2), size = 0.75, color = "black")+
        geom_path(data = crv3, aes(UMAP_1, UMAP_2), size = 0.75, color = "black")+
        guides(color = guide_legend(override.aes = list(size = 3)))

p2 <- ggplot(ggdf, aes(UMAP_1, UMAP_2, color = Lineage1))+
        geom_point()+
        theme_void()+
        theme(legend.text=element_text(size=10),
              legend.title=element_blank())+
        xlab("UMAP 1")+
        ylab("UMAP 2")+
        viridis::scale_color_viridis(option = "rocket")+
        guides(color = guide_legend(override.aes = list(size = 3)))

pdf(glue::glue("{outdir}/slingshot_inn.pdf"), width = 12, height = 5)
p1 | p2
dev.off()

