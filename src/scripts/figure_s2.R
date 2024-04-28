# Info --------------------------------------------------------------------

# Anders E.
# Taylor lab
# April 5, 2024


# Notes -------------------------------------------------------------------

# Figure 3N: Replace 3N with dot plots and create a corresponding UMAP (featureplots) in Supplemental Figure 2A
# Supplemental Figure 2B: stacked bar plots stratifying cell types across x axis with stacked layers indicating the relative contributions of G1-S-G2-M phase per cells (from Figure 4).
# Supplemental Figure 2C and 2D: Timepoints at which PRDM13+ KI67+ and SKOR2+  KI67+ cells are identified.


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

# import dev mouse data
mv <- readRDS(paste0(yml$mv_dir, "/original_manuscript_mouse.rds"))
mm_lin <- qs::qread("out/linnarsson_hindbrain_mouse.qs")

# import dev human data
ald <- readRDS(glue::glue("{yml$ald_path}/cbl_integrated_cleanCC_210111.rds"))
br <- qs::qread(glue::glue("{yml$lin_dir}/hb_cb_final_v3.qs"))
cs20 <- qs::qread("out/cs20s_annot.qs")

# import adult human data
lake <- qs::qread("src/data/lake_cerebellum.qs")
sil <- qs::qread(glue::glue("{yml$public_scrnaseq_dir}/siletti_2022_biorxiv/sil_cerebellum.qs"))


# Panel A -------------------------------------------------------------
# featureplots to correspond to 3N

a_feats <- c("PVALB", "SST", "PAX2")
DefaultAssay(ald) <- DefaultAssay(br) <- DefaultAssay(lake) <- "RNA"
p1 <- FeaturePlot(ald, a_feats, raster = T, order = T, ncol = 1, pt.size = 2.5) & NoAxes() & NoLegend()
p2 <- FeaturePlot(br, a_feats, raster = T, order = T, ncol = 1, pt.size = 2.5) & NoAxes() & NoLegend()
p3 <- FeaturePlot(lake, a_feats, raster = T, order = T, ncol = 1, pt.size = 2.5) & NoAxes() & NoLegend()

pdf("out/feats_pvalb_sst_pax2.pdf", width = 9, height = 9)
p1 | p2 | p3
dev.off()


# Panel B -------------------------------------------------------------

gg <- cs20@meta.data %>% 
    dplyr::select(celltype_fine, Phase) %>% 
    dplyr::filter(grepl("^PC", celltype_fine)) %>% 
    ggplot(aes(celltype_fine, fill = Phase))+
    geom_bar(position="fill", stat="count")+
    theme_classic()+
    scale_fill_manual(values = phase_pal)+
    ylab("Proportion of cells")+
    xlab("")+
    theme(axis.text.x = element_text(angle = 45, hjust = 1))


pdf("out/stacked_bar_phase.pdf", width = 5, height = 4)
print(gg)
dev.off()


# Panel C/D -------------------------------------------------------------
# made in src/scripts/figure_4.R


