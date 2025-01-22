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


# Bonus for revisions 2024-10-28 -----------------------------------------

# import hallmark apoptosis gene sets
read_hallmark <- function(x){
    x %>% 
        readr::read_tsv() %>% 
        dplyr::filter(STANDARD_NAME == "GENE_SYMBOLS") %>% 
        dplyr::pull(HALLMARK_APOPTOSIS) %>% 
        stringr::str_split(",") %>% 
        unlist() %>% 
        return()
}

apop_hs <- read_hallmark("src/metadata/HALLMARK_APOPTOSIS.v2024.1.Hs.tsv")
apop_mm <- read_hallmark("src/metadata/HALLMARK_APOPTOSIS.v2024.1.Mm.tsv")

# re-import integrated human dataset
int <- qs::qread(glue::glue("out/integration/pc/final_object/integrated.qs"))

# apply module scores
int %<>% AddModuleScore(list(apop_hs), name = "hallmark_apoptosis")
mv %<>% AddModuleScore(list(apop_mm), name = "hallmark_apoptosis")
mm_lin %<>% AddModuleScore(list(apop_mm), name = "hallmark_apoptosis")
ald %<>% AddModuleScore(list(apop_hs), name = "hallmark_apoptosis")
br %<>% AddModuleScore(list(apop_hs), name = "hallmark_apoptosis")
cs20 %<>% AddModuleScore(list(apop_hs), name = "hallmark_apoptosis")
lake %<>% AddModuleScore(list(apop_hs), name = "hallmark_apoptosis")
sil %<>% AddModuleScore(list(apop_hs), name = "hallmark_apoptosis")

# clean int object
int$study_id <- int$sample_id %>% stringr::str_replace('ts', 'this study') %>% snakecase::to_sentence_case()

# plot on int object
int1 <- FeaturePlot(int, "hallmark_apoptosis1", raster = T) & NoAxes() & NoLegend() & ggtitle("Apoptosis Score (Hallmark)")
int2 <- DimPlot(int, group.by = "celltype_factor", raster = T, label = T, cols = tableau10) & NoAxes() & NoLegend() & ggtitle("")

pdf("sandbox/apop.pdf", w = 12, h = 6)
int1 | int2
dev.off()

# plot int violins
subset_vln <- function(y, z){
    celltype_lvls <- c("NE", "PC 1", "PC 2", "PC 3", "PC 4", "PIP 1", "PIP 2", "IN 1", "IN 2", "IN 3")
    ct_vln_pal <- tableau10 %>% magrittr::set_names(celltype_lvls)
    y %>% 
        subset(subset = study_id == z) %>% 
        VlnPlot("hallmark_apoptosis1", group.by = "celltype_factor", pt.size = 0, cols = ct_vln_pal) & 
        NoLegend() & 
        theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5)) & 
        ggtitle(z) & 
        ylab("Apoptosis Score (Hallmark)")
}

plot_list <- purrr::map(unique(int$study_id), ~subset_vln(int, .x))

pdf("sandbox/apop_vln.pdf", w = 6, h = length(unique(int$sample_id)) * 3)
cowplot::plot_grid(plotlist = plot_list, ncol = 1)
dev.off()


pdf("sandbox/apop_vln_mv.pdf", w = 15, h = 6)
VlnPlot(mv, "hallmark_apoptosis1", group.by = "Orig_ann", pt.size = 0) & 
        NoLegend() & 
        theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1)) & 
        ylab("Apoptosis Score (Hallmark)")&
        ggtitle("")
dev.off()



