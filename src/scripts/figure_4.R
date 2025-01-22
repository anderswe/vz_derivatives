# Info --------------------------------------------------------------------

# Anders E.
# Taylor lab
# April 24, 2023


# Notes -------------------------------------------------------------------

# 
# original plots by Liam H.
# 


# Libraries ---------------------------------------------------------------
# module load R/4.2.1_arrow
# R
.libPaths(c("src/envs/r_vz_deriv/renv/library/R-4.2/x86_64-pc-linux-gnu", .libPaths()))
library(Seurat)
library(magrittr)
library(ggplot2)
library(patchwork)
library(yaml)
library(slingshot)
library(tradeSeq)
options(future.globals.maxSize = 1.2 * 1024 ^ 3)
yml <- yaml::read_yaml("config.yaml")
source("src/scripts/utils.R")


# Panels A-G prep ---------------------------------------------------------------------

# load
# so <- readRDS(glue::glue("{yml$lh_dev_dir}/seurat_objs/CS_10X_merge_HOXless_CS20s_SCTccmtreg_PC10_final_nomicroglia.RDS")) %>% 
#     FindNeighbors(dims = 1:20, reduction = "pca", verbose = T) %>%
#     RunUMAP(dims = 1:20, reduction = "pca", min.dist = 0.1, spread = 1.2) %>% 
#     FindClusters(resolution = 2, verbose = T) # res = 1
# qs::qsave(so, "out/cs20s.qs")

# pdf('sandbox/feats_foxp2_cs20.pdf', h = 5, w = 25)
# FeaturePlot(cs20, c("FOXP2", "CALB2", "SKOR2", "PRDM13", "MKI67"), raster = T, ncol = 5) & NoAxes() & NoLegend()
# dev.off()

so <- qs::qread("out/cs20s.qs")

# convert to h5ad for scFates
# Embeddings(so, reduction = "umap") %>% 
#   as.data.frame() %>% 
#   cbind(so@meta.data) %>% 
#   tibble::rownames_to_column("sample_barcode") %>% 
#   data.table::fwrite("out/metadata_cs20s.csv")
# SeuratDisk::SaveH5Seurat(so, filename = "out/cs20s.h5Seurat")
# SeuratDisk::Convert("out/cs20s.h5Seurat", dest = "h5ad")


# dimplot
rast <- F
dp1 <- DimPlot(so, group.by = c("Cell_type"), raster = rast, cols = tableau20) & NoAxes() & ggtitle("")
dp2 <- DimPlot(so, group.by = c("seurat_clusters"), raster = rast, label = T, cols = colorRampPalette(tableau20)(length(unique(so$seurat_clusters)))) & NoAxes() & ggtitle("")
dp3 <- DimPlot(so, group.by = "orig.ident", raster = rast, cols = RColorBrewer::brewer.pal(9, "Spectral")[c(4,7,9)])+NoAxes()+ggtitle("")+theme(legend.position = c(0.7, 0.2))
dp4 <- DimPlot(so, group.by = "Phase", raster = rast, cols = phase_pal) & NoAxes() & ggtitle("") & theme(legend.position = c(0.7, 0.2))

pdf("out/dimplot_cs20s.pdf", h = 20, w = 6)
(dp1 / dp2 / dp3 / dp4) + plot_layout(ncol = 1)
dev.off()


fp_genes <- genes %>% .[. %in% rownames(so)] %>% c(., "FOXP2", "FOXP1", "SORCS3", "NXPH1", "NES", "OPCML") %>% sort()
pdf("out/featureplot_cs20s.pdf", h = 110, w = 25)
FeaturePlot(so, fp_genes, order = F, raster = T, ncol = 5) & NoAxes() & NoLegend()
dev.off()


# panel A (annotation) --------------------------------------------------

so@meta.data %<>% dplyr::mutate(celltype_fine = dplyr::case_when(
  seurat_clusters %in% c(32,15) ~ "NE", # NFIX+, PLP1+, SALL3/4+
  seurat_clusters %in% c(7,14) ~ "VZ (Prolif.)", # PTF1A+, PIEZO2+, and S/G2M
  seurat_clusters %in% c(30,19,13,17,22) ~ "PC (PRDM13+)", # PRDM13+
  seurat_clusters %in% c(16,23,5) ~ "PC (SKOR2+)", # SKOR2+
  seurat_clusters %in% c(4,3,12,33,0,9) ~ "PC (ESRRB+)", # ESRRB+
  seurat_clusters %in% c(2,6) ~ "PC (PCP4+)", # PCP4+ 25???
  seurat_clusters %in% c(25) ~ "PC (GRIN2A+)", # GRIN2A+
  seurat_clusters %in% c(8,24,27,29) ~ "GABA CN (SOX14+)", # SOX14+
  seurat_clusters %in% c(26) ~ "RL (LMX1A+)", 
  seurat_clusters %in% c(21) ~ "RL deriv. (ATOH1+)",
  seurat_clusters %in% c(20,11,1,10,28,18,31) ~ "Gluta. CN", # LHX9+, PAX5+
) %>% factor(levels = rev(c("NE", "VZ (Prolif.)", "PC (PRDM13+)", "PC (SKOR2+)", "PC (ESRRB+)", "PC (PCP4+)", "PC (GRIN2A+)", "GABA CN (SOX14+)", "RL (LMX1A+)", "RL deriv. (ATOH1+)", "Gluta. CN"))))

set.seed(667)
ct_pal_new <- tableau20[c(11,2,8,3,14,12,10,19,1,17,5)] %>% magrittr::set_names(levels(so$celltype_fine))

dp_annot <- DimPlot(so, group.by = "celltype_fine", cols = ct_pal_new) + NoAxes() + ggtitle("") + NoLegend() #label = T, repel = T, label.size = 3
pdf("out/annotation_dimplots_cs20s_recolor.pdf", h = 4, w = 10)
dp_annot | dp3 | dp4
dev.off()

pdf("out/annotation_dimplots_cs20s_recolor_fine_annotation_only.pdf", h = 4, w = 4)
dp_annot
dev.off()

# PC subset
ss <- so %>% subset(subset = celltype_fine %in% c("PC (SKOR2+)", "PC (ESRRB+)", "PC (PCP4+)", "PC (GRIN2A+)")) %>% 
  SCTransform(vars.to.regress = "percent.mt") %>% 
  RunPCA() %>% 
  harmony::RunHarmony(group.by.vars = "orig.ident",
             reduction = "pca",
             assay.use = "RNA",
             theta = 2) %>%
  FindNeighbors(dims = 1:20, reduction = "harmony", verbose = T) %>%
  RunUMAP(dims = 1:20, reduction = "harmony") %>% 
  FindClusters(resolution = 0.2, verbose = T)

pdf("out/featureplot_cs20s_pc_subset.pdf", h = 110, w = 25)
FeaturePlot(ss, fp_genes, order = F, raster = T, ncol = 5) & NoAxes() & NoLegend()
dev.off()

dp1 <- DimPlot(ss, group.by = "celltype_fine", label = T, repel = T, raster = T) & NoAxes()
dp2 <- DimPlot(ss, group.by = "seurat_clusters", label = T, repel = T, raster = T) & NoAxes()
dp3 <- DimPlot(ss, group.by = "orig.ident", label = T, repel = T, raster = T) & NoAxes()
pdf("out/dimplot_cs20s_pc_subset.pdf", h = 5, w = 15)
dp1 | dp2 | dp3
dev.off()

mks_5 <- FindMarkers(ss, ident.1 = "5", only.pos = TRUE)
mks_4 <- FindMarkers(ss, ident.1 = "4", only.pos = TRUE)

head(mks_5, 20)
pdf("out/cs20s_markers_4.pdf", w = 20, h = 25)
FeaturePlot(so, features = rownames(mks_4)[1:20], raster = T) & NoAxes() & NoLegend()
dev.off()
pdf("out/cs20s_markers_5.pdf", w = 20, h = 25)
FeaturePlot(so, features = rownames(mks_5)[1:20], raster = T) & NoAxes() & NoLegend()
dev.off()
pdf("out/cs20s_pc_subset_markers_4.pdf", w = 20, h = 25)
FeaturePlot(ss, features = rownames(mks_4)[1:20], raster = T) & NoAxes() & NoLegend()
dev.off()
pdf("out/cs20s_pc_subset_markers_5.pdf", w = 20, h = 25)
FeaturePlot(ss, features = rownames(mks_5)[1:20], raster = T) & NoAxes() & NoLegend()
dev.off()

# recluster with harmony
hmny <- so %>% 
  harmony::RunHarmony(group.by.vars = "orig.ident",
             reduction = "pca",
             assay.use = "RNA",
             theta = 2) %>%
  FindNeighbors(dims = 1:20, reduction = "harmony", verbose = T) %>%
  RunUMAP(dims = 1:20, reduction = "harmony", spread = 1.2, min.dist = 0.2) %>% 
  FindClusters(resolution = 0.5, verbose = T)

pdf("out/dim_cs20s_hmny.pdf", h = 5, w = 15)
DimPlot(hmny, group.by = c("seurat_clusters", "celltype_fine", "orig.ident"), raster = T, ncol = 3) & NoAxes() & NoLegend()
dev.off()

pdf("out/featureplot_cs20s_hmny.pdf", h = 110, w = 25)
FeaturePlot(hmny, fp_genes, order = F, raster = T, ncol = 5) & NoAxes() & NoLegend()
dev.off()

# markers
# future::plan("sequential")
# Idents(so) <- so$celltype_fine
# mks <- FindAllMarkers(so, only.pos = T) %>% dplyr::mutate(unique = pct.1 - pct.2)
# qs::qsave(mks, "out/mks_cs20s.qs")
mks <- qs::qread("out/mks_cs20s.qs")

so@meta.data %<>% dplyr::mutate(celltype_broad = dplyr::case_when(
  grepl("LMX1A", celltype_fine) ~ "RL (LMX1A+)",
  grepl("ATOH", celltype_fine) ~ "RL deriv. (ATOH1+)",
  grepl("Glut", celltype_fine) ~ "Gluta. CN",
  grepl("GABA", celltype_fine) ~ "GABA CN",
  grepl("ESRRB|PCP|GRIN", celltype_fine) ~ "Late PC",
  grepl("SKOR2|PRDM13", celltype_fine) ~ "Early PC",
  grepl("VZ", celltype_fine) ~ "VZ",
  grepl("NE", celltype_fine) ~ "NE",
))

qs::qsave(so, "out/cs20s_annot.qs")
so <- qs::qread("out/cs20s_annot.qs")

ct_broad_pal <- c(tableau10[c(3,2,4,1,5,7)], tableau20[2], tableau10[6]) %>% magrittr::set_names(unique(so$celltype_broad))
dp_annot_broad <- DimPlot(so, group.by = "celltype_broad", cols = ct_broad_pal) + NoAxes() + ggtitle("") + NoLegend() #label = T, repel = T, label.size = 3
dp_annot_broad2 <- DimPlot(so, group.by = "celltype_broad", cols = ct_broad_pal, label = T, repel = T, label.size = 3) + NoAxes() + ggtitle("") + NoLegend()
bp_annot_broad <- so@meta.data %>% 
  dplyr::select(celltype_broad, orig.ident) %>% 
  table() %>% 
  as.data.frame() %>% 
  dplyr::rename(ct = 1, age = 2, freq = 3) %>% 
  ggplot(aes(fill=ct, y=freq, x=age)) + 
    geom_bar(position="fill", stat="identity", color = "black")+
    theme_classic()+
    scale_fill_manual(values = ct_broad_pal)+
    labs(y = "Proportion", x = "", fill = "")
pdf("out/annotation_dimplots_broad_cs20s.pdf", h = 8, w = 8)
(dp_annot_broad | dp_annot_broad2) / (bp_annot_broad + patchwork::plot_spacer() + patchwork::plot_layout(widths = c(1,2)))
dev.off()

# quick features
pdf("out/lmx1a_cs20s.pdf", h = 4, w = 4)
FeaturePlot(so, "LMX1A") & NoAxes()
dev.off()


# NE markers
Idents(so) <- so$celltype_fine
nemks <- FindMarkers(so, ident.1 = "NE") %>% dplyr::mutate(gene = rownames(.))
qs::qsave(nemks, glue::glue("out/ne_markers.qs"))
nemks <- qs::qread("out/ne_markers.qs")

so@meta.data %<>% dplyr::mutate(celltype_fine2 = dplyr::case_when(
    celltype_fine == "NE" ~ "NE",
    celltype_fine %in% c("VZ (Prolif.)", "RL") ~ "VZ & RL",
    TRUE ~ "Other"
))

Idents(so) <- so$celltype_fine2
nemks3 <- FindMarkers(so, ident.1 = "NE", ident.2 = "VZ & RL") %>% dplyr::mutate(gene = rownames(.))


do_volcano <- function(markers, lfc_thresh = 0.5, padj_thresh = 1e-75, label1 = "label1", label2 = "label2", goi = NULL){
  markers %>% 
    dplyr::filter(p_val_adj != 1 & !grepl("^ENSG", gene)) %>% 
    dplyr::mutate(color = dplyr::case_when(
      p_val_adj > padj_thresh | abs(avg_log2FC) < lfc_thresh ~ "insig",
      avg_log2FC > lfc_thresh ~ label1,
      avg_log2FC < -lfc_thresh ~ label2)) %>% 
    { if(!is.null(goi)) dplyr::mutate(., label = ifelse(gene %in% goi, gene, NA)) else dplyr::mutate(., label = ifelse(color != "insig", gene, NA)) } %>% 
    dplyr::mutate(neg_log_padj = -log10(p_val_adj) %>% ifelse(. == Inf, 305, .)) %>% # assign Inf to 305
    ggplot(aes(avg_log2FC, neg_log_padj, color = color, label = label))+
    geom_point(size = 3)+
    theme_classic()+
    scale_color_manual(breaks = c(label1, label2), values = c(fave_pal[c(3,1)], "lightgrey"), na.value="lightgrey", name = "")+
    ggrepel::geom_text_repel(size = 3, segment.size = 0.1, 
                             box.padding = 1, color = "black", 
                             segment.color = "black", min.segment.length = 0, 
                             max.overlaps = Inf)+
    labs(x = bquote(~log[2]~FC), y = bquote(~-log[10]~p[adj]), color = "")+
    guides(colour = guide_legend(override.aes = list(size=3)))
}

pdf("out/volcano_ne_vs_rl&vz_condensed.pdf", w = 4, h = 3)
do_volcano(nemks3, goi = c("PTPRZ1", "ZIC1", "ERBB4", "ATXN1", "WNT7B", "DCC", "NFIA", "LRP2"), label1 = "NE", label2 = "VZ & RL")
dev.off()

pdf("out/volcano_ne_vs_rl&vz_expanded.pdf", w = 10, h = 9)
do_volcano(nemks3, label1 = "NE", label2 = "VZ & RL")
dev.off()

pdf("out/volcano_ne_vs_all_condensed.pdf", w = 4, h = 3)
do_volcano(nemks, goi = c("WNT7B", "SHROOM", "PTPRZ1", "SMOC1", "LRP2", "CTNNA2", "EBF1", "EBF3", "LHX1", "LRP2"), label1 = "NE", label2 = "Other")
dev.off()

pdf("out/volcano_ne_vs_all_expanded.pdf", w = 20, h = 18)
do_volcano(nemks, label1 = "NE", label2 = "Other")
dev.off()




# ISH marker plots
pdf("out/ne_mks.pdf", w = 10, h = 4)
FeaturePlot(so, c("LRP2", "SMOC2", "WNT7B"), raster = T, order = T, ncol = 3) & NoAxes() & NoLegend()
dev.off()





# panel B (dotplot) --------------------------------------------------

dp_genes <- c("NOTCH1", "NES", "WLS", "TCF7L1", "PARD3", "SOX2", "NFIA", "NFIX", "PLEKHA7",
  "PTF1A", "MKI67", "PAX6", "PRDM13", "GSX1", "OLIG2", "LHX5", "PRKG2", "LHX1", "SKOR2", "ESRRB",
  "PCP4", "NTN4", "GRIN2A", "SOX14", "PEX5L", "PAX2", "LMX1A", "BARHL1", "ATOH1", "MEIS2", "LHX9", "LHX2", "PAX5")

dpal <- RColorBrewer::brewer.pal(9, "RdYlBu")[c(1,5,9)]

dotp <- DotPlot(so, features = dp_genes, assay = "SCT", group.by = "celltype_fine", cols = "RdYlBu") +
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
  guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white")))+
  labs(x="", y="")+
  theme(axis.text.x = element_text(angle=45, vjust = 1, hjust = 1, face = "italic"))

pdf("out/dotplot.pdf", h = 5, w = 11)
dotp
dev.off()


# panel C (featureplots) --------------------------------------------------
feats <- c("PTF1A", "PRDM13", "LHX1", "SKOR2", "ESRRB", "PCP4", "GRIN2A", "PARD3B", "OPCML")
DefaultAssay(so) <- "SCT"

# happen to have same scale
pdf(glue::glue("out/featureplots_{paste(feats, collapse = '_')}_cs20s.pdf"), h = 9, w = 9)
FeaturePlot(so, feats, cols = c("cadetblue1", "deeppink3"), order = F, raster = F, ncol = 3, pt.size = 0.4) & NoAxes() & theme(plot.title = element_text(size = 12, face = "bold.italic")) & labs(color = "Expr.\n(SCT)")
dev.off()



# panel E (pseudotime) --------------------------------------------------

# import pseudotime from scFates
pt <- data.table::fread("out/metadata_scfates.csv") %>% 
  dplyr::select(sample_barcode = V1, pt = t) %>% 
  tibble::column_to_rownames("sample_barcode")

# import pt-associated genes from scFates
# gaba_genes <- data.table::fread("out/scfates_de_branches.csv") %>% 
#   dplyr::rename(gene = V1) %>% 
#   dplyr::filter(branch == 8) %>% # select gaba
#   dplyr::slice_max(up_A, n = 80) %>% 
#   dplyr::pull(gene)

# add md
so %<>% AddMetaData(pt)

# quick plot
pdf("out/dimplot_pt.pdf", h = 5, w = 5)
FeaturePlot(so, "pt", cols = viridis::viridis(100))+NoAxes()+ggtitle("")
dev.off()

# subset to gaba lineage minus gaba ntz
gaba_cells <- so@meta.data %>% dplyr::filter(grepl("NE|PC", celltype_broad)) %>% rownames()
ss <- subset(so, cells = gaba_cells)

# subset to corresponding marker genes
gaba_genes <- mks %>% 
  dplyr::filter(grepl("NE|PC|VZ", cluster)) %>% 
  dplyr::group_by(cluster) %>% 
  dplyr::filter(pct.2 < 0.2) %>% 
  dplyr::slice_min(p_val_adj, n = 20) %>% 
  dplyr::slice_max(avg_log2FC, n = 20) %>% 
  dplyr::arrange(factor(cluster, levels = c("NE", "VZ (Prolif.)", "PC (PRDM13+)", "PC (SKOR2+)", "PC (ESRRB+)", "PC (GRIN2A+)", "PC (PCP4+)"))) %>% 
  dplyr::pull(gene) %>% 
  unique() %>% 
  R.utils::insert(34, "PTF1A") %>% 
  rev()

# cut pt into bins
ss@meta.data %<>% dplyr::mutate(pt_bin = cut(pt, 50, labels = paste0("bin", 1:50)))

# get weighted mean pseudotime to assign gene order in heatmap
wmp <- apply(ss@assays$SCT@data[gaba_genes,], 1, function(x){
    x <- minmax_scale(x) # so all weights are >= 0
    non0 <- which(x!=0)
    wex <- weighted.mean(ss$pt[non0], x[non0])
    return(wex)}) %>% sort()

# get mean expression by bin
xp_df <- ss@assays$SCT[gaba_genes, ] %>%
  as.matrix() %>% 
  t() %>% 
  tibble::as_tibble(rownames = "sample_barcode") %>% 
  dplyr::mutate(pt_bin = plyr::mapvalues(sample_barcode, rownames(ss@meta.data), ss@meta.data$pt_bin)) %>% 
  tidyr::pivot_longer(c(-pt_bin, -sample_barcode), names_to = "gene", values_to = "sct") %>% 
  dplyr::group_by(pt_bin, gene) %>% 
  dplyr::summarise(mean_sct = mean(sct)) %>% 
  dplyr::group_by(gene) %>% 
  dplyr::mutate(scaled_mean_sct = minmax_scale(mean_sct),
    gene = factor(gene, levels = gaba_genes)) #rev(names(wmp))


# density plot
dens <- ggplot(ss@meta.data, aes(x = pt, fill = celltype_fine), colour = "black")+
  geom_density(alpha = 0.8, adjust = 3)+
  scale_fill_manual(values = ct_pal_new)+
  theme_void()+
  theme(legend.position = "none")+
  labs(x = "", y = "")

# heatmap
h <- ggplot(xp_df, aes(x = as.numeric(pt_bin), y = gene, fill = scaled_mean_sct))+
  geom_tile()+
  scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(9, "RdYlBu")), breaks = c(0,1))+
  theme_void()+
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(face = "italic", size = 2),
        legend.position = "bottom",
        legend.key.size = unit(5, "pt"))+
  labs(x = "", y = "", fill = "Scaled\nmean\nSCT")

pdf("out/heatmap_recolor.pdf", h = 4, w = 3)
dens / h + patchwork::plot_layout(heights = c(1,4))
dev.off()

# dimplot for legend
so@meta.data %<>% dplyr::mutate(ct_fine = factor(celltype_fine, levels = rev(levels(so$celltype_fine))))
dp_annot_legend <- DimPlot(so, group.by = "ct_fine", cols = ct_pal) + NoAxes() + ggtitle("") #label = T, repel = T, label.size = 3
pdf("out/annotation_legend_cs20s.pdf", h = 4, w = 13)
dp_annot_legend | dp3 | dp4
dev.off()

# tradeseq ###
# tso <- ss %>% subset(features = gaba_genes)
# tso_cts <- tso %>% GetAssayData("SCT", slot = "counts") %>% as.matrix()
# tso_pt <- tso$pt
# tso_cw <- rep(1, length(tso$pt))

# # fit gam
# sce <- fitGAM(counts = tso_cts, pseudotime = tso_pt, cellWeights = tso_cw, nknots = 6, verbose = T)

# # smooth
# ydf <- tradeSeq::predictSmooth(sce, gene = gaba_genes, nPoints = 50, tidy = T) %>% 
#   tibble::as_tibble() %>% 
#   dplyr::group_by(gene) %>% 
#   dplyr::mutate(yhat_scaled = minmax_scale(yhat),
#     gene = factor(gene, levels = gaba_genes))

# # heatmap
# y <- ggplot(ydf, aes(x = time, y = gene, fill = yhat_scaled))+
#   geom_tile()+
#   scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(9, "RdYlBu")))+
#   theme_void()+
#   theme(axis.text.x = element_blank(),
#         axis.text.y = element_text(face = "italic", size = 4))+
#   labs(x = "", y = "", fill = "Smoothed\nExpr.")

# pdf("out/heatmap_tradeseq.pdf", h = 6, w = 5)
# y
# dev.off()
###

# Bonus (mouse - Vladoiu) -------------------------------------------------------------------

mv_dir <- paste0(yml$mv_dir, "/original_manuscript_mouse.rds")
mv <- readRDS(mv_dir)

mv@meta.data %<>% dplyr::mutate(age = factor(Donor, levels = c(glue::glue("E{c(10,12,14,16,18)}"), glue::glue("P{c(0,5,7,14)}"))))
mv1 <- DimPlot(mv, group.by = "age", raster = T, cols = colorRampPalette(RColorBrewer::brewer.pal(9, "Spectral"))(length(unique(mv$age)))) + NoAxes() + ggtitle("")
mv2 <- DimPlot(mv, group.by = "Orig_ann", raster = T, label = T, cols = colorRampPalette(tableau20)(length(unique(mv$Orig_ann)))) + NoAxes() + ggtitle("")
pdf("out/mv.pdf", h = 20, w = 50)
mv1 | mv2
dev.off()

# convert to mouse
orth <- data.table::fread("src/metadata/new_ortho.csv")
hox_mm <- orth %>% dplyr::filter(grepl("Hox", mouse_gene_name)) %>% dplyr::pull(mouse_gene_name)
mm_genes <- plyr::mapvalues(genes, orth$gene_name, orth$mouse_gene_name, warn_missing = FALSE) %>% c(., hox_mm) %>% unique() %>% sort() %>% .[. %in% rownames(mv)]

pdf("out/featureplots_maria_mm.pdf", h = 80, w = 20)
FeaturePlot(mv, mm_genes, order = T, raster = T) & NoAxes() & NoLegend()
dev.off()

# featurescatter loop
fs_out <- "out/featurescatters"
dir.create(fs_out, showWarnings = F, recursive = T)
DefaultAssay(mv) <- "RNA"
iters <- c("Pax2", "Skor2", "Mki67", "Sox2") %>% combn(m = 2) %>% as.data.frame() %>% as.list()
for(i in iters){
  
  pdf(glue::glue("{fs_out}/mm_maria_{paste(i, collapse = '_')}.pdf"), h = 5, w = 5)
  print(FeatureScatter(mv, i[1], i[2], group.by = "orig.ident", slot = "counts", raster = T, pt.size = 0) + NoLegend() + ggtitle("") + geom_jitter())
  dev.off()

}


# Bonus (mouse - Linnarsson) -------------------------------------------------------------------


# mouse_dir <- paste0(dirname(yml$lin_dir), "/linnarsson_pb")
# mm <- readRDS(glue::glue("{mouse_dir}/linnarsson_pb_20210823.rds"))

# # initial dimplot
# mm@meta.data %<>% dplyr::mutate(age_day = Age %>% stringr::str_replace("e", "") %>% as.numeric())
# pdf("out/dimplot_mm.pdf", h = 5, w = 5)
# DimPlot(mm, group.by = "age_day", cols = colorRampPalette(RColorBrewer::brewer.pal(9, "Spectral"))(length(unique(mm$Age))), raster = T) + NoAxes() + ggtitle("") + labs(color = "Age\n(Emb. Day)")
# dev.off()

# # subset
# dims <- 20
# res <- 0.2
# ssmm <- mm %>% 
#   subset(subset = Tissue == "Hindbrain") %>% 
#   SCTransform(vars.to.regress = "percent.mt") %>% 
#   RunPCA() %>% 
#   FindNeighbors(dims = 1:dims) %>%
#   FindClusters(resolution = res) %>% 
#   RunUMAP(dims = 1:dims)
# qs::qsave(ssmm, "out/linnarsson_hindbrain_mouse.qs")
ssmm <- qs::qread("out/linnarsson_hindbrain_mouse.qs")

m1 <- DimPlot(ssmm, group.by = "age_day", cols = colorRampPalette(RColorBrewer::brewer.pal(9, "Spectral"))(length(unique(ssmm$Age))), raster = T) + NoAxes() + ggtitle("") + labs(color = "Age")
m2 <- DimPlot(ssmm, group.by = "Phase", cols = phase_pal, raster = T) + NoAxes() + ggtitle("") + labs(color = "Phase")
m3 <- DimPlot(ssmm, group.by = "seurat_clusters", cols = tableau20, raster = T) + NoAxes() + ggtitle("") + labs(color = "Cluster")
m4 <- DimPlot(ssmm, group.by = "Subclass", raster = T, cols = colorRampPalette(tableau20)(length(unique(ssmm$Subclass)))) + NoAxes() + ggtitle("") + labs(color = "Subclass")
m5 <- DimPlot(ssmm, group.by = "Subclass", raster = T, cols = colorRampPalette(tableau20)(length(unique(ssmm$Subclass))), label = T) + NoAxes() + ggtitle("") + NoLegend()
m6 <- DimPlot(ssmm, group.by = "Class", raster = T, cols = tableau20, label = T) + NoAxes() + ggtitle("") + NoLegend()

pdf("out/dimplot_mm_hb.pdf", h = 5, w = 15)
m1 | m2 | m3
dev.off()

pdf("out/subclass_mm_hb.pdf", h = 5, w = 20)
m4
dev.off()

pdf("out/subclass_mm_hb_labelled.pdf", h = 20, w = 40)
m5 | m6
dev.off()


# convert to mouse
orth <- data.table::fread("src/metadata/new_ortho.csv")
hox_mm <- orth %>% dplyr::filter(grepl("Hox", mouse_gene_name)) %>% dplyr::pull(mouse_gene_name)
mm_genes <- plyr::mapvalues(genes, orth$gene_name, orth$mouse_gene_name, warn_missing = FALSE) %>% c(., hox_mm) %>% unique() %>% sort() %>% .[. %in% rownames(ssmm)]


pdf("out/featureplots_linnarsson_mm.pdf", h = 80, w = 20)
FeaturePlot(ssmm, mm_genes, order = T, raster = T) & NoAxes() & NoLegend()
dev.off()

s_mm <- Seurat::cc.genes.updated.2019$s.genes %>% plyr::mapvalues(orth$gene_name, orth$mouse_gene_name, warn_missing = FALSE)
g2m_mm <- Seurat::cc.genes.updated.2019$g2m.genes %>% plyr::mapvalues(orth$gene_name, orth$mouse_gene_name, warn_missing = FALSE)
ssmm %<>% Seurat::CellCycleScoring(s.features = s_mm, g2m.features = g2m_mm, set.ident = F)

# featurescatter loop
fs_out <- "out/featurescatters"
DefaultAssay(ssmm) <- "RNA"
iters <- c("Pax2", "Skor2", "Mki67", "Sox2") %>% combn(m = 2) %>% as.data.frame() %>% as.list()
for(i in iters){
  
  pdf(glue::glue("{fs_out}/mm_lin_{paste(i, collapse = '_')}.pdf"), h = 5, w = 5)
  print(FeatureScatter(ssmm, i[1], i[2], group.by = "orig.ident", slot = "counts", raster = T, pt.size = 0) + NoLegend() + ggtitle("") + geom_jitter())
  dev.off()

}





# FeatureScatters ----------------------------

# outdir
fs_out <- "out/featurescatters"

# import dev mouse data
mv_dir <- paste0(yml$mv_dir, "/original_manuscript_mouse.rds")
mv <- readRDS(mv_dir)
mm_lin <- qs::qread("out/linnarsson_hindbrain_mouse.qs")

# import dev human data
so <- readRDS(glue::glue("{yml$ald_path}/cbl_integrated_cleanCC_210111.rds"))
lin <- qs::qread(glue::glue("{yml$lin_dir}/hb_cb_final_v3.qs"))
cs20 <- qs::qread("out/cs20s.qs")

# import adult human data
lake <- qs::qread("src/data/lake_cerebellum.qs")
sil <- qs::qread(glue::glue("{yml$public_scrnaseq_dir}/siletti_2022_biorxiv/sil_cerebellum.qs"))


# feature selection
feats <- c("MKI67", "SKOR2") #c("MKI67", "PRDM13") #c("SKOR2", "MKI67")
mmfeats <- c("Mki67", "Skor2") #c("Mki67", "Prdm13")

# assign groups
bc_list <- purrr::map(c(lake, sil), \(x){ #c(so, lin, cs20, mv, mm_lin, lake, sil)

  if("Prdm13" %in% rownames(x@assays$RNA@counts)){ feats_loop <- mmfeats } else { feats_loop <- feats }

  x@assays$RNA@counts[feats_loop, ] %>% 
    as.matrix() %>% 
    t() %>% 
    as.data.frame() %>% 
    dplyr::filter(dplyr::if_all(dplyr::everything(), ~ . != 0)) %>% 
    rownames() %>% 
    return()

}) %>% magrittr::set_names(c("lake", "sil")) #c("so", "lin", "cs20", "mv", "mm_lin", "lake", "sil")


so@meta.data %<>% dplyr::mutate(feat_status = ifelse(rownames(.) %in% bc_list[["so"]], "double_pos", "not"))
lin@meta.data %<>% dplyr::mutate(feat_status = ifelse(rownames(.) %in% bc_list[["lin"]], "double_pos", "not"))
cs20@meta.data %<>% dplyr::mutate(feat_status = ifelse(rownames(.) %in% bc_list[["cs20"]], "double_pos", "not"))
mv@meta.data %<>% dplyr::mutate(feat_status = ifelse(rownames(.) %in% bc_list[["mv"]], "double_pos", "not"))
mm_lin@meta.data %<>% dplyr::mutate(feat_status = ifelse(rownames(.) %in% bc_list[["mm_lin"]], "double_pos", "not"))
lake@meta.data %<>% dplyr::mutate(feat_status = ifelse(rownames(.) %in% bc_list[["lake"]], "double_pos", "not"))
sil@meta.data %<>% dplyr::mutate(feat_status = ifelse(rownames(.) %in% bc_list[["sil"]], "double_pos", "not"))

# stacked barplots by timepoint
sbp1 <- ggplot(so@meta.data %>% dplyr::filter(feat_status == "double_pos"), aes(fill = feat_status, x = age)) + geom_bar(position="stack", stat="count") + scale_fill_manual(values = "deeppink3") + theme_classic() + theme(legend.position = "none") + ylab("Num. cells") + xlab("") + ggtitle("Aldinger et al. 2021")
sbp2 <- ggplot(lin@meta.data %>% dplyr::filter(feat_status == "double_pos"), aes(fill = feat_status, x = Age)) + geom_bar(position="stack", stat="count") + scale_fill_manual(values = "deeppink3") + theme_classic() + theme(legend.position = "none") + ylab("Num. cells") + xlab("") + ggtitle("Braun et al. 2023")
sbp3 <- ggplot(cs20@meta.data %>% dplyr::filter(feat_status == "double_pos"), aes(fill = feat_status, x = orig.ident)) + geom_bar(position="stack", stat="count") + scale_fill_manual(values = "deeppink3") + theme_classic() + theme(legend.position = "none") + ylab("Num. cells") + xlab("") + ggtitle("This study")
sbp4 <- ggplot(mv@meta.data %>% dplyr::filter(feat_status == "double_pos"), aes(fill = feat_status, x = Donor)) + geom_bar(position="stack", stat="count") + scale_fill_manual(values = "deeppink3") + theme_classic() + theme(legend.position = "none") + ylab("Num. cells") + xlab("") + ggtitle("Vladoiu et al. 2019")
sbp5 <- ggplot(mm_lin@meta.data %>% dplyr::filter(feat_status == "double_pos"), aes(fill = feat_status, x = age_day)) + geom_bar(position="stack", stat="count") + scale_fill_manual(values = "deeppink3") + theme_classic() + theme(legend.position = "none") + ylab("Num. cells") + xlab("") + ggtitle("La Manno et al. 2021")
sbp6 <- ggplot(lake@meta.data %>% dplyr::filter(feat_status == "double_pos"), aes(fill = feat_status, x = age)) + geom_bar(position="stack", stat="count") + scale_fill_manual(values = "deeppink3") + theme_classic() + theme(legend.position = "none") + ylab("Num. cells") + xlab("") + ggtitle("Lake et al. 2017")
sbp7 <- ggplot(sil@meta.data %>% dplyr::filter(feat_status == "double_pos"), aes(fill = feat_status, x = Age)) + geom_bar(position="stack", stat="count") + scale_fill_manual(values = "deeppink3") + theme_classic() + theme(legend.position = "none") + ylab("Num. cells") + xlab("") + ggtitle("Siletti et al. 2023")

pdf(glue::glue("{fs_out}/stacked_barplots_by_age_{paste(feats, collapse = '_')}.pdf"), h = 15, w = 7)
sbp3 / sbp1 / sbp2 / sbp4 / sbp5
dev.off()

pdf(glue::glue("{fs_out}/stacked_barplots_by_age_lake_siletti_{paste(feats, collapse = '_')}.pdf"), h = 6, w = 7)
sbp6 / sbp7
dev.off()

# adult featureplots
pdf(glue::glue("sandbox/lake_skor2_mki67.pdf"), w = 10, h = 5)
FeaturePlot(lake, c("SKOR2", "MKI67"), order = F, raster = T) & NoAxes()
dev.off()
pdf(glue::glue("sandbox/siletti_skor2_mki67.pdf"), w = 10, h = 5)
FeaturePlot(sil, c("SKOR2", "MKI67"), order = F, raster = T) & NoAxes()
dev.off()


# feature scatterplots
jit <- TRUE
pt_size <- 0.5
grp <- "feat_status"
pal <- c("double_pos" = "deeppink3", "not" = "cadetblue1")
DefaultAssay(so) <- DefaultAssay(lin) <- DefaultAssay(cs20) <- DefaultAssay(mv) <- DefaultAssay(mm_lin) <- "RNA"
fs1 <- FeatureScatter(lin, feats[1], feats[2], group.by = grp, raster = T, slot = "counts", pt.size = pt_size, jitter = jit, cols = pal) + NoLegend() + ggtitle("Braun et al. 2022")
fs2 <- FeatureScatter(so, feats[1], feats[2], group.by = grp, raster = T, slot = "counts", pt.size = pt_size, jitter = jit, cols = pal) + NoLegend() + ggtitle("Aldinger et al. 2021")
fs5 <- FeatureScatter(cs20, feats[1], feats[2], group.by = grp, raster = T, slot = "counts", pt.size = pt_size, jitter = jit, cols = pal) + NoLegend() + ggtitle("This study")
fs3 <- FeatureScatter(mv, mmfeats[1], mmfeats[2], group.by = grp, raster = T, slot = "counts", pt.size = pt_size, jitter = jit, cols = pal) + NoLegend() + ggtitle("Vladoiu et al. 2019")
fs4 <- FeatureScatter(mm_lin, mmfeats[1], mmfeats[2], group.by = grp, raster = T, slot = "counts", pt.size = pt_size, jitter = jit, cols = pal) + NoLegend() + ggtitle("La Manno et al. 2021")


# print
pdf(glue::glue("{fs_out}/all_{paste(feats, collapse = '_')}.pdf"), h = 5, w = 25)
print(fs5 | fs1 | fs2 | fs3 | fs4)
dev.off()

# dimplots
fd1 <- DimPlot(so, group.by = grp, reduction = "umap", raster = T, cols = pal) + NoAxes() + NoLegend() + ggtitle("Aldinger et al. 2021")
fd2 <- DimPlot(lin, group.by = grp, reduction = "umap", raster = T, cols = pal) + NoAxes() + NoLegend() + ggtitle("Braun et al. 2022")
fd5 <- DimPlot(cs20, group.by = grp, reduction = "umap", raster = T, cols = pal) + NoAxes() + NoLegend() + ggtitle("This study")
fd3 <- DimPlot(mv, group.by = grp, reduction = "umap", raster = T, cols = pal) + NoAxes() + NoLegend() + ggtitle("Vladoiu et al. 2019")
fd4 <- DimPlot(mm_lin, group.by = grp, reduction = "umap", raster = T, cols = pal) + NoAxes() + NoLegend() + ggtitle("La Manno et al. 2021")

# print
pdf(glue::glue("{fs_out}/all_dimplots_{paste(feats, collapse = '_')}.pdf"), h = 5, w = 25)
print(fd5 | fd1 | fd2 | fd3 | fd4)
dev.off()

feats <- c("SKOR2", "MKI67", "FOXP2", "CALB1")
# featureplots
DefaultAssay(so) <- "SCT"
fp1 <- FeaturePlot(so, feats, reduction = "umap", raster = T, ncol = 1) & NoAxes() & NoLegend()# & ggtitle(glue::glue("Aldinger et al. 2021"))
fp2 <- FeaturePlot(lin, feats, reduction = "umap", raster = T, ncol = 1) & NoAxes() & NoLegend()# & ggtitle("Braun et al. 2022")
fp5 <- FeaturePlot(cs20, feats, reduction = "umap", raster = T, ncol = 1) & NoAxes() & NoLegend()# & ggtitle("This study")
fp6 <- FeaturePlot(lake, feats, reduction = "umap", raster = T, ncol = 1) & NoAxes() & NoLegend() #& ggtitle("Lake et al. 2017")
fp7 <- FeaturePlot(sil, feats, reduction = "umap", raster = T, ncol = 1) & NoAxes() & NoLegend() #& ggtitle("Siletti et al. 2023")
fp3 <- FeaturePlot(mv, mmfeats, reduction = "umap", raster = T, ncol = 1) & NoAxes() & NoLegend() #& ggtitle("Vladoiu et al. 2019")
fp4 <- FeaturePlot(mm_lin, mmfeats, reduction = "umap", raster = T, ncol = 1) & NoAxes() & NoLegend()#& ggtitle("La Manno et al. 2021")

# print
pdf(glue::glue("{fs_out}/all_featureplots_{paste(feats, collapse = '_')}.pdf"), h = 20, w = 25)
print(fp5 | fp1 | fp2 | fp6 | fp7)
dev.off()

# vlnplots
vp1 <- VlnPlot(so, feats, reduction = "umap", raster = T, ncol = 1) & NoAxes() & NoLegend()# & ggtitle(glue::glue("Aldinger et al. 2021"))
vp2 <- VlnPlot(lin, feats, reduction = "umap", raster = T, ncol = 1) & NoAxes() & NoLegend()# & ggtitle("Braun et al. 2022")
vp5 <- VlnPlot(cs20, feats, reduction = "umap", raster = T, ncol = 1) & NoAxes() & NoLegend()# & ggtitle("This study")

# stacked barplot quantification
bpdf <- purrr::map(c(so, lin, cs20, mv, mm_lin), \(x){
  
  # human vs mouse
  if("Prdm13" %in% rownames(x@assays$RNA@counts)){ feats_loop <- mmfeats } else { feats_loop <- feats }

  # subset to feats of interest
  tmp <- x@assays$RNA@counts[feats_loop, ] %>% 
    as.matrix() %>% 
    t() %>% 
    as.data.frame()
  
  # get barcodes
  tmp_double_pos <- tmp %>% dplyr::filter(dplyr::if_all(dplyr::everything(), ~ . != 0)) %>% rownames()
  tmp_feat1_pos <- tmp %>% dplyr::filter(!!rlang::sym(feats_loop[1]) > 0) %>% rownames()
  tmp_feat2_pos <- tmp %>% dplyr::filter(!!rlang::sym(feats_loop[2]) > 0) %>% rownames()
  tmp_feat1_single_pos <- setdiff(tmp_feat1_pos, tmp_double_pos)
  tmp_feat2_single_pos <- setdiff(tmp_feat2_pos, tmp_double_pos)
  
  # get study name
  study_name <- dplyr::case_when(
    "01k" %in% x@meta.data$orig.ident ~ "ald",
    "CS21" %in% x@meta.data$orig.ident ~ "cs20",
    "10X104" %in% x@meta.data$orig.ident ~ "lin",
    "10X_MC" %in% x@meta.data$orig.ident ~ "mv",
    "10X14" %in% x@meta.data$orig.ident ~ "mm_lin")

  # build df of barcode lengths
  data.frame(group = c("double_pos", feats_loop[1], feats_loop[2]),
              size = c(length(tmp_double_pos), length(tmp_feat1_single_pos), length(tmp_feat2_single_pos))) %>% 
              dplyr::mutate(study = study_name, group = toupper(group)) %>% return()
  
}) %>% data.table::rbindlist() %>% 
dplyr::mutate(study = factor(study, levels = c("cs20", "lin", "ald", "mm_lin", "mv")),
              group = factor(group, levels = c(feats[2], feats[1], "DOUBLE_POS")))

# stacked barplot
ggbp <- ggplot(bpdf, aes(study, size, fill = group))+
  geom_bar(position = "fill", stat = "identity", color = "black")+ #position = "stack" 
  scale_fill_manual(values = c("cadetblue1", "cadetblue4", "deeppink3"),  #c("cadetblue1", "cadetblue4", "deeppink3") # c("cadetblue1", "steelblue3", "#D4A6C8")
                  labels = c(glue::glue("{feats[2]}+"), glue::glue("{feats[1]}+"), glue::glue("{feats[1]}+ {feats[2]}+")))+
  scale_x_discrete(labels=c("This study", "Braun", "Aldinger", "La Manno", "Vladoiu"))+
  theme_classic()+
  labs(y = "Cells (Prop.)", x = "", fill = "")+
  theme(axis.text.x=element_text(colour="black", size = 10), #, angle = 45, hjust = 1, vjust = 1
    axis.text.y=element_text(colour="black"),
    legend.key.size = unit(5, "pt"))


pdf(glue::glue("{fs_out}/all_{paste(feats, collapse = '_')}_barplot.pdf"), h = 2, w = 6)
print(ggbp)
dev.off()

# stacked barplot - absolute
if("SKOR2" %in% feats){
  bppal <-  c("cadetblue1", "cadetblue4", "deeppink3")} else {
    bppal <- c("lightskyblue", "cadetblue4", "#D4A6C8")
  }
ggbpa <- ggplot(bpdf, aes(study, size, fill = group))+
  geom_bar(position = "stack", stat = "identity", color = "black")+ #position = "stack" 
  scale_fill_manual(values = bppal,
                  labels = c(glue::glue("{feats[2]}+"), glue::glue("{feats[1]}+"), glue::glue("{feats[1]}+ {feats[2]}+")))+
  scale_x_discrete(labels=c("This study", "Braun", "Aldinger", "La Manno", "Vladoiu"))+
  theme_classic()+
  labs(y = "Cells (n)", x = "", fill = "")+
  theme(axis.text.x=element_text(colour="black", size = 10), #, angle = 45, hjust = 1, vjust = 1
    axis.text.y=element_text(colour="black"),
    legend.position = c(0.8, 0.8))


pdf(glue::glue("{fs_out}/absolute_barplot_{paste(feats, collapse = '_')}.pdf"), h = 4, w = 6)
print(ggbpa)
dev.off()

# stacked barplot - without mki67+ single pos
ggbp_single <- ggplot(bpdf %>% dplyr::filter(group != "MKI67"), aes(study, size, fill = group))+
  geom_bar(position = "fill", stat = "identity", color = "black")+ #position = "stack" 
  scale_fill_manual(values = c("cadetblue1", "deeppink3"),  #c("cadetblue1", "cadetblue4", "deeppink3") # c("cadetblue1", "steelblue3", "#D4A6C8")
                  labels = c(glue::glue("{feats[2]}+"), glue::glue("{feats[1]}+ {feats[2]}+")))+
  scale_x_discrete(labels=c("This study", "Braun", "Aldinger", "La Manno", "Vladoiu"))+
  theme_classic()+
  labs(y = "Cells (Prop.)", x = "", fill = "")+
  theme(axis.text.x=element_text(colour="black", size = 8, angle = 45, hjust = 1, vjust = 1), #, angle = 45, 
    axis.text.y=element_text(colour="black"),
    legend.key.size = unit(5, "pt"))


pdf(glue::glue("{fs_out}/all_{paste(feats, collapse = '_')}_barplot_no_mki67_single_pos.pdf"), h = 6, w = 8)
print(ggbp_single)
dev.off()

# prop test table
bp_mvh <- bpdf %>%
  dplyr::filter(group != "MKI67") %>% 
  dplyr::mutate(species = ifelse(study %in% c("ald", "cs20", "lin"), "Human", "Mouse")) %>% 
  dplyr::group_by(group, species) %>% 
  dplyr::summarise(sum_size = sum(size))

# prop test skor2
prop.test(x = c(3520, #human double pos
                8), #mouse double pos
          n = c(3520 + 58368, #human all skor2+
                8 + 1142), #mouse all skor2+
          alternative = "two.sided",
          correct = FALSE) #since no cells expected value <= 5 

# prop test prdm13
prop.test(x = c(4198, #human double pos
                374), #mouse double pos
          n = c(4198 + 13009, #human all prdm13+
                374 + 1161), #mouse all prdm13+
          alternative = "two.sided",
          correct = FALSE) #since no cells expected value <= 5 


pdf(glue::glue("{fs_out}/all_{paste(feats, collapse = '_')}_barplot_no_mki67_single_pos_mouse_v_human.pdf"), h = 2, w = 3)
print(ggbp_single_mvh)
dev.off()

# pax6 and ptf1a
mm1 <- FeaturePlot(mv, c("Pax6", "Ptf1a"), raster = T, ncol = 1) & NoAxes() & NoLegend()
mm2 <- FeaturePlot(mm_lin, c("Pax6", "Ptf1a"), raster = T, ncol = 1) & NoAxes() & NoLegend()
pdf("sandbox/pax6_ptf1a_mm.pdf")
mm1 | mm2
dev.off()

# timepoint breakdown---
csmd <- cs20@meta.data %>% tibble::rownames_to_column("time_bc") %>% tidyr::separate(time_bc, c("time", "bc"), sep = "_")
table(so$feat_status, so$age)
table(lin$feat_status, lin$Age)
table(csmd$feat_status, csmd$time)
table(mv$feat_status, mv$Donor)
table(mm_lin$feat_status, mm_lin$age_day)



# Braun subset ----------------------------------------------------------------

# import
lin <- qs::qread(glue::glue("{yml$lin_dir}/hb_cb_final_v3.qs"))

# check
pdf("sandbox/lin_dimplot.pdf", h = 5, w = 5)
DimPlot(lin, group.by = "seurat_clusters", raster = T, label = T, cols = c(tableau20, "lightgrey", "black")) + NoAxes() + ggtitle("")
dev.off()

ss_cells <- lin@meta.data %>% dplyr::filter(seurat_clusters %in% c(0,1,6,14)) %>% rownames() %>% sample(2000)
dims <- 20
res <- 0.2
ss <- lin %>% subset(cells = ss_cells) %>% 
  AddMetaData(lin@meta.data) %>% 
  NormalizeData(normalization.method = "LogNormalize", scale.factor = 10000) %>% 
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>% 
  ScaleData(vars.to.regress = "cc_difference") %>% 
  RunPCA() %>% 
  harmony::RunHarmony(group.by.vars = "Donor",
             reduction = "pca",
             assay.use = "RNA",
             theta = 2) %>%
  FindNeighbors(dims = 1:dims, reduction = "harmony") %>%
  FindClusters(resolution = res) %>% 
  RunUMAP(dims = 1:dims, reduction = "harmony")
pdf("sandbox/lin_ss_dimplot.pdf", h = 5, w = 10)
DimPlot(ss, group.by = c("seurat_clusters", "orig.ident"), raster = T, cols = c(tableau20, "lightgrey", "black")) & NoAxes() & ggtitle("")
dev.off()
pdf("sandbox/lin_ss_fp.pdf", h = ceiling(length(genes) / 5) * 4, w = 20)
FeaturePlot(ss, genes, raster = T) & NoAxes() & NoLegend()
dev.off()
pdf("sandbox/lin_ss_fp.pdf", h = 4, w = 8)
FeaturePlot(ss, c("PAX2", "SKOR2"), raster = T) & NoAxes() & NoLegend()
dev.off()
pdf("sandbox/lin_pax2_skor2_orderF.pdf", h = 4, w = 8)
FeaturePlot(lin, c("PAX2", "SKOR2"), raster = T, order = F) & NoAxes() & NoLegend()
dev.off()

# Aldinger subset -------------------------------------------------------------

# ald
ald <- readRDS(glue::glue("{yml$ald_path}/cbl_integrated_cleanCC_210111.rds"))
DefaultAssay(ald) <- "SCT"
pdf("sandbox/ald_pax2_skor2_orderF.pdf", h = 4, w = 8)
FeaturePlot(ald, c("PAX2", "SKOR2"), raster = T, order = F) & NoAxes() & NoLegend()
dev.off()

