# Info --------------------------------------------------------------------
# Anders E.
# May 24, 2024
# Taylor lab

# Notes -------------------------------------------------------------------
# 
# original script by Liam H.
# modified by me (Anders) March 19, 2023
# running locally
# following along this vignette:
# https://bioconductor.org/packages/devel/workflows/vignettes/GeoMxWorkflows/inst/doc/GeomxTools_RNA-NGS_Analysis.html
# also check out: https://davislaboratory.github.io/standR/articles/standR_introduction.html


# Libraries ---------------------------------------------------------------
library(GeomxTools)
library(GeoDiff)
library(NanoStringNCTools)
library(lme4)
library(lmerTest)
library(ggplot2)
library(plotly)
library(gplots)
library(RColorBrewer)
library(plyr)
library(dplyr)
library(viridis)
library(Biobase)
library(reshape2)
library(EnvStats)
library(ggiraph)
library(scales) 
library(reshape2)
library(cowplot) 
library(umap)
library(clipr)
library(pheatmap) 
library(ggrepel) 
library(ggsignif) 
library(magrittr)
library(ComplexHeatmap)
source("src/scripts/utils.R")

# Inputs ------------------------------------------------------------------

# outdir
outdir <- "out/geomx"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# geomx dir
geomx_dir <- "../g4mb_multiome/src/metadata/from_lh_20230316"


# pkc
pkc_dir <- glue::glue("{geomx_dir}/PKC")
# unzip(list.files(path = pkc_dir, pattern = ".zip$", full.names = TRUE), exdir = pkc_dir)
pkc_files <- list.files(path = pkc_dir, pattern = ".pkc$", full.names = TRUE)

# annot
annot_file <- glue::glue("{geomx_dir}/M-726 All Data WTA.xlsx")
annots <- annot_file %>% 
  readxl::read_excel() %>% 
  dplyr::mutate(low_anno = full_anno %>%  
                  stringr::str_replace("^CS-|^CS23-", "") %>% 
                  stringr::str_replace("-Medial$|-Lateral$", "") %>% 
                  stringr::str_replace("1|2", "") %>% 
                  stringr::str_replace("-LMX1A|-EOMES-KI-", "") %>% 
                  stringr::str_replace("RLSVZ", "RL-SVZ") %>% 
                  stringr::str_replace("SVZ-Ki67", "SVZ"))

# get CS medial samples
selects <- annots %>%
  dplyr::mutate(select = ifelse(Position == "Medial" & grepl("CS", Age), TRUE, FALSE)) %>% 
  dplyr::filter(select) %>%
  dplyr::pull(DCC_name)


# dcc
dcc_files <- glue::glue("{geomx_dir}/DCC/{selects}") %>% as.character()


# Import ------------------------------------------------------------------

# Create object
dat <-readNanoStringGeoMxSet(dccFiles = dcc_files,
                             pkcFiles = pkc_files[[1]],
                             phenoDataFile = annot_file,
                             phenoDataSheet = "SegmentProperties",
                             phenoDataDccColName = "DCC_name")


# QC ----------------------------------------------------------------------

# add pseudocount + 1
dat <- shiftCountsOne(dat, useDALogic = TRUE)


# Probe QC 
dat <- setBioProbeQCFlags(dat, 
                          qcCutoffs = list(minProbeRatio = 0.1,
                                           percentFailGrubbs = 20), 
                          removeLocalOutliers = TRUE)

ProbeQCResults <- fData(dat)[["QCFlags"]]

# Define QC table for Probe QC
qc_df <- data.frame(Passed = sum(rowSums(ProbeQCResults[, -1]) == 0),
                    Global = sum(ProbeQCResults$GlobalGrubbsOutlier),
                    Local = sum(rowSums(ProbeQCResults[, -2:-1]) > 0
                                & !ProbeQCResults$GlobalGrubbsOutlier))

# Filter to passing probes
ProbeQCPassed <- subset(dat, 
                        fData(dat)[["QCFlags"]][,c("LowProbeRatio")] == FALSE &
                          fData(dat)[["QCFlags"]][,c("GlobalGrubbsOutlier")] == FALSE)
dat <- ProbeQCPassed 



# Aggregate counts --------------------------------------------------------

# agg
target_dat <- aggregateCounts(dat)


# QC,  cont. --------------------------------------------------------------

# modules
modules <- annotation(dat) %>% stringr::str_replace(".pkc", "")

# LOQ---
cutoff <- 2
minLOQ <- 2
LOQ <- data.frame(row.names = colnames(target_dat))

# assign minLOQ per segment
for(module in modules) {
  
  vars <- paste0(c("NegGeoMean_", "NegGeoSD_"), module)
  
  if(all(vars[1:2] %in% colnames(pData(target_dat)))) {
    
    LOQ[, module] <- pmax(minLOQ, pData(target_dat)[, vars[1]] * pData(target_dat)[, vars[2]] ^ cutoff)
    
  }
}

# store
pData(target_dat)$LOQ <- LOQ


# LOQ matrix
LOQ_Mat <- c()
for(module in modules) {
  ind <- fData(target_dat)$Module == module
  Mat_i <- t(esApply(target_dat[ind, ], MARGIN = 1,
                     FUN = function(x) {
                       x > LOQ[, module]
                     }))
  LOQ_Mat <- rbind(LOQ_Mat, Mat_i)
}

# reorder
LOQ_Mat <- LOQ_Mat[fData(target_dat)$TargetName, ]

# genes detected per segment
pData(target_dat)$GenesDetected <- colSums(LOQ_Mat, na.rm = TRUE)
pData(target_dat)$GeneDetectionRate <- pData(target_dat)$GenesDetected / nrow(target_dat)

# Determine detection thresholds: 1%, 5%, 10%, 15%, >15%
pData(target_dat)$DetectionThreshold <- 
  cut(pData(target_dat)$GeneDetectionRate,
      breaks = c(0, 0.01, 0.05, 0.1, 0.15, 1),
      labels = c("<1%", "1-5%", "5-10%", "10-15%", ">15%"))

# stacked bar plot of different cut points (1%, 5%, 10%, 15%)
qc1 <- ggplot(pData(target_dat),
              aes(x = DetectionThreshold)) +
  geom_bar(aes(fill = Position)) +
  geom_text(stat = "count", aes(label = ..count..), vjust = -0.5) +
  theme_classic() +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  scale_fill_manual(values = tableau10[c(1,3,4)]) +
  labs(x = "Gene Detection Rate",
       y = "Segments, #",
       fill = "Segment Type")

# violinplot
qc2 <- ggplot(pData(target_dat), aes(Position, GeneDetectionRate, color = Position))+
  geom_violin(width=0.8)+
  geom_boxplot(width=0.2, color="black", alpha = 0.2, outlier.shape = NA)+
  geom_jitter(width=0.1, alpha = 1)+
  ggpubr::stat_compare_means(method = "t.test", comparisons = list(c("Lateral", "Medial"), c("Lateral", "RL"), c("Medial", "RL")))+
  theme_classic()+
  scale_color_manual(values = tableau10[c(1,3,4)])+
  ylab('Gene Detection Rate')+
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.text.x = element_text(colour = "black"))

# Gene detection rate barplot
LOQ_Mat <- LOQ_Mat[, colnames(target_dat)]
fData(target_dat)$DetectedSegments <- rowSums(LOQ_Mat, na.rm = TRUE)
fData(target_dat)$DetectionRate <-
  fData(target_dat)$DetectedSegments / nrow(pData(target_dat))


plot_detect <- data.frame(Freq = c(1, 5, 10, 20, 30, 50))
plot_detect$Number <-
  unlist(lapply(c(0.01, 0.05, 0.1, 0.2, 0.3, 0.5),
                function(x) {sum(fData(target_dat)$DetectionRate >= x)}))
plot_detect$Rate <- plot_detect$Number / nrow(fData(target_dat))
rownames(plot_detect) <- plot_detect$Freq

qc3 <- ggplot(plot_detect, aes(x = as.factor(Freq), y = Rate, fill = Rate)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = formatC(Number, format = "d", big.mark = ",")),
            vjust = 1.6, color = "black", size = 4) +
  scale_fill_gradient2(low = "orange2", mid = "lightblue",
                       high = "dodgerblue3", midpoint = 0.65,
                       limits = c(0,1),
                       labels = scales::percent) +
  theme_bw() +
  scale_y_continuous(labels = scales::percent, limits = c(0,1),
                     expand = expansion(mult = c(0, 0))) +
  labs(x = "% of Segments",
       y = "Genes Detected, % of Panel > LOQ")+
  theme(legend.position = "none")

# plot
pdf(glue::glue("{outdir}/qc_gene_detection.pdf"), h = 3.5, w = 10)
qc1 + qc2 + qc3 + patchwork::plot_layout(widths = c(1,2,2))
dev.off()


# threshold for gene detection
gene_thresh <- 0.5 # i.e. 50%

# manually retain negative control probe
negativeProbefData <- subset(fData(target_dat), CodeClass == "Negative")
neg_probes <- unique(negativeProbefData$TargetName)

# how many features retained, exactly?
dim(target_dat[fData(target_dat)$DetectionRate >= gene_thresh |
                 fData(target_dat)$TargetName %in% neg_probes, ])

# filter to gene_thresh
target_dat <- target_dat[fData(target_dat)$DetectionRate >= 0.5 |
                           fData(target_dat)$TargetName %in% neg_probes, ]


# Normalization QC -----------------------------------------------------------
# NB: there are two slides total, named:
# Slide4; 3_Med

# Graph Q3 value vs negGeoMean of Negatives
ann_of_interest <- "SlideName"
Stat_data <- data.frame(row.names = colnames(exprs(target_dat)),
                        Segment = colnames(exprs(target_dat)),
                        Annotation = pData(target_dat)[, ann_of_interest],
                        Q3 = unlist(apply(exprs(target_dat), 2,
                                          quantile, 0.75, na.rm = TRUE)),
                        NegProbe = exprs(target_dat)[neg_probes, ])
Stat_data_m <- reshape2::melt(Stat_data, measure.vars = c("Q3", "NegProbe"),
                              variable.name = "Statistic", value.name = "Value")

plt1 <- ggplot(Stat_data_m,
               aes(x = Value, fill = Statistic)) +
  geom_histogram(bins = 40) + theme_bw() +
  scale_x_continuous(trans = "log2") +
  facet_wrap(~Annotation, nrow = 1) +
  scale_fill_brewer(palette = 3, type = "qual") +
  labs(x = "Counts", y = "Segments, #")

plt2 <- ggplot(Stat_data,
               aes(x = NegProbe, y = Q3, color = Annotation)) +
  geom_abline(intercept = 0, slope = 1, lty = "dashed", color = "darkgray") +
  geom_point() + guides(color = "none") + theme_bw() +
  scale_x_continuous(trans = "log2") +
  scale_y_continuous(trans = "log2") +
  theme(aspect.ratio = 1) +
  labs(x = "Negative Probe GeoMean, Counts", y = "Q3 Value, Counts")

plt3 <- ggplot(Stat_data,
               aes(x = NegProbe, y = Q3 / NegProbe, color = Annotation)) +
  geom_hline(yintercept = 1, lty = "dashed", color = "darkgray") +
  geom_point() + theme_bw() +
  scale_x_continuous(trans = "log2") +
  scale_y_continuous(trans = "log2") +
  theme(aspect.ratio = 1) +
  labs(x = "Negative Probe GeoMean, Counts", y = "Q3/NegProbe Value, Counts")

btm_row <- cowplot::plot_grid(plt2, plt3, nrow = 1, labels = c("B", ""),
                              rel_widths = c(0.43,0.57))

# plot
pdf(glue::glue("{outdir}/qc_normalization_plots.pdf"), h = 5, w = 7)
cowplot::plot_grid(plt1, btm_row, ncol = 1, labels = c("A", ""))
dev.off()



# Normalization -----------------------------------------------------------

target_dat <- normalize(target_dat , data_type = "RNA",
                        norm_method = "quant",
                        desiredQuantile = .75,
                        toElt = "q_norm")

boxplot(exprs(target_dat)[,1:10],
        col = "#9EDAE5", main = "Raw Counts",
        log = "y", names = 1:10, xlab = "Segment",
        ylab = "Counts, Raw")
boxplot(assayDataElement(target_dat[,1:10], elt = "q_norm"),
        col = "#2CA02C", main = "Q3 Norm Counts",
        log = "y", names = 1:10, xlab = "Segment",
        ylab = "Counts, Q3 Normalized")

# save
counts_outdir <- glue::glue("{outdir}/counts")
dir.create(counts_outdir, recursive = T, showWarnings = F)

qs::qsave(target_dat, glue::glue("{counts_outdir}/geomx_target_dat.qs"))
qs::qsave(list("counts" = exprs(target_dat),
               "norm_counts" = assayDataElement(target_dat, elt = "q_norm")),
          glue::glue("{counts_outdir}/geomx_counts.qs"))

counts <- qs::qread(glue::glue("{counts_outdir}/geomx_counts.qs"))
target_dat <- qs::qread(glue::glue("{counts_outdir}/geomx_target_dat.qs"))



# Heatmap CS medial ------------------------------------------------------------------


# calc CV function
calc_CV <- function(x) {sd(x) / mean(x)}

# Low res anno: 
pData(target_dat)$low_anno <- pData(target_dat)$full_anno %>% 
  stringr::str_replace("CS23-SVZ", "CS23-oSVZ") %>% 
  stringr::str_replace("^CS-|^CS23-", "") %>% 
  stringr::str_replace("-Medial$|-Lateral$", "") %>% 
  stringr::str_replace("1|2", "") %>% 
  stringr::str_replace("-LMX1A|-EOMES-KI-", "") %>% 
  stringr::str_replace("RLSVZ", "RL-SVZ") %>% 
  stringr::str_replace("SVZ-Ki67", "SVZ")


# log2
assayDataElement(object = target_dat, elt = "log_q") <- assayDataApply(target_dat, 2, FUN = log, base = 2, elt = "q_norm")

# subset
selects <- pData(target_dat) %>% dplyr::mutate(select = ifelse(Position == "Medial" & 
                                                               grepl("CS", Age), TRUE, FALSE) &
                                                               grepl("VZ", low_anno) & 
                                                               !grepl("BS", low_anno) & 
                                                               !grepl("Ki67", full_anno)) %>% 
  dplyr::pull(select)
ss <- target_dat %>% subset(select = selects)

# goi
ss_goi <- assayDataApply(ss, elt = "log_q", MARGIN = 1, calc_CV) %>% 
  tibble::enframe() %>% 
  dplyr::slice_max(value, n = 250) %>% 
  dplyr::pull(name)

# norm expr matrix
ss_mat <- assayDataElement(ss[ss_goi,], elt = "log_q")


# heatmap color bounds
low_bound <- min(ss_mat, na.rm = TRUE)
mid_bound <- median(as.matrix(ss_mat), na.rm = TRUE) # or mean()
quant_3rd <- quantile(ss_mat, 0.75, na.rm = TRUE)
high_bound <- max(ss_mat, na.rm = TRUE)

# top annotation df
ss_ta_df <- pData(ss) %>% dplyr::select(Age, Region = low_anno) %>% dplyr::mutate(Age = stringr::str_replace(Age, "CS18", "CS19"),
                                                                                  Region = factor(Region, levels = c("VZ", "iSVZ", "oSVZ")))

# top annot palettes
age_pal <- data.frame(sort(unique(ss_ta_df$Age)), brewer.pal(9, "BuPu")[c(2,4,5)]) %>% tibble::deframe()
region_pal <- data.frame(unique(ss_ta_df$Region) %>% sort(), brewer.pal(9, "BuGn")[c(2,4,6)]) %>% tibble::deframe()




# top annot
ss_ta <- ComplexHeatmap::HeatmapAnnotation(df = ss_ta_df,
                                           col = list(Age = age_pal, Region = region_pal), #Laterality = lat_pal, 
                                           annotation_name_gp = grid::gpar(fontsize = 8),
                                           simple_anno_size = unit(0.25, "cm"),
                                           annotation_name_side = "left")


column_dend <- as.dendrogram(hclust(dist(t(ss_mat))))


# heatmap
ss_hm <- ComplexHeatmap::Heatmap(as.matrix(ss_mat),
                                 name = "Norm. expr.",
                                 top_annotation = ss_ta,
                                 col = circlize::colorRamp2(c(mid_bound, quant_3rd, high_bound), RColorBrewer::brewer.pal(9, "GnBu")[c(1,5,9)]),
                                 # column_order = ord,
                                 show_row_names = FALSE,
                                 show_column_names = FALSE,
                                 # row_names_gp = grid::gpar(fontsize = 6),
                                 cluster_rows = TRUE,
                                 cluster_columns = column_dend,
                                 clustering_distance_rows = "pearson",
                                 clustering_method_rows = "ward.D2",
                                 show_column_dend = TRUE,
                                 show_row_dend = TRUE,
                                 border = TRUE,
                                 width = unit(5, "cm"),
                                 height = unit(15, "cm"),
                                 row_title = " ",
                                 column_title = " ",
                                 row_gap = unit(0, "cm"),
                                 column_gap = unit(0, "cm"),
                                 heatmap_legend_param = list(direction = "vertical",
                                                             at = c(floor(mid_bound), ceiling(high_bound)),
                                                             labels = c(floor(mid_bound), ceiling(high_bound)),
                                                             legend_width = unit(1.75, "cm"),
                                                             title_position = "leftcenter-rot")) #leftcenter-rot

pdf(glue::glue("{outdir}/cs_medial_heatmap.pdf"), w = 4, h = 8)
ss_hm %<>% ComplexHeatmap::draw(merge_legend = TRUE, heatmap_legend_side = "right", annotation_legend_side = "right")
dev.off()

# save excel of gene list
ss_goi[ComplexHeatmap::row_order(ss_hm)] %>% 
  as.data.frame() %>% 
  writexl::write_xlsx(glue::glue("{outdir}/geomx_heatmap_gene_order.xlsx"))

# modified heatmap --------------------------------------------------------

# omitting CS23 and averaging CS19/20
sss_samples <- ss_ta_df %>% dplyr::filter(Age %in% c("CS19", "CS20")) %>% rownames()

# subset
sss_selects <- pData(target_dat) %>% dplyr::mutate(select = ifelse(rownames(.) %in% sss_samples, TRUE, FALSE)) %>% dplyr::pull(select)
sss <- target_dat %>% subset(select = sss_selects)

# need to pick top var
sss_goi <- assayDataApply(sss, elt = "log_q", MARGIN = 1, calc_CV) %>%
  tibble::enframe() %>%
  dplyr::slice_max(value, n = 250) %>%
  dplyr::pull(name)

# or top expressed
# sss_goi <- sss_mat %>% 
#   as.data.frame() %>% 
#   tibble::rownames_to_column("gene") %>% 
#   dplyr::rowwise() %>% 
#   dplyr::mutate(avg = mean(VZ, iSVZ1, iSVZ2, oSVZ)) %>%
#   dplyr::ungroup() %>% 
#   dplyr::slice_max(avg, n = 250) %>% 
#   dplyr::pull(gene)

# dict
sss_dict <- pData(ss) %>% 
  dplyr::filter(rownames(.) %in% sss_samples) %>% 
  dplyr::mutate(refined_low_anno = dplyr::case_when(
    low_anno == "iSVZ" & grepl("iSVZ1", full_anno) ~ "iSVZ1",
    low_anno == "iSVZ" & grepl("iSVZ2", full_anno) ~ "iSVZ2",
    TRUE ~ low_anno
  )) %>% tibble::rownames_to_column("name") %>% 
  dplyr::select(name, Region = refined_low_anno)


# norm expr matrix...summarised (mean) by region
sss_mat <- assayDataElement(sss[sss_goi,], elt = "log_q") %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column("gene") %>% 
  tidyr::pivot_longer(-gene) %>% 
  dplyr::left_join(sss_dict, by = "name") %>% 
  dplyr::group_by(gene, Region) %>% 
  dplyr::summarise(mean_xp = mean(value)) %>% 
  tidyr::pivot_wider(id_cols = gene, names_from = Region, values_from = mean_xp) %>% 
  tibble::column_to_rownames("gene") %>% 
  as.matrix()
  
# or redo
# assigning each gene to region based on max expression
# then subsetting to top variable genes per group
per_group <- F
if(per_group){
  tmp <- assayDataElement(sss, elt = "log_q") %>% 
    as.data.frame() %>% 
    tibble::rownames_to_column("gene") %>% 
    tidyr::pivot_longer(-gene) %>% 
    dplyr::left_join(sss_dict, by = "name") %>% 
    dplyr::group_by(gene, Region) %>% 
    dplyr::summarise(mean_xp = mean(value))
  
  max_grps <- tmp %>% 
    dplyr::group_by(gene) %>% 
    dplyr::slice_max(mean_xp)
  grp_vars <- tmp %>% 
    dplyr::group_by(gene) %>% 
    dplyr::summarise(grp_var = var(mean_xp))
  joined <- max_grps %>% 
    dplyr::left_join(grp_vars, by = "gene")
  grp_var_goi <- joined %>% 
    dplyr::group_by(Region) %>% 
    dplyr::slice_max(grp_var, n = 50) %>% 
    dplyr::pull(gene)
  
  sss_mat <- assayDataElement(sss[grp_var_goi,], elt = "log_q") %>% 
    as.data.frame() %>% 
    tibble::rownames_to_column("gene") %>% 
    tidyr::pivot_longer(-gene) %>% 
    dplyr::left_join(sss_dict, by = "name") %>% 
    dplyr::group_by(gene, Region) %>% 
    dplyr::summarise(mean_xp = mean(value)) %>% 
    tidyr::pivot_wider(id_cols = gene, names_from = Region, values_from = mean_xp) %>% 
    tibble::column_to_rownames("gene") %>% 
    as.matrix()
}


# heatmap color bounds
low_bound <- min(sss_mat, na.rm = TRUE)
mid_bound <- median(as.matrix(sss_mat), na.rm = TRUE) # or mean()
quant_3rd <- quantile(sss_mat, 0.75, na.rm = TRUE)
high_bound <- max(sss_mat, na.rm = TRUE)

# top annotation df
sss_ta_df <- data.frame(Region = c("VZ", "iSVZ1", "iSVZ2", "oSVZ") %>% factor(levels = c("VZ", "iSVZ1","iSVZ2", "oSVZ"))) %>% magrittr::set_rownames(.$Region)

# top annot
new_regionpal <- data.frame(levels(sss_ta_df$Region), brewer.pal(9, "BuGn")[c(2,4,5,6)]) %>% tibble::deframe()
sss_ta <- ComplexHeatmap::HeatmapAnnotation(df = sss_ta_df,
                                           col = list(Region = new_regionpal),
                                           annotation_name_gp = grid::gpar(fontsize = 8),
                                           simple_anno_size = unit(0.25, "cm"),
                                           annotation_name_side = "left")

# reorder
column_dend <- as.dendrogram(hclust(dist(t(sss_mat)))) %>% dendextend::rotate(c(2:4,1))

# break genes into groups
groups <- sss_mat %>% 
  as.data.frame() %>% 
  dplyr::rowwise() %>% 
  dplyr::mutate(max_column = names(.)[which.max(dplyr::c_across())]) %>% 
  as.data.frame() %>% 
  magrittr::set_rownames(rownames(sss_mat))

# heatmap
highlight_genes <- c('CASZ1', 'SOX3', 'NOTCH1', 'CDK6', 'HES5', 'SOX9', 'SOX2', 'MKI67', 'CENPF', 'ARHGAP11B', 'ASPM', 'SHROOM3', 'SLC1A3', 'PLP1', 'SMOC2', 'CRABP1', 'LHX5', 'EBF2', 'GRIK3', 'DCX', 'LHX1', 'NCAM1', 'POU2F2', 'POU3F1', 'INSM1', 'OLIG3', 'PAX5', 'MEIS2', 'SHANK3', 'SPANXA2', 'CCNA1', 'CALB1', 'KCNC3', 'EBF3', 'GRIA2', 'EBF1', 'STXBP1', 'CACNA1B', 'EN1', 'CBLN1', 'KCNIP1', 'SFRP1', 'RORA', 'KCNMA1', 'DAB1', 'RBFOX1', 'RUNX1T1', 'NRXN1', 'CALB2', 'SNCA', 'KCNJ3', 'KBTBD4', 'NKX1-2', 'NRXN1', 'FOXP2', 'UNCX', 'ASIC1', "NOTCH2NLA", 'PTPRZ1') %>% unique()
ra <- rowAnnotation(foo = anno_mark(at = which(rownames(sss_mat) %in% highlight_genes),
                                    labels = rownames(sss_mat)[rownames(sss_mat) %in% highlight_genes],
                                    labels_gp = gpar(fontsize = 8), 
                                    padding = unit(0.5, "mm")))

sss_hm <- ComplexHeatmap::Heatmap(as.matrix(sss_mat),
                                 name = "Norm. expr.",
                                 top_annotation = sss_ta,
                                 col = circlize::colorRamp2(c(mid_bound, quant_3rd, high_bound), RColorBrewer::brewer.pal(9, "GnBu")[c(1,5,9)]),
                                 row_split = factor(groups$max_column, levels = c("VZ", "iSVZ1", "iSVZ2", "oSVZ")),
                                 row_gap = unit(1, "mm"),
                                 row_title = "%s",
                                 row_title_rot = 0,
                                 row_title_gp = gpar(fontsize = 8),
                                 # column_order = ord,
                                 show_row_names = FALSE,
                                 # row_names_gp = gpar(fontsize = 4),
                                 show_column_names = FALSE,
                                 # cluster_rows = FALSE,
                                 cluster_columns = column_dend,
                                 # clustering_distance_rows = "pearson",
                                 # clustering_method_rows = "ward.D2",
                                 cluster_row_slices = FALSE,
                                 show_column_dend = TRUE,
                                 show_row_dend = FALSE,
                                 border = TRUE,
                                 width = unit(3, "cm"),
                                 height = unit(15, "cm"),
                                 column_title = " ",
                                 column_gap = unit(0, "cm"),
                                 right_annotation = ra,
                                 heatmap_legend_param = list(direction = "vertical",
                                                             at = c(floor(mid_bound), ceiling(high_bound)),
                                                             labels = c(floor(mid_bound), ceiling(high_bound)),
                                                             legend_width = unit(1.75, "cm"),
                                                             title_position = "leftcenter-rot")) #leftcenter-rot


pdf(glue::glue("{outdir}/cs_medial_heatmap_aggr_split_isvz.pdf"), w = 4, h = 8)
sss_hm %<>% ComplexHeatmap::draw(merge_legend = TRUE, heatmap_legend_side = "right", annotation_legend_side = "right")
dev.off()


sss_goi[ComplexHeatmap::row_order(sss_hm)] %>% 
  as.data.frame() %>% 
  writexl::write_xlsx(glue::glue("{outdir}/geomx_heatmap_lumped_gene_order.xlsx"))


# DE prep ----------------------------------------------------------------------


# Log2 transform
assayDataElement(object = ss, elt = "log_q") <- assayDataApply(ss, 2, FUN = log, base = 2, elt = "q_norm")


# get top variable genes--- 
CV_dat <- assayDataApply(ss, elt = "log_q", MARGIN = 1, calc_CV)


# sneak peek
CV_dat %>% 
  tibble::enframe() %>% 
  dplyr::slice_max(value, n = 10)

# Identify genes in the top 3rd of the CV values
GOI <- names(CV_dat)[CV_dat > quantile(CV_dat, 0.8)]



# Inter-region DE ---------------------------------------------------------

pData(target_dat) %<>% dplyr::mutate(refined_low_anno = dplyr::case_when(
  low_anno == "iSVZ" & grepl("iSVZ1", full_anno) ~ "iSVZ1",
  low_anno == "iSVZ" & grepl("iSVZ2", full_anno) ~ "iSVZ2",
  TRUE ~ low_anno),
  Age = stringr::str_replace(Age, "CS18", "CS19"))

# WITHIN-SLIDE analysis:
# Volcano plots for VZ vs iSVZ1, VZ vs iSVZ combined, VZ vs oSVZ


# subset
de_dat <- target_dat %>% 
  subset(select = (pData(target_dat)[["Age"]] %in% c("CS19", "CS20")) & 
           (pData(target_dat)[["refined_low_anno"]] %in% c("VZ", "iSVZ1", "iSVZ2", "oSVZ")) & 
           (pData(target_dat)[["Position"]] %in% c("Medial")))


# convert test variable to factor
pData(de_dat)$testClass <- factor(pData(de_dat)$refined_low_anno, levels = c("VZ", "iSVZ1", "iSVZ2", "oSVZ"))

# run LMM:
tictoc::tic()
mixedOutmc <- mixedModelDE(de_dat,
                           elt = "log_q",
                           modelFormula = ~ testClass + (1 + testClass | SlideName),
                           groupVar = "testClass",
                           nCores = 6,
                           multiCore = TRUE)
tictoc::toc()
qs::qsave(mixedOutmc, glue::glue("{outdir}/mixedOutmc.qs"))
mixedOutmc <- qs::qread(glue::glue("{outdir}/mixedOutmc.qs"))

# format results as data.frame
r_test <- do.call(rbind, mixedOutmc["lsmeans", ])
tests <- rownames(r_test)
r_test <- as.data.frame(r_test)
r_test$Contrast <- tests

# use lapply in case you have multiple levels of your test factor to
# correctly associate gene name with its row in the results table
r_test$Gene <- unlist(lapply(colnames(mixedOutmc), rep, nrow(mixedOutmc["lsmeans", ][[1]])))
r_test$Subset <- "CS-VZ"
r_test$FDR <- p.adjust(r_test$`Pr(>|t|)`, method = "fdr")
de_df <- r_test[, c("Gene", "Subset", "Contrast", "Estimate", "Pr(>|t|)", "FDR")] %>% 
  tibble::as_tibble() %>% 
  dplyr::rename(p = `Pr(>|t|)`) %>% 
  janitor::clean_names()


# Volcano plots ---

# setup
padj_thresh <- 0.1
est_thresh <- 0.25
n <- 50

for(i in unique(de_df$contrast)){
 
  # parse
  comps <- i %>% stringr::str_split(" - ") %>% unlist()
  right <- comps[1]
  left <- comps[2]
  
  # add color and labels
  vdf <- de_df %>% 
    dplyr::filter(contrast == i) %>% 
    dplyr::mutate(color = dplyr::case_when(
      fdr > padj_thresh | abs(estimate) < est_thresh ~ "insig",
      estimate > est_thresh ~ right,
      estimate < -est_thresh ~ left),
      label = ifelse(color != "insig", gene, NA))
  
  # quick save
  vdf %>% dplyr::arrange(fdr) %>% writexl::write_xlsx(glue::glue("{outdir}/degs_volcano_{i}.xlsx"))
  
  # x bound
  xbound <- max(abs(vdf$estimate)) + 0.25
  
  # plot
  v1 <- ggplot(vdf, aes(estimate, -log10(fdr), color = color, label = label))+
    geom_point(size = 3)+
    theme_classic()+
    scale_color_manual(breaks = c(right, left), values = c(tableau10[c(1,3)], "lightgrey"), na.value="lightgrey")+
    ggrepel::geom_text_repel(size = 3, segment.size = 0.1, 
                             box.padding = 0.4, color = "black", 
                             segment.color = "black", min.segment.length = 0, 
                             max.overlaps = 10)+
    labs(x = bquote(~log[2]~FC), y = bquote(~-log[10]~p[adj]), color = "")+
    xlim(-xbound, xbound)+
    # theme(legend.position = c(0.1,0.8))+
    guides(colour = guide_legend(override.aes = list(size=3)))
  
  # save
  pdf(glue::glue("{outdir}/volcano_{i}.pdf"), h = 8, w = 8)
  print(v1)
  dev.off()
  
}




# Inter-region DE modified ------------------------------------------------

# convert test variable to factor
pData(de_dat)$testClass <- factor(pData(de_dat)$low_anno, levels = c("VZ", "iSVZ", "oSVZ"))

# run LMM:
tictoc::tic()
mixedOutmc <- mixedModelDE(de_dat,
                           elt = "log_q",
                           modelFormula = ~ testClass + (1 + testClass | SlideName),
                           groupVar = "testClass",
                           nCores = 6,
                           multiCore = TRUE)
tictoc::toc()
qs::qsave(mixedOutmc, glue::glue("{outdir}/mixedOutmc_aggr.qs"))

# format results as data.frame
r_test <- do.call(rbind, mixedOutmc["lsmeans", ])
tests <- rownames(r_test)
r_test <- as.data.frame(r_test)
r_test$Contrast <- tests

# use lapply in case you have multiple levels of your test factor to
# correctly associate gene name with its row in the results table
r_test$Gene <- unlist(lapply(colnames(mixedOutmc), rep, nrow(mixedOutmc["lsmeans", ][[1]])))
r_test$Subset <- "CS-VZ"
r_test$FDR <- p.adjust(r_test$`Pr(>|t|)`, method = "fdr")
de_df <- r_test[, c("Gene", "Subset", "Contrast", "Estimate", "Pr(>|t|)", "FDR")] %>% 
  tibble::as_tibble() %>% 
  dplyr::rename(p = `Pr(>|t|)`) %>% 
  janitor::clean_names()


# Volcano plots ---

# setup
padj_thresh <- 0.1
est_thresh <- 0.25
n <- 50

for(i in unique(de_df$contrast)){
  
  # parse
  comps <- i %>% stringr::str_split(" - ") %>% unlist()
  right <- comps[1]
  left <- comps[2]
  
  # add color and labels
  vdf <- de_df %>% 
    dplyr::filter(contrast == i) %>% 
    dplyr::mutate(color = dplyr::case_when(
      fdr > padj_thresh | abs(estimate) < est_thresh ~ "insig",
      estimate > est_thresh ~ right,
      estimate < -est_thresh ~ left),
      label = ifelse(color != "insig", gene, NA))
  
  # x bound
  xbound <- max(abs(vdf$estimate)) + 0.25
  
  # plot
  v1 <- ggplot(vdf, aes(estimate, -log10(fdr), color = color, label = label))+
    geom_point(size = 3)+
    theme_classic()+
    scale_color_manual(breaks = c(right, left), values = c(tableau10[c(1,3)], "lightgrey"), na.value="lightgrey")+
    ggrepel::geom_text_repel(size = 3, segment.size = 0.1, 
                             box.padding = 0.4, color = "black", 
                             segment.color = "black", min.segment.length = 0, 
                             max.overlaps = 10)+
    labs(x = bquote(~log[2]~FC), y = bquote(~-log[10]~p[adj]), color = "")+
    xlim(-xbound, xbound)+
    # theme(legend.position = c(0.1,0.8))+
    guides(colour = guide_legend(override.aes = list(size=3)))
  
  # save
  pdf(glue::glue("{outdir}/volcano_aggr_{i}.pdf"), h = 8, w = 8)
  print(v1)
  dev.off()
  
}





# Same-region heatmap ----------------------------------------------------------



# averaging by age and region
pData(target_dat) %<>% dplyr::mutate(svz_lumped = dplyr::case_when(
  low_anno == "VZ" ~ "VZ",
  grepl("Ki67", low_anno) ~ "ki67_SVZ",
  grepl("SVZ", low_anno) ~ "SVZ",
  TRUE ~ low_anno),
  age_lumped = dplyr::case_when(
    Age %in% c("CS19", "CS20") ~ "CS19/20",
    TRUE ~ Age))


for(i in c("VZ", "SVZ")){
  
  # subset
  selects <- pData(target_dat) %>% dplyr::mutate(select = ifelse(age_lumped %in% c("CS19/20", "CS23") & svz_lumped == i, TRUE, FALSE)) %>% dplyr::pull(select)
  s <- target_dat %>% subset(select = selects)
  
  # dict
  s_dict <- pData(s) %>% 
    tibble::rownames_to_column("name") %>% 
    dplyr::select(name, Age = age_lumped)
  
  
  # norm expr matrix...summarised (mean) by region
  mat <- assayDataElement(s, elt = "log_q") %>% 
    as.data.frame() %>% 
    tibble::rownames_to_column("gene") %>% 
    tidyr::pivot_longer(-gene) %>% 
    dplyr::left_join(s_dict, by = "name") %>% 
    dplyr::group_by(gene, Age) %>% 
    dplyr::summarise(mean_xp = mean(value)) %>% 
    tidyr::pivot_wider(id_cols = gene, names_from = Age, values_from = mean_xp) %>% 
    tibble::column_to_rownames("gene") %>% 
    as.matrix()
  
  
  # goi
  diffs <- mat %>% as.data.frame() %>% dplyr::mutate(diff = CS23 - `CS19/20`)
  goi <- c(diffs %>% dplyr::slice_max(diffs, n = 250) %>% rownames(),
           diffs %>% dplyr::slice_min(diffs, n = 250) %>% rownames())
  
  # heatmap color bounds
  low_bound <- min(mat, na.rm = TRUE)
  mid_bound <- median(as.matrix(mat), na.rm = TRUE) # or mean()
  quant_3rd <- quantile(mat, 0.75, na.rm = TRUE)
  high_bound <- max(mat, na.rm = TRUE)
  
  # top annotation df
  ta_df <- data.frame(Age = c("CS19/20", "CS23") %>% factor(levels = c("CS19/20", "CS23"))) %>% magrittr::set_rownames(.$Age)
  
  # top annot
  agepal <- data.frame(levels(ta_df$Age), brewer.pal(9, "BuGn")[c(2,6)]) %>% tibble::deframe()
  ta <- ComplexHeatmap::HeatmapAnnotation(df = ta_df,
                                              col = list(Region = agepal),
                                              annotation_name_gp = grid::gpar(fontsize = 8),
                                              simple_anno_size = unit(0.25, "cm"),
                                              annotation_name_side = "left")
  
  # reorder
  column_dend <- as.dendrogram(hclust(dist(t(mat)))) %>% rotate(c(2:4,1))
  
  # heatmap
  hm <- ComplexHeatmap::Heatmap(as.matrix(mat[goi,]),
                                    name = "Norm. expr.",
                                    top_annotation = ta,
                                    col = circlize::colorRamp2(c(mid_bound, quant_3rd, high_bound), RColorBrewer::brewer.pal(9, "GnBu")[c(1,5,9)]),
                                    # column_order = ord,
                                    show_row_names = FALSE,
                                    show_column_names = FALSE,
                                    cluster_rows = TRUE,
                                    cluster_columns = TRUE,
                                    clustering_distance_rows = "pearson",
                                    clustering_method_rows = "ward.D2",
                                    show_column_dend = TRUE,
                                    show_row_dend = TRUE,
                                    border = TRUE,
                                    width = unit(3, "cm"),
                                    height = unit(15, "cm"),
                                    row_title = " ",
                                    column_title = " ",
                                    row_gap = unit(0, "cm"),
                                    column_gap = unit(0, "cm"),
                                    heatmap_legend_param = list(direction = "vertical",
                                                                at = c(floor(mid_bound), ceiling(high_bound)),
                                                                labels = c(floor(mid_bound), ceiling(high_bound)),
                                                                legend_width = unit(1.75, "cm"),
                                                                title_position = "leftcenter-rot")) #leftcenter-rot
  
  pdf(glue::glue("{outdir}/cs_medial_heatmap_age_split_{i}.pdf"), w = 4, h = 8)
  hm %<>% ComplexHeatmap::draw(merge_legend = TRUE, heatmap_legend_side = "right", annotation_legend_side = "right")
  dev.off()
  
}



