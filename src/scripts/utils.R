# Info --------------------------------------------------------------------
# Anders E.
# Apr 24, 2023
# Taylor lab

# Libraries ---------------------------------------------------------------
require(magrittr)

# Palettes ---------------------------------------------------------------
colors_dutch <- c('#FFC312','#C4E538','#12CBC4','#FDA7DF','#ED4C67','#F79F1F','#A3CB38','#1289A7','#D980FA','#B53471','#EE5A24','#009432','#0652DD','#9980FA','#833471','#EA2027','#006266','#1B1464','#5758BB','#6F1E51')
colors_spanish <- c('#40407a','#706fd3','#f7f1e3','#34ace0','#33d9b2','#2c2c54','#474787','#aaa69d','#227093','#218c74','#ff5252','#ff793f','#d1ccc0','#ffb142','#ffda79','#b33939','#cd6133','#84817a','#cc8e35','#ccae62')
custom_colors <- c(colors_dutch, colors_spanish)
fave_pal <- c("#FB6467FF", "grey", "#197EC0FF")
phase_pal <- c("#2B83BA", "#FDAE61", "#D7191C")
tableau10 <- c("#4E79A7", "#F28E2B", "#E15759", "#76B7B2", "#59A14F", "#EDC948", "#B07AA1", "#FF9DA7", "#9C755F", "#BAB0AC") # from {ggthemes}
tableau20 <- c('#4E79A7', '#A0CBE8', '#F28E2B', '#FFBE7D', '#59A14F', '#8CD17D', '#B6992D', '#F1CE63', '#499894', '#86BCB6', '#E15759', '#FF9D9A', '#79706E', '#BAB0AC', '#D37295', '#FABFD2', '#B07AA1', '#D4A6C8', '#9D7660', '#D7B5A6')

# Gene list ---------------------------------------------------------------
genes <- c("PTF1A", "PRDM13", "PAX2", "LHX5", "LHX9",
              "PAX5", "LHX1", "SKOR2", "ESRRB", "PCP4", 
              "SST", "LMX1A", "BARHL1", "ATOH1", "TBR1",
              "SOX14", "AIF1", "CSF1R", "CX3CR1", "ERG", "OTX2",
               "ARHGAP11B", "SMOC2", "PAX3", "CD4", "GRIN2B", "GRIN2A",
                "PAX6", "LHX2", "OLIG2", "GSX1", "WNT1", "FGF8", 
                "EOMES","PLP1", "NFIX", "PRTG", "CRYAB", "NFIA", 
                "NFIB", 'TCF7L1', 'TCF7L2', 'WLS', 'MKI67', "NKX2-1", "SHH", "XIST", "SHOX2",
                "MECOM", "SALL3", "SALL4", "LINC01965", "AL157944.1", "EFNA5", "ELP4", "GATA3",
                "SHOX", "TFAP2D", "ONECUT1", "DSCAM", "ZFHX3", "MIR137HG", "BNC2", "SOX2", "RBBP8", "CALB1", "HES5",
                "TTYH1", "SOX5", "NOTCH1", "PRDM5", "PVT1", "PARD3B", "PIEZO2", "KIRREL2", "ZIC5", "HES3", "MIR9-3HG", "ZIC4",
                "WNT3A", "SOX9", "VIM", 'PDGFRA', 'PECAM1', 'COL1A1', 'CD19', 'MS4A1', "SOX10", "PAX7", "SNAI1", "NEUROG1", "NEUROG2", 
                "NEUROG3", "EBF1", "EBF2", "EN2", "SHROOM3", "CNPY1", "IRF2BPL", "HES1", "ZIC1", "ZIC3", 
                "PVALB") %>% unique() %>% sort()


# Functions ---------------------------------------------------------------

# quick minmax scaling function
minmax_scale <- function(x, na_rm = TRUE){ 
  return((x - min(x, na.rm = na_rm)) / (max(x, na.rm = na_rm) - min(x, na.rm = na_rm)))
}


# featureplots
fp <- function(so, label, dir = "tmp", reduction = "wnn.umap.recip"){
  plot_list <- purrr::map(sort(genes), function(x){
    return(FeaturePlot(so, x, reduction = reduction, order = T, raster = T)+NoAxes()+NoLegend())
  })
  pdf(glue::glue("{dir}/featureplots_{label}.pdf"), w = 15, h = (length(genes)/2))
  print(cowplot::plot_grid(plotlist = plot_list, ncol = 6))
  dev.off()
}