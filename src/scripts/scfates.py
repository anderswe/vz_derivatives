# Info --------------------------------------------------------------------
# Anders E.
# Apr 24, 2023
# Taylor lab

# Notes --------------------------------------------------------------------
# 
# aim is to 1. generate a tree structure
# 2. assign pseudotime to each cell
# 

# Modules & settings --------------------------------------------------------------
import sys
import os
from anndata import AnnData
import numpy as np
import pandas as pd
import scanpy as sc
import scFates as scf
import palantir
import matplotlib.pyplot as plt
import seaborn

sc.settings.verbosity = 3
sc.settings.logfile = sys.stdout
seaborn.reset_orig()
sc.set_figure_params()
scf.set_figure_pubready()
outdir = "out"
os.makedirs(outdir, exist_ok=True)


# Import and clean --------------------------------------------------------------
ad_path = "out/cs20s.h5ad"
ad = sc.read_h5ad(ad_path)

# add embeddings
embs = pd.read_csv("out/metadata_cs20s.csv").set_index('sample_barcode')[['UMAP_1', 'UMAP_2']]
ad.obsm['X_umap'] = embs
ad.obsm['X_umap'].index.name = None
ad.obsm["X_umap"] = ad.obsm["X_umap"].to_numpy() # convert UMAP to numpy array

# subset out gaba ntz
non_gaba_ntz_cells = ad.obs[ad.obs['SCT_snn_res.1'] != 3].index

ad = ad[non_gaba_ntz_cells, :]

# Commented out since using pre-existing umap --------------------------------------------------------------
# norm_df = sc.pp.normalize_per_cell(np.asarray(ad.X.todense()),copy=True)
# norm_df = palantir.preprocess.log_transform(norm_df)

# adata=sc.AnnData(norm_df)
# sc.pp.highly_variable_genes(adata, n_top_genes=1500, flavor='cell_ranger')
# sc.pp.pca(adata)
# pca_projections = pd.DataFrame(adata.obsm["X_pca"],index=adata.obs_names)

# dm_res = palantir.utils.run_diffusion_maps(pca_projections)
# ms_data = palantir.utils.determine_multiscale_space(dm_res,n_eigs=4)

# # generate neighbor draph in multiscale diffusion space
# adata.obsm["X_palantir"]=ms_data.values
# sc.pp.neighbors(adata,n_neighbors=30,use_rep="X_palantir")


# sc.set_figure_params()
# sc.pl.draw_graph(adata,color="CD34",color_map="RdBu_r")

# Run --------------------------------------------------------------

# get sigma
sig = scf.tl.explore_sigma(ad,Nodes=50,use_rep="X_umap",seed=42,plot=True)
plt.tight_layout()
plt.savefig((outdir + '/sigma.png'), dpi = 300)
plt.clf()

# fit tree
scf.tl.tree(ad,
            Nodes=50,
            use_rep="X_umap",
            method="ppt",
            ppt_nsteps=200,
            ppt_sigma=sig,
            ppt_lambda=100,
            seed=42,plot=True)

# plot tree over embedding
scf.pl.graph(ad)
plt.tight_layout()
plt.savefig((outdir + '/tree.png'), dpi = 300)
plt.clf()

# assign root
root = 17
scf.tl.root(ad,root)

# calculate pseudotime
scf.tl.pseudotime(ad,n_jobs=10,n_map=100,seed=666)

# plot pseudotime
sc.set_figure_params(frameon=False,dpi_save=300)
scf.pl.trajectory(ad) #, color_cells = 'clusters', palette = colors_dutch
plt.tight_layout()
plt.savefig((outdir + '/tree_pseudotime_root_' + str(root) + '.png'))
plt.clf()

# save pseudotime
ad.obs.to_csv(outdir + "/metadata_scfates.csv")

# dendrogram
scf.tl.dendrogram(ad)
sc.set_figure_params(frameon=False,dpi_save=300)
scf.pl.dendrogram(ad,color="t",show_info=False,cmap="viridis",title = " ")
plt.tight_layout()
plt.savefig((outdir + '/dendrogram.png'), dpi = 300)
plt.clf()

# test gene associations over pseudotime
scf.tl.test_association(ad, n_jobs=4, fdr_cut = 0.0001, A_cut = 1) # defaults: A_cut = 1, fdr_cut = 0.0001
sc.set_figure_params()
scf.pl.test_association(ad)
plt.tight_layout()
plt.savefig((outdir + '/gene_associations_root.png'), dpi = 300)
plt.clf()

# save before fitting
ad.write((outdir + "/ad_scfates_pre_fit.h5ad"))
# ad = sc.read_h5ad((outdir + "/ad_scfates.h5ad"))


# fit significant genes --> subsets to signif only
scf.tl.fit(ad, n_jobs=10, gamma=1.5)
ad.var.to_csv(outdir + "/lineage_drivers_scfates.csv")

# save after fitting
ad.write((outdir + "/ad_scfates_post_fit.h5ad"))
# ad = sc.read_h5ad((outdir + "/ad_scfates_post_fit.h5ad"))

# test DE between branches
scf.tl.test_fork(ad,root_milestone=str(root),milestones=["8", "24"],n_jobs=10, rescale=True)


# assign branches to DE genes
scf.tl.branch_specific(ad,root_milestone=str(root),milestones=["8", "24"], effect=0.001) # basically as lenient as possible

# save
ad.uns['17->8<>24']['fork'].rename(columns={"8": "gaba", "24": "glut"}).to_csv(outdir + "/scfates_de_branches.csv")


# Interactive session follow-up --------------------------------------------------------------

g1 = scf.pl.trends(ad,
                 root_milestone=str(root),
                 milestones=["8", "24"],
                 branch="8",
                 plot_emb=False,ordering="max",return_genes=True)
plt.tight_layout()
plt.savefig((outdir + '/gene_trends.png'), dpi = 300)
plt.clf()
