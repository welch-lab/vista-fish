import matplotlib.pyplot as plt
import seaborn as sns
import scanpy as sc
import squidpy as sq
import numpy as np
import pandas as pd
import anndata as ad
import scvelo as scv
import pyensembl
import os

sc.settings.verbosity = 3 
sc.logging.print_header()
sc.settings.set_figure_params(dpi=100, facecolor="white")

xenium_dir = "../cell_feature_matrix"
zarr_path = "../Xenium.zarr"

adata = sc.read_10x_mtx(
    xenium_dir,             
    var_names='gene_symbols', 
    cache=True            
)

adata.var

adata.obs['total_counts'] = adata.X.sum(axis=1).A1

adata_filtered = adata[adata.obs['total_counts'] >= 500].copy()
adata_filtered

adata_filtered.layers["counts"] = adata_filtered.X.copy()

df = pd.read_csv(
    "../cells.csv.gz",
index_col=0)

print(adata_filtered.shape)

print("df index example:", df.index[:5])  
print("adata_filtered.obs_names example:", adata_filtered.obs_names[:5]) 

matching_ids = set(df.index) & set(adata_filtered.obs_names)
print("Number of matching cell IDs:", len(matching_ids))


df_filtered = df.loc[df.index.isin(adata_filtered.obs_names)]
df_filtered = df_filtered.loc[adata_filtered.obs_names]  


df_filtered.set_index(adata_filtered.obs_names, inplace=True)
adata_filtered.obs = df_filtered.copy()
adata_filtered.obsm["spatial"] = adata_filtered.obs[["x_centroid", "y_centroid"]].copy().to_numpy()

adata_filtered.obs

sc.pp.calculate_qc_metrics(adata_filtered, percent_top=(10, 20, 50, 150), inplace=True)

fig, axs = plt.subplots(1, 3, figsize=(15, 4))
axs[0].set_title("Total transcripts per cell")
sns.histplot(
    adata_filtered.obs["total_counts"],
    kde=False,
    ax=axs[0],
)
axs[1].set_title("Unique transcripts per cell")
sns.histplot(
    adata_filtered.obs["n_genes_by_counts"],
    kde=False,
    ax=axs[1],
)
axs[2].set_title("Area of segmented cells")
sns.histplot(
    adata_filtered.obs["cell_area"],
    kde=False,
    ax=axs[2],
)

adata_filtered.layers["counts"] = adata_filtered.X.copy()

sc.pp.normalize_total(adata_filtered, inplace=True)
sc.pp.log1p(adata_filtered)
sc.pp.pca(adata_filtered)
adata_filtered

sc.pl.pca_variance_ratio(adata_filtered, n_pcs=50, log=True)
sc.pp.neighbors(adata_filtered)
sc.tl.umap(adata_filtered)
adata_filtered

sc.tl.leiden(adata_filtered, resolution=0.6)

sc.pl.umap(adata_filtered, color=["leiden"])

sc.pl.umap(
    adata_filtered,
    color=[
        "total_counts",
        "n_genes_by_counts",
        "leiden",
    ],
    wspace=0.4,
)

clustering = ['gene_names']
sc.pl.dotplot(adata_filtered, clustering, groupby = 'leiden')

cluster_annotations = {
    '0': 'cluster 0', 
    '1': 'Pcluster 1', 
    '2': 'cluster 2', 
    '3': 'cluster 3', 
    '4': 'cluster 4', 
    '5': 'cluster 5', 
} 

adata_filtered.obs['cell_type'] = adata_filtered.obs['leiden'].map(cluster_annotations).astype('category')

sc.pl.umap(adata_filtered, color='cell_type', legend_loc='right margin', 
           legend_fontsize=8, legend_fontoutline=2, frameon=False, title='')

gene = "IDH1"
print(gene in adata_filtered.var_names)
sq.pl.spatial_scatter(
    adata_filtered,
    library_id = "spatial",
    shape = None,
    color = [gene],
    wspace= 0.4,
)







