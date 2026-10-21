#!/usr/bin/env python
# coding: utf-8

# In[1]:


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


# In[2]:


sc.settings.verbosity = 3  # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.logging.print_header()
sc.settings.set_figure_params(dpi=100, facecolor="white")


# In[3]:


xenium_dir = "../cell_feature_matrix"
zarr_path = "../Xenium.zarr"


# In[4]:


adata = sc.read_10x_mtx(
    xenium_dir,             # the folder with matrix.mtx, barcodes.tsv, features.tsv
    var_names='gene_symbols',  # use gene names instead of IDs (can also use 'gene_ids')
    cache=True            # optional: speeds up future loading
)


# In[5]:


adata.var


# In[6]:


adata.obs['total_counts'] = adata.X.sum(axis=1).A1


# In[7]:


adata_filtered = adata[adata.obs['total_counts'] >= 500].copy()
adata_filtered


# In[8]:


adata_filtered.layers["counts"] = adata_filtered.X.copy()


# In[9]:


df = pd.read_csv(
    "../cells.csv.gz",
index_col=0)


# In[10]:


df


# In[11]:


print(adata_filtered.shape)


# In[12]:


print("df index example:", df.index[:5])  
print("adata_filtered.obs_names example:", adata_filtered.obs_names[:5]) 


# In[13]:


matching_ids = set(df.index) & set(adata_filtered.obs_names)
print("Number of matching cell IDs:", len(matching_ids))


# In[14]:


df_filtered = df.loc[df.index.isin(adata_filtered.obs_names)]
df_filtered = df_filtered.loc[adata_filtered.obs_names]  


# In[15]:


df_filtered.set_index(adata_filtered.obs_names, inplace=True)
adata_filtered.obs = df_filtered.copy()
adata_filtered.obsm["spatial"] = adata_filtered.obs[["x_centroid", "y_centroid"]].copy().to_numpy()


# In[16]:


adata_filtered.obs


# In[17]:


sc.pp.calculate_qc_metrics(adata_filtered, percent_top=(10, 20, 50, 150), inplace=True)


# In[18]:


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


# In[19]:


adata_filtered.layers["counts"] = adata_filtered.X.copy()


# In[20]:


sc.pp.normalize_total(adata_filtered, inplace=True)


# In[21]:


sc.pp.log1p(adata_filtered)


# In[22]:


sc.pp.pca(adata_filtered)


# In[23]:


adata_filtered


# In[24]:


sc.pl.pca_variance_ratio(adata_filtered, n_pcs=50, log=True)


# In[25]:


sc.pp.neighbors(adata_filtered)


# In[26]:


sc.tl.umap(adata_filtered)


# In[27]:


adata_filtered


# In[42]:


sc.tl.leiden(adata_filtered, resolution=0.6)


# In[43]:


sc.pl.umap(adata_filtered, color=["leiden"])


# In[44]:


sc.pl.umap(
    adata_filtered,
    color=[
        "total_counts",
        "n_genes_by_counts",
        "leiden",
    ],
    wspace=0.4,
)


# In[45]:


clustering = ['POU4F1','GPM6A','PHOX2B','PRPH','VIM','CD99','NTNG2']
sc.pl.dotplot(adata_filtered, clustering, groupby = 'leiden')

cluster_annotations = {
    '0': 'POU4F1 Arsenite-treated', 
    '1': 'POU4F1/NTNG2', 
    '2': 'CNTNAP3B/POU4F1/GPM6A', 
    '3': 'PHOX2B/PRPH', 
    '4': 'CD99/VIM Arsenite-treated', 
    '5': 'CD99/VIM', 
} 

adata_filtered.obs['cell_type'] = adata_filtered.obs['leiden'].map(cluster_annotations).astype('category')

sc.pl.umap(adata_filtered, color='cell_type', legend_loc='right margin', 
           legend_fontsize=8, legend_fontoutline=2, frameon=False, title='')


# In[46]:


gene = "IDH1"

print(gene in adata_filtered.var_names)

sq.pl.spatial_scatter(
    adata_filtered,
    library_id = "spatial",
    shape = None,
    color = [gene],
    wspace= 0.4,
)







