import matplotlib.pyplot as plt
import seaborn as sns
import scanpy as sc
import squidpy as sq
import numpy as np
import pandas as pd
import anndata as ad
import scvelo as scv
import pyensembl
from scipy.io import mmread
import itertools
import os

sc.settings.set_figure_params(dpi=300)

data_dir = "/Users/javidghaemmaghami/Desktop/Welch_Lab/Xenium Project/gRNA Detection_round2/data_dir"
adata = sc.read_10x_mtx(
    data_dir,
    var_names='gene_symbols',
    cache=True
)

adata.obs['total_counts'] = adata.X.sum(axis=1).A1
num_cells_800 = (adata.obs['total_counts'] >= 800).sum()
num_cells_800
adata_filtered = adata[adata.obs['total_counts'] >= 800].copy()
adata_filtered.layers["counts"] = adata_filtered.X.copy()
sc.pp.normalize_total(adata_filtered)
sc.pp.log1p(adata_filtered)
sc.pp.scale(adata_filtered, max_value=10)
sc.tl.pca(adata_filtered)
sc.pl.pca_variance_ratio(adata_filtered, n_pcs=50, log=True)
sc.pp.neighbors(adata_filtered,n_neighbors=5,n_pcs=30, metric='cosine')
sc.tl.umap(adata_filtered)
sc.tl.leiden(adata_filtered, flavor="igraph")

cluster_rename_map = {
    '0': 'GPM6A',
    '1': 'GPM6A', 
    '2': 'GPM6A',    
    '3': 'GPM6A',
    '4': 'GPM6A',
    '5': 'POU4F1',
    '6': 'POU4F1',  
    '7': 'POU4F1',
    '8': 'CD99',
    '9': 'POU4F1',
    '10': 'GPM6A',
    '11': 'CD99',
    '12': 'PHOX2B',
    '13': 'CD99',
    '14': 'POU4F1',
    '15': 'CD99',
    '16': 'PHOX2B',
    '17': 'POU4F1',
    '18': 'CD99',
    '19': 'PHOX2B',
}

adata_filtered.obs['leiden_renamed_2'] = adata_filtered.obs['leiden'].replace(cluster_rename_map)
sc.pl.umap(adata_filtered, color='leiden_renamed_2')

# Guides
sgRNA_genes = [gene for gene in adata_filtered.var_names if '-' in gene]
sgRNA_to_target = {sg: sg.split('-')[0] for sg in sgRNA_genes}
sgRNA_expr = adata_filtered[:, sgRNA_genes].X
if not isinstance(sgRNA_expr, np.ndarray):
    sgRNA_expr = sgRNA_expr.toarray()
sgRNA_df = pd.DataFrame(
    sgRNA_expr,
    index=adata_filtered.obs_names,
    columns=sgRNA_genes
)
target_df = sgRNA_df.groupby(sgRNA_to_target, axis=1).sum()
target_binary = (target_df > 0).astype(int)
for target in target_binary.columns:
    adata_filtered.obs[target] = target_binary[target]
gRNA_cols = [col for col in adata_filtered.obs.columns if col not in ['total_counts', 'Negative_control'] and 
             np.isin(adata_filtered.obs[col].unique(), [0, 1]).all()]
unperturbed_mask = (adata_filtered.obs[gRNA_cols].sum(axis=1) == 0)
adata_filtered.obs['unperturbed'] = unperturbed_mask.astype(str)


# Plotting Guide Probes
grnas = [
    'AFF1','BCL6','BNC2','CARHSP1','CREB5','EBF1','EBF3','HMGB2',
    'ID2','ID4','ISL1','JUNB','KLF3','LHX9','MEIS2','NFE2L2','NFIB',
    'ONECUT1','ONECUT2','PROX1','RUNX1','SIX1','SOX5','ZBTB7C'
]
expected_probes = [f"{g}-{i}" for g in grnas for i in range(1,6)]

probes = [p for p in expected_probes if p in adata_filtered.var_names]

if not probes:
    raise ValueError("No probes found in var_names! Check naming scheme.")
mat = adata_filtered[:, probes].X

if sparse.issparse(mat):
    mat = mat.toarray()
counts = np.sum(mat > 0, axis=0)

probe_counts = pd.Series(counts, index=probes).sort_values(ascending=False)

plt.figure(figsize=(12, 6))
probe_counts.plot(kind='bar')
plt.ylabel('Number of cells')
plt.title('Cells per probe (nonzero counts)')
plt.xticks(rotation=90, ha='right')
plt.tight_layout()
plt.show()

# Removing outlier probes
probes_to_remove = ['ID2-4', 'SOX5-3', 'NFIB-2', 'ID2-2']
existing = [p for p in probes_to_remove if p in adata_filtered.var_names]
missing  = set(probes_to_remove) - set(existing)
if missing:
    print(f"These probes were not found/skipping: {missing}")
X = adata_filtered[:, existing].X  # shape = (n_cells, n_probes)

if sparse.issparse(X):
    mat = X.toarray()
else:
    mat = np.asarray(X)
bool_mat = mat > 0
counts_per_cell = bool_mat.sum(axis=1)  
mask_exclusive = counts_per_cell == 1
adata = adata_filtered[~mask_exclusive].copy()

# Plotting cells on UMAP
xenium_cell_ids = ['boljjmgm-1','bolicola-1', 'chhohmfi-1', 'bolonabg-1']   
cells_found = [cell for cell in xenium_cell_ids if cell in adata_filtered.obs_names]

if len(cells_found) == 0:
else:
    subset_adata = adata_filtered[cells_found, :]
    print(subset_adata.obsm['X_umap'])
    
    plt.figure(figsize=(8, 6))
    sc.pl.umap(adata_filtered, show=False, color='leiden_renamed_2',  title='Cell Type Cluster UMAP')
    
    selected_umap = subset_adata.obsm['X_umap']
    plt.scatter(selected_umap[:, 0], selected_umap[:, 1],
                s=4, c='red', edgecolors='black')
    plt.show()
