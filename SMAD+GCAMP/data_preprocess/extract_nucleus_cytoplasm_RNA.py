import numpy as np
import pandas as pd
import anndata as ad
from shapely.geometry import Point, Polygon
import geopandas as gpd
import scanpy as sc
from scipy.sparse import csr_matrix

transcripts_path = "path/to/file/transcripts.parquet"
df_transcripts = pd.read_parquet(transcripts_path)

#filter the nucleus and cytoplasm transcripts
filtered_transcripts_nucleus = df_transcripts[
    (df_transcripts['cell_id'] != 'UNASSIGNED') & (df_transcripts['is_gene'] == True)
    & (df_transcripts['overlaps_nucleus'] == 1) & (df_transcripts['qv'] > 20)
]
filtered_transcripts_cyto = df_transcripts[
    (df_transcripts['cell_id'] != 'UNASSIGNED') & (df_transcripts['is_gene'] == True)
    & (df_transcripts['overlaps_nucleus'] == 0) & (df_transcripts['qv'] > 20)
]

#extract nucleus transcripts in cells
counts_nucleus = (
    filtered_transcripts_nucleus.groupby(["cell_id", "feature_name"])
    .size()
    .reset_index(name="count")
)
count_matrix_nucleus = counts_nucleus.pivot(index="cell_id", columns="feature_name", values="count").fillna(0)
X = csr_matrix(count_matrix_nucleus.values)
obs = pd.DataFrame(index=count_matrix_nucleus.index)
var = pd.DataFrame(index=count_matrix_nucleus.columns)
adata_nucleus = ad.AnnData(X=X, obs=obs, var=var)
adata.write_h5ad("nucleus_RNA_counts_SMAD.h5ad")

#extract cytoplasm transcripts in cells
counts_cyto = (
    filtered_transcripts_cyto.groupby(["cell_id", "feature_name"])
    .size()
    .reset_index(name="count")
)
count_matrix_cyto = counts_cyto.pivot(index="cell_id", columns="feature_name", values="count").fillna(0)
X = csr_matrix(count_matrix_cyto.values)
obs = pd.DataFrame(index=count_matrix_cyto.index)
var = pd.DataFrame(index=count_matrix_cyto.columns)
adata_cyto = ad.AnnData(X=X, obs=obs, var=var)
adata.write_h5ad("cytoplasm_RNA_counts_SMAD.h5ad")
