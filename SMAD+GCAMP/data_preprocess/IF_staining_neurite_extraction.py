import cv2
import numpy as np
import pandas as pd
import scipy.io
import anndata
from skimage.transform import warp
from skimage.draw import polygon
from PIL import Image
import tifffile as tiff
import seaborn as sns
import matplotlib.pyplot as plt
import dask.array as da
import imageio
import imagecodecs
import math
import skimage
import os
from skimage import transform
import scanpy as sc
from scipy.stats import ttest_ind

#transform the cell boundry using affine matrix
affine_matrix = np.loadtxt('path/to/file/matrix.csv', delimiter=',')
affine_matrix_inv = np.linalg.inv(affine_matrix)
csv_file_path = 'path/to/file/cell_boundaries.csv.gz'
cell_boundaries_csv = pd.read_csv(csv_file_path)
cell_boundaries_csv['vertex_x_coord'] = cell_boundaries_csv['vertex_x']/0.2125
cell_boundaries_csv['vertex_y_coord'] = cell_boundaries_csv['vertex_y']/0.2125
homogeneous_coordinates = np.hstack([cell_boundaries_csv[['vertex_x_coord', 'vertex_y_coord']].values, 
                                     np.ones((cell_boundaries_csv.shape[0], 1))])
transformed_coordinates = homogeneous_coordinates.dot(affine_matrix_inv.T)
cell_boundaries_csv['x_transformed'] = transformed_coordinates[:, 0]
cell_boundaries_csv['y_transformed'] = transformed_coordinates[:, 1]
cell_boundaries_csv['x_transformed'] = round(cell_boundaries_csv['x_transformed']).astype(int)
cell_boundaries_csv['y_transformed'] = round(cell_boundaries_csv['y_transformed']).astype(int)

#load RNA expression data
adata = sc.read_10x_h5(
    filename="path/to/file/cell_feature_matrix.h5"
)

#load transcripts
transcripts_path = "path/to/file/transcripts.parquet"
df_transcripts = pd.read_parquet(transcripts_path)
filtered_transcripts = df_transcripts[
    (df_transcripts['cell_id'] == 'UNASSIGNED') & (df_transcripts['is_gene'] == True)
    & (df_transcripts['qv'] > 20)
]
filtered_transcripts['x_coord'] = filtered_transcripts['x_location']/0.2125
filtered_transcripts['y_coord'] = filtered_transcripts['y_location']/0.2125
filtered_transcripts['x'] = round(filtered_transcripts['x_coord']).astype(int)
filtered_transcripts['y'] = round(filtered_transcripts['x_coord']).astype(int)

#get the high intensity pixels in Axon staining
w1c2_image = np.load('path/to/axon/image/W1_C2.npy')
cell_dict = cell_boundaries_csv.groupby("cell_id")[["x_transformed", "y_transformed"]].apply(lambda x: x.values.tolist())
for cell_id, vertices in cell_dict.items():
    polygon = np.array(vertices, dtype=np.int32)
    cv2.fillPoly(w1c2_image, [polygon], 0)
threshold=110
rows, cols = np.where(w1c2_image > threshold)
coords_xy = np.column_stack((cols,rows))
coords_df = pd.DataFrame(coords_xy[:,:2], columns=['x_coord', 'y_coord'])
homogeneous_coordinates = np.hstack([coords_df[['x_coord', 'y_coord']].values, np.ones((coords_df.shape[0], 1))])
transformed_coordinates = homogeneous_coordinates.dot(affine_matrix.T)
coords_df['x_transformed'] = transformed_coordinates[:, 0]
coords_df['y_transformed'] = transformed_coordinates[:, 1]
coords_df['x'] = round(coords_df['x_transformed']).astype(int)
coords_df['y'] = round(coords_df['y_transformed']).astype(int)
axon_transcripts = pd.merge(filtered_transcripts, coords_df, on=['x', 'y'], how='inner')

#run DEG for Axon staining
axon_grouped = axon_transcripts.groupby("feature_name").size().reset_index(name="axon_count")
cell_grouped = pd.DataFrame(np.array(adata.X.sum(axis=0)).flatten(), index=adata.var.index, columns=["cell_count"])
cell_grouped.reset_index(inplace=True)
cell_grouped.rename(columns={"index": "feature_name"}, inplace=True)
merged = pd.merge(axon_grouped, cell_grouped, on="feature_name", how="outer").fillna(0)
merged["normalized_axon_count"] = (merged["axon_count"] / merged["axon_count"].sum()) * 1e6
merged["normalized_cell_count"] = (merged["cell_count"] / merged["cell_count"].sum()) * 1e6
merged["p_value"] = merged.apply(
    lambda row: ttest_ind(
        np.repeat(row["normalized_axon_count"], int(row["axon_count"])),
        np.repeat(row["normalized_cell_count"], int(row["cell_count"])),
        equal_var=False
    )[1] if row["axon_count"] > 1 and row["cell_count"] > 1 else 1.0,
    axis=1
)
merged["adjusted_p_value"] = merged["p_value"] * len(merged)
merged["adjusted_p_value"] = merged["adjusted_p_value"].clip(upper=1)
merged["log2_fold_change"] = np.log2((merged["normalized_axon_count"] + 1) / (merged["normalized_cell_count"] + 1))
fc_threshold = 1
p_threshold = 0.05
merged["expression_category"] = merged.apply(
    lambda row: "Highly expressed in axons"
    if row["log2_fold_change"] > fc_threshold and row["adjusted_p_value"] < p_threshold
    else "Highly expressed in cells"
    if row["log2_fold_change"] < -fc_threshold and row["adjusted_p_value"] < p_threshold
    else "Not significant",
    axis=1,
)
print(merged["expression_category"].value_counts())
merged.to_csv("axon_cell_DEGs_results.csv", index=False)

#get the high intensity pixels in Dendrite staining
w1c1_image = np.load('path/to/file/W1_C1.npy')
for cell_id, vertices in cell_dict.items():
    polygon = np.array(vertices, dtype=np.int32)
    cv2.fillPoly(w1c1_image, [polygon], 0)
threshold=200
rows, cols = np.where(w1c1_image > threshold)
coords_xy = np.column_stack((cols,rows))
coords_df = pd.DataFrame(coords_xy[:,:2], columns=['x_coord', 'y_coord'])
homogeneous_coordinates = np.hstack([coords_df[['x_coord', 'y_coord']].values, np.ones((coords_df.shape[0], 1))])
transformed_coordinates = homogeneous_coordinates.dot(affine_matrix.T)
coords_df['x_transformed'] = transformed_coordinates[:, 0]
coords_df['y_transformed'] = transformed_coordinates[:, 1]
coords_df['x'] = round(coords_df['x_transformed']).astype(int)
coords_df['y'] = round(coords_df['y_transformed']).astype(int)
dendrites_transcripts = pd.merge(filtered_transcripts, coords_df, on=['x', 'y'], how='inner')

#run DEG for Dendrite staining
dendrite_grouped = dendrites_transcripts.groupby("feature_name").size().reset_index(name="dendrite_count")
cell_grouped = pd.DataFrame(np.array(adata.X.sum(axis=0)).flatten(), index=adata.var.index, columns=["cell_count"])
cell_grouped.reset_index(inplace=True)
cell_grouped.rename(columns={"index": "feature_name"}, inplace=True)
merged = pd.merge(dendrite_grouped, cell_grouped, on="feature_name", how="outer").fillna(0)
merged["normalized_dendrite_count"] = (merged["dendrite_count"] / merged["dendrite_count"].sum()) * 1e6
merged["normalized_cell_count"] = (merged["cell_count"] / merged["cell_count"].sum()) * 1e6
merged["p_value"] = merged.apply(
    lambda row: ttest_ind(
        np.repeat(row["normalized_dendrite_count"], int(row["dendrite_count"])),
        np.repeat(row["normalized_cell_count"], int(row["cell_count"])),
        equal_var=False
    )[1] if row["dendrite_count"] > 1 and row["cell_count"] > 1 else 1.0,
    axis=1
)
merged["adjusted_p_value"] = merged["p_value"] * len(merged)
merged["adjusted_p_value"] = merged["adjusted_p_value"].clip(upper=1)
merged["log2_fold_change"] = np.log2((merged["normalized_dendrite_count"] + 1) / (merged["normalized_cell_count"] + 1))
fc_threshold = 1
p_threshold = 0.05
merged["expression_category"] = merged.apply(
    lambda row: "Highly expressed in dendrites"
    if row["log2_fold_change"] > fc_threshold and row["adjusted_p_value"] < p_threshold
    else "Highly expressed in cells"
    if row["log2_fold_change"] < -fc_threshold and row["adjusted_p_value"] < p_threshold
    else "Not significant",
    axis=1,
)
print(merged["expression_category"].value_counts())
merged.to_csv("dendrite_cell_DEGs_results.csv", index=False)
