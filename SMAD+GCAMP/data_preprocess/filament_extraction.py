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

#load filament image
img = io.imread('path/to/file/morphology_focus/morphology_focus_0000.ome.tif')
print(img.shape)
filament_image = img[3,:,:].copy()

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
filtered_transcripts['y'] = round(filtered_transcripts['y_coord']).astype(int)

#extract filament transcripts
threshold=1000
rows, cols = np.where(filament_image > threshold)
coords_xy = np.column_stack((cols,rows))
coords_df = pd.DataFrame(coords_xy[:,:2], columns=['x', 'y'])
filament_transcripts = pd.merge(filtered_transcripts, coords_df, on=['x', 'y'], how='inner')

#run DEG
adata = sc.read_10x_h5(filename="path/to/file/cell_feature_matrix.h5")
filament_grouped = filament_transcripts.groupby("feature_name").size().reset_index(name="filament_count")
cell_grouped = pd.DataFrame(np.array(adata.X.sum(axis=0)).flatten(), index=adata.var.index, columns=["cell_count"])
cell_grouped.reset_index(inplace=True)
cell_grouped.rename(columns={"index": "feature_name"}, inplace=True)
merged = pd.merge(filament_grouped, cell_grouped, on="feature_name", how="outer").fillna(0)
merged["normalized_filament_count"] = (merged["filament_count"] / merged["filament_count"].sum()) * 1e6
merged["normalized_cell_count"] = (merged["cell_count"] / merged["cell_count"].sum()) * 1e6
merged["p_value"] = merged.apply(
    lambda row: ttest_ind(
        np.repeat(row["normalized_filament_count"], int(row["filament_count"])),
        np.repeat(row["normalized_cell_count"], int(row["cell_count"])),
        equal_var=False
    )[1] if row["filament_count"] > 1 and row["cell_count"] > 1 else 1.0,
    axis=1
)
merged["adjusted_p_value"] = merged["p_value"] * len(merged)
merged["adjusted_p_value"] = merged["adjusted_p_value"].clip(upper=1)
merged["log2_fold_change"] = np.log2((merged["normalized_filament_count"] + 1) / (merged["normalized_cell_count"] + 1))
fc_threshold = 1
p_threshold = 0.05
merged["expression_category"] = merged.apply(
    lambda row: "Highly expressed in filaments"
    if row["log2_fold_change"] > fc_threshold and row["adjusted_p_value"] < p_threshold
    else "Highly expressed in cells"
    if row["log2_fold_change"] < -fc_threshold and row["adjusted_p_value"] < p_threshold
    else "Not significant",
    axis=1,
)
merged.to_csv("filament_cell_DEGs_results.csv", index=False)

