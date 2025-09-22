import cv2
import numpy as np
import pandas as pd
import scipy.io
import anndata as ad
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

adata_xenium = ad.read_h5ad('path/to/file/adata.h5ad')
cell_boundaries_csv = pd.read_csv('path/to/file/cell_boundaries.csv.gz')
#in full-resolution (series 1) Xenium image, the pixel size is 0.2125 microns. Transform the vertices positions to pixel units
cell_boundaries_csv['vertex_x_coord'] = cell_boundaries_csv['vertex_x']/0.2125
cell_boundaries_csv['vertex_y_coord'] = cell_boundaries_csv['vertex_y']/0.2125
#register gcamp iamges with Xenium using Xenium explorer, get the affine matrix
affine_matrix = np.loadtxt('path/to/file/matrix.csv', delimiter=',')
affine_matrix_inv = np.linalg.inv(affine_matrix)
homogeneous_coordinates = np.hstack([cell_boundaries_csv[['vertex_x_coord', 'vertex_y_coord']].values, 
                                     np.ones((cell_boundaries_csv.shape[0], 1))])
transformed_coordinates = homogeneous_coordinates.dot(affine_matrix_inv.T)
cell_boundaries_csv['x_transformed'] = transformed_coordinates[:, 0]
cell_boundaries_csv['y_transformed'] = transformed_coordinates[:, 1]
cell_boundaries_csv['x_transformed'] = round(cell_boundaries_csv['x_transformed']).astype(int)
cell_boundaries_csv['y_transformed'] = round(cell_boundaries_csv['y_transformed']).astype(int)

#get the size of the stitched image using an example image
gcamp_image = np.load('stitched_images_gcamp/Z001.npy')
cell_boundaries_csv['x_transformed'] = cell_boundaries_csv['x_transformed'].clip(lower=0, upper=gcamp_image.shape[1]-1)
cell_boundaries_csv['y_transformed'] = cell_boundaries_csv['y_transformed'].clip(lower=0, upper=gcamp_image.shape[0]-1)

#load the stitched images
image_dir = 'stitched_images_gcamp/'
image_names = os.listdir(image_dir)
image_names.sort()
#calculate the gcamp signal mean intensity of each cell within the cell boundary
mean_intensities_ls = []
for name in image_names:
    name_path = image_dir + name
    gcamp_image = np.load(name_path)
    cell_ids = []
    mean_intensities = []
    for cell_id, group in cell_boundaries_csv.groupby('cell_id'):
        cell_ids.append(cell_id)
        x_coords = group['x_transformed'].values
        y_coords = group['y_transformed'].values
        
        rr, cc = polygon(y_coords, x_coords, gcamp_image.shape)
        mean_intensity = gcamp_image[rr, cc].mean()
        mean_intensities.append(mean_intensity)
    mean_intensities_ls.append(mean_intensities)
mean_intensities_np = np.array(mean_intensities_ls)
np.save('mean_intensities_SMAD.npy', mean_intensities_np)
np.save("cell_ids_SMAD.npy", cell_ids)


