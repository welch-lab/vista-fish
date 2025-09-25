from os import path
import numpy as np
import pandas as pd
import m2stitch
import os
import cv2
import tifffile as tiff
import imageio
import matplotlib.pyplot as plt
from skimage.transform import resize
from skimage import exposure

image_dir = 'path/to/IF_staining/images/'
image_names = os.listdir(image_dir)
image_names.sort()
file_names = [f"path/to/IF_staining/images/{s}" for s in image_names]

overlap_rate = 0.1 #ovelap rate of IF staining images
col_num = 12 #number of columns of for image stitching
row_num = 22 #number of rows of for image stitching
image_size = 2000 #size of each image, for this example, each image is in 2000*2000 pixels

col = list(range(col_num))*row_num
row = [i for i in range(row_num) for _ in range(col_num)]
stitched_props = {
    'row': row,
    'col': col
}
result_df = pd.DataFrame(stitched_props)
x_pos = []
y_pos = []
for i in range(col_num*row_num):
    x_pos.append(int(col[i]*image_size*(1-overlap_rate)))
    y_pos.append(int(row[i]*image_size*(1-overlap_rate)))
result_df['y_pos'] = y_pos
result_df['x_pos'] = x_pos
result_df["y_pos2"] = result_df["y_pos"] - result_df["y_pos"].min()
result_df["x_pos2"] = result_df["x_pos"] - result_df["x_pos"].min()

#Axon staining image
tif_files_W1_C1 = [path for path in file_names if "W0001" in path and "C1" in path]
tif_files_W1_C1.sort()
images = []
for tiff_file_path in tif_files_W1_C1:
    image = imageio.imread(tiff_file_path)
    images.append(image)
images = np.stack(images)
size_y = images.shape[1]
size_x = images.shape[2]
stitched_image_size = (
    result_df["y_pos2"].max() + size_y,
    result_df["x_pos2"].max() + size_x,
)
stitched_image = np.zeros_like(images, shape=stitched_image_size)
for i, row in result_df.iterrows():
    stitched_image[
        row["y_pos2"] : row["y_pos2"] + size_y,
        row["x_pos2"] : row["x_pos2"] + size_x,
    ] = images[i]
result_image_file_path = 'stitched_images_IF_staining/' + 'W1_C1' + '.npy'
tiff.imwrite(result_image_file_path, stitched_image)

#Dendrite staining image
tif_files_W1_C2 = [path for path in file_names if "W0001" in path and "C2" in path]
tif_files_W1_C2.sort()
images = []
for tiff_file_path in tif_files_W1_C2:
    image = imageio.imread(tiff_file_path)
    images.append(image)
images = np.stack(images)
size_y = images.shape[1]
size_x = images.shape[2]
stitched_image_size = (
    result_df["y_pos2"].max() + size_y,
    result_df["x_pos2"].max() + size_x,
)
stitched_image = np.zeros_like(images, shape=stitched_image_size)
for i, row in result_df.iterrows():
    stitched_image[
        row["y_pos2"] : row["y_pos2"] + size_y,
        row["x_pos2"] : row["x_pos2"] + size_x,
    ] = images[i]
result_image_file_path = 'stitched_images_IF_staining/' + 'W1_C2' + '.npy'
tiff.imwrite(result_image_file_path, stitched_image)

