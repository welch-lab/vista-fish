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

image_dir = 'path/to/perturbation_BFP/images/'
image_names = os.listdir(image_dir)
image_names.sort()
file_names = [f"path/to/perturbation_BFP/images/{s}" for s in image_names]

#calculate the median intensity of each pixel position
images = np.array([imageio.imread(file) for file in file_names])
median_intensity_image = np.median(images, axis=0)

#subtract pixel medians and clip for all images
for image_name in image_names:
    file_name = 'path/to/perturbation_BFP/images/' + image_name
    image = imageio.imread(file_name)
    image = np.clip(np.subtract(image, median_image), 0, 500)
    result_image_file_path = 'median_subtracted/' + image_name + '.npy'
    np.save(result_image_file_path, image)
    
#modify the overlap regions of median subtracted images
image_dir = 'median_subtracted/'
image_names = os.listdir(image_dir)
image_names.sort()

overlap_rate = 0.02 #ovelap rate of BFP images
col_num = 15 #number of columns of for image stitching
row_num = 32 #number of rows of for image stitching
image_size = 2000 #size of each image, for this example, each image is in 2000*2000 pixels

reshaped_images = []
for name in image_names:
    tiff_file_path = image_dir + name
    image = np.load(tiff_file_path)
    reshaped_images.append(image)
images = np.stack(reshaped_images)

for i in range(images.shape[0]):
    image_name = image_names_sub[i]
    image = images[i].copy()
    new_image = images[i].copy()
    #right overlap
    if (i+1)%col_num != 0:
        image_right = images[i+1].copy()
        overlap_region = [image[:,image_size*(1-overlap_rate):],image_right[:,:image_size*overlap_rate]]
        stacked_overlap_region = np.stack(overlap_region)
        median_overlap_region = np.median(stacked_overlap_region, axis=0)
        new_image[:,image_size*(1-overlap_rate):] = np.ceil(median_overlap_region).astype(int)
    #left overlap
    if i%col_num != 0:
        image_left = images[i-1].copy()
        overlap_region = [image[:,:image_size*overlap_rate],image_left[:,image_size*(1-overlap_rate):]]
        stacked_overlap_region = np.stack(overlap_region)
        median_overlap_region = np.median(stacked_overlap_region, axis=0)
        new_image[:,:image_size*overlap_rate] = np.ceil(median_overlap_region).astype(int)
    #down overlap
    if i < col_num*(row_num-1):
        image_down = images[i+col_num].copy()
        overlap_region = [image[image_size*(1-overlap_rate):,:],image_down[:image_size*overlap_rate,:]]
        stacked_overlap_region = np.stack(overlap_region)
        median_overlap_region = np.median(stacked_overlap_region, axis=0)
        new_image[image_size*(1-overlap_rate):,:] = np.ceil(median_overlap_region).astype(int)
    #top overlap
    if i >= col_num:
        image_top = images[i-col_num].copy()
        overlap_region = [image[:image_size*overlap_rate,:],image_top[image_size*(1-overlap_rate):,:]]
        stacked_overlap_region = np.stack(overlap_region)
        median_overlap_region = np.median(stacked_overlap_region, axis=0)
        new_image[:image_size*overlap_rate,:] = np.ceil(median_overlap_region).astype(int)
    #down right
    if i < col_num*(row_num-1) and (i+1)%col_num != 0:
        image_dr = images[i+col_num+1].copy()
        overlap_region = [image[image_size*(1-overlap_rate):,image_size*(1-overlap_rate):],image_right[image_size*(1-overlap_rate):,:image_size*overlap_rate],image_down[:image_size*overlap_rate,image_size*(1-overlap_rate):],
                         image_dr[:image_size*overlap_rate,:image_size*overlap_rate]]
        stacked_overlap_region = np.stack(overlap_region)
        median_overlap_region = np.median(stacked_overlap_region, axis=0)
        new_image[image_size*(1-overlap_rate):,image_size*(1-overlap_rate):] = np.ceil(median_overlap_region).astype(int)
    #down left
    if i < col_num*(row_num-1) and i%col_num != 0:
        image_dl = images[i+col_num-1].copy()
        overlap_region = [image[image_size*(1-overlap_rate):,:image_size*overlap_rate],image_left[image_size*(1-overlap_rate):,image_size*(1-overlap_rate):],image_down[:image_size*overlap_rate,:image_size*overlap_rate],
                         image_dl[:image_size*overlap_rate,image_size*(1-overlap_rate):]]
        stacked_overlap_region = np.stack(overlap_region)
        median_overlap_region = np.median(stacked_overlap_region, axis=0)
        new_image[image_size*(1-overlap_rate):,:image_size*overlap_rate] = np.ceil(median_overlap_region).astype(int)
    #top right
    if i >= col_num and (i+1)%col_num != 0:
        image_tr = images[i-col_num+1].copy()
        overlap_region = [image[:image_size*overlap_rate,image_size*(1-overlap_rate):],image_right[:image_size*overlap_rate,:image_size*overlap_rate],image_top[image_size*(1-overlap_rate):,image_size*(1-overlap_rate):],
                         image_tr[image_size*(1-overlap_rate):,:image_size*overlap_rate]]
        stacked_overlap_region = np.stack(overlap_region)
        median_overlap_region = np.median(stacked_overlap_region, axis=0)
        new_image[:image_size*overlap_rate,image_size*(1-overlap_rate):] = np.ceil(median_overlap_region).astype(int)
    #top left
    if i >= col_num and i%col_num != 0:
        image_tl = images[i-col_num-1].copy()
        overlap_region = [image[:image_size*overlap_rate,:image_size*overlap_rate],image_left[:image_size*overlap_rate,image_size*(1-overlap_rate):],image_top[image_size*(1-overlap_rate):,:image_size*overlap_rate],
                         image_tl[image_size*(1-overlap_rate):,image_size*(1-overlap_rate):]]
        stacked_overlap_region = np.stack(overlap_region)
        median_overlap_region = np.median(stacked_overlap_region, axis=0)
        new_image[:image_size*overlap_rate,:image_size*overlap_rate] = np.ceil(median_overlap_region).astype(int)
    result_image_file_path = 'overlap_modified_median_subtracted/' + image_name + '.npy'
    np.save(result_image_file_path, new_image)

#stitch images
image_dir = 'overlap_modified_median_subtracted/'
image_names = os.listdir(image_dir)
image_names.sort()
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

image_names_sub = [s for s in image_names if time_point in s]
reshaped_images = []
for name in image_names_sub:
    tiff_file_path = image_dir + name
    image = np.load(tiff_file_path)
    reshaped_images.append(image)
images = np.stack(reshaped_images)

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
result_image_file_path = 'stitched_images_BFP/' + 'stitched_BFP.npy'
np.save(result_image_file_path, stitched_image)

    
