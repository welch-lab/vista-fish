This folder contains the code for SMAD+GCAMP data analysis

GCaMP_image_preprocess/image_stitching_gcamp.py: stitch GCaMP images for all 40 time frames
GCaMP_image_preprocess/gcamp_signal.py: calculate the mean GCaMP intensity of cells using the stitched GCaMP image and the affine matrix that we get from image registration using Xenium Explorer

data_preprocess/SMAD_preprocessing.py: cell annotation of the SMAD sample
data_preprocess/image_stitching_IF.py: stitch IF staining images for both dendrite and axon stainings
data_preprocess/IF_staining_neurite_extraction.py: extract pixels that have higher intensities than the threshold from dendrite and axon staining images, match the pixels with transcripts, run DEG analysis for the RNA expression in neurites compared to RNA expression in cells
data_preprocess/filament_extraction.py: extract pixels that have higher intensities than the threshold from the filament layer of the Xenium image, match the pixels with transcripts, run DEG analysis for the RNA expression in the filament compared to RNA expression in cells
data_preprocess/extract_nucleus_cytoplasm_RNA.py: for each cell, calculate the RNA count in the nucleus and cytoplasm; the counts will be used in RNA velocity analysis

Figure2_velocity.ipynb: run RNA velocity analysis using nucleus and cytoplasm RNA counts, generate the plots used in Figure 2 and S2
Figure2_4_suppFig_boxplots_umaps.ipynb: generate UMAPs and boxplots in Figure 2, 4, and S3
Figure3_Cluster_by_gcamp.ipynb: calculate GCaMP signal features for cells, cluster cells using k-means based on GCaMP signal features, and generate plots in Figure 3
Figure4_LASSO.ipynb: prepare the genes used in LASSO regression, and check the performance of LASSO regression prediction
Figure4_LASSO.html: run LASSO regression using R
Figure5_neurite_DEG.html: run GO analysis for the DEGs we found in neurites, and compare the UTR length of neurite DEGs and cell DEGs
Figure5_filament_DEG.html: run GO analysis for the DEGs we found in filament, and compare the UTR length of filament DEGs and cell DEGs
FigureS3_neighboring_cells.ipynb: check the composition of cell types for neighboring cells, generate the plot used in Figure S3
