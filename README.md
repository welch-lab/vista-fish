# VISTA-FISH: Video Imaging with Spatial-Temporal Analysis by Fluorescent In Situ Hybridization
Contains all scripts for datasets analysis from the VISTA-FISH protocol, as stated in the paper.

Scripts are grouped by datasets and functions. SMAD+GCAMP contains all preprocessing, clustering, Calcium imaging, RNA velocity, and other analyses conducted on the SMAD neuron dataset. The iNeuron+Lysosome folder contains code that generated all of the data and figures in Figure 6 of the paper, such as the lysosome tracking and morphologies, along with the preprocessing and clustering of the perturbed iNeuron dataset. 

For each folder, users should begin with preprocessing, dimensionality reduction, and clustering of the dataset with the provided script and then move on to the analysis for each folder:
1) SMAD+GCAMP: GCaMP_image_preprocess contains files on the calcium imaging/gcamp signal processing script as well as the stitching of the images. Data_preprocess contains scripts for preprocessing of the SMAD dataset as well as the Immunofluorescence (IF) extractions and the nucleus and cytoplasm RNA extractions. Contains scripts for generating Figures 2-5.
2) iNeuron+Lysosome: BFP_image_preprocess contains the script for the Xenium BFP channel image stitching. TrackMate_code contains scripts on utilizing TrackMate in detecting and analyzing lysosomal movements. Contains scripts for preprocessing the perturbed iNeuron dataset as well as generating the plots and graphs in Figure 6 of the paper. Brightfield-segmentation folder contains all scripts and instructions pertaining to the segmentations conducted with Cellpose and Cellprofiler on the perturbed iNeurons.
3) Arsenite-treated iNeuron+Lysosome+Stress Granules: 


    
