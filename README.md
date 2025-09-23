# VISTA-FISH
Contains all scripts for datasets analysis from the VISTA-FISH protocol, as stated in the paper.

Scripts are grouped by datasets and functions. SMAD+GCAMP contains all preprocessing, clustering, Calcium imaging, RNA velocity, and other analysis that was conducted on the SMAD neuron dataset. Brightfield-segmentation folder contains all scripts and instructions pertaining to the segmentations conducted with Cellpose and Cellprofiler on the perturbed iNeurons. The iNeuron+Lysosome folder contains code that generated all of the data and figures in Figure 6 of the paper such as the lysosome tracking and morphologies, along with the preprocessing and clustering of the perturbed iNeuron dataset. 

For each folder, users should begin with preprocessing, dimentionality reduction, and clustering of the dataset with the provided script and then move onto the analysis for each folder:
1) SMAD+GCAMP:GCaMP_image_preprocess contains files on the calcium imaging/gcamp signal processing script as well as the stitching of the images. Data_preprocess contains scripts for preprocessing of the SMAD dataset as well as the IF extractions and the nucleus and cytoplasm RNA                 extractions. Contains scripts for generating Figures 3 and 4
2) brightfield-segmentation: Contains all scripts, github page, demos, pipelines for the IF cell segmentations conducted with Cellpose and Cellprofiler
3) iNeuron+Lysosome: BFP_image_preprocess contains the script on the xenium IF image stitching. TrackMate_code contains scripts on utilizing trackmate in detecting and analyzing lysosomal movements. Contains scripts on preprocessing of the perturbed iNeuron dataset as well as generating the plots and graphs in figure 6 of the paper


    
