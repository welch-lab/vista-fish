# VISTA-FISH
Contains all scripts for datasets analysis from the VISTA-FISH protocol, as stated in the paper.

Scripts are grouped by datasets and functions. SMAD+GCAMP contains all preprocessing, clustering, Calcium imaging, RNA velocity, and other analysis that was conducted on the SMAD neuron dataset. Brightfield-segmentation folder contains all scripts and instructions pertaining to the segmentations conducted with Cellpose and Cellprofiler on the perturbed iNeurons. The iNeuron+Lysosome folder contains code that generated all of the data and figures in Figure 6 of the paper such as the lysosome tracking and morphologies, along with the preprocessing and clustering of the perturbed iNeuron dataset. 

For each folder, users should begin with preprocessing, dimentionality reduction, and clustering of the dataset with the provided script and then move onto the analysis for each folder:
1) SMAD+GCAMP:
    
