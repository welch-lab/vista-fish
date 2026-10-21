This folder contains the code for iNeuron+Lysosome data analysis

BFP_image_preprocess/image_stitching_BFP.py: stitch BFP channel images

TrackMate_code/trackmate_jython.sh & TrackMate_code/trackmate_jython_all.py: track lysosome movement using TrackMate

TrackMate_code/match_track_with_cell.ipynb: match TrackMate tracks with cell ID

brightfield-segmentation: run CellPose SAM and CellProfiler to segment neurons with neurite structures, and calculate morphology features of the new neuron segmentations

Perturbed_iNeuron_Preprocessing.py: cell annotation of the iNeuron sample

Perturbed_iNeuron_DEGs.ipynb: find the DEGs of the iNeuron sample

Figure6_BFP.ipynb: compare the mean BFP intensity of gene-perturbed cells with unperturbed cells

Figure6_lysosome_track_metrics.ipynb: compare the lysosome track metrics of each gene knockdown group with the negative control group 

Figure6_lysosome_kymograph.ipynb: plot kymograph of lysosome movement for example cells in gene knockdown groups with significant gene perturbation effect, compared to the negative control group

assign_transcripts_to_cellpose_masks.ipynb: assign transcripts to the neuron segmentation result we get from CellPose SAM, generate the RNA count matrix for the new neuron segmentation dataset

Figure6_morphology_metrics.ipynb: compare the morphology metrics of each gene knockdown group with the negative control group

Figure6_draw_morphology.ipynb: generate images of soma and neurite boundary outlines overlaid on brightfield images

Figure6_BFP.ipynb, Figure6_lysosome_track_metrics.ipynb, Figure6_morphology_metrics.ipynb, and Perturbed_iNeuron_DEGs.ipynb call a cell perturbed for a gene if any of that gene's sgRNA probes has a raw count > 1, excluding the two most abundant probes genome-wide (likely off-target/background).

Figure6_lysosome_track_metrics.ipynb and Figure6_morphology_metrics.ipynb also report each metric's effect size as (perturbed mean - control mean) / pooled SD, alongside percent change.
