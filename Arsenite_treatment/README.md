# Arsenite Treatment + CRISPRi screening + live-imaging (lysosome movement) + Stress Granule Immunostaining + Xenium iNeuron Analysis

Analysis code for a pooled CRISPRi screen in arsenite-treated iNeurons, with sgRNA detection by Xenium in situ imaging and stress granule immunostaining.

##
Brief description of the experiment: iPSC-derived neurons (iNeurons), CRISPRi gRNAs targeting 24 transcription factors, arsenite treatment, and lysosome movement

## Contents

- `Fig7_4-6.py` — sgRNA outlier detection, pool cross-specificity (sensitivity/specificity), and raw expected-vs-unexpected guide counts.
- `Image_Processing/` — brightfield/IHC image stitching and FOV-overlap handling.
- `run_4well_pipeline_lysosome.py` — assigns TrackMate lysosome tracks to Xenium-segmented cells per well (own affine matrix per well, FOV-by-FOV point-in-polygon against each cell's soma mask), producing the `combined_<well>.csv` / `combined_all_wells.csv` tables the notebooks below read.
- `qc_well_grna_analysis.ipynb` — per-cell/per-well QC, affine-based well assignment, gRNA probe identification and per-cell gRNA calling, per-gRNA/per-gene cell counts, and a knockdown-validation check (gRNA vs. its own target transcript, for the 5 genes with a plain expression probe in the panel).
- `bfp_intensity_gRNA_analysis.ipynb` — per-cell BFP intensity extraction from the stitched, affine-registered live-imaging channel; off-target guide screening; and the BFP+/BFP- intensity cutoff found via ROC/Youden's J against the `sgRNA_total_count > 1` perturbation label.
- `joint_regression_nc_baseline_with_kymographs.ipynb` — TrackMate track-count QC, a joint treatment x gene-perturbation regression (one model per lysosome motility metric, population restricted to a `Negative_control` baseline), effect-size bar plots and a combined significance heatmap, and BNC2 kymographs built from tracks sampled within each group's own metric percentile range (no minimum track-length cutoff).
