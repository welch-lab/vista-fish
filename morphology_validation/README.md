# Morphology Validation: Local Cell Density Confound Check

Validation analysis testing whether soma and neurite morphology differences between CRISPR-perturbed and negative-control iNeurons (CellProfiler features on Cellpose-segmented brightfield images) are confounded by local cell density, for two gene perturbations.

For each cell, local density is the number of other QC-passed cells within a 150 px radius of its centroid in the same field of view, independent of BFP status. Cells are progressively filtered by excluding the most locally dense cells (top 5/10/15/20/30/40/50/60/70%) and the perturbation-vs-control comparison (t-test and Mann-Whitney U, BH-corrected) is repeated at each threshold, to check that the effect sizes and significance calls are not an artifact of crowding.

## Contents

- `density_filtering_onecut1_vs_control.ipynb` — OneCut1 vs. negative control.
- `density_filtering_isl1_vs_control.ipynb` — ISL1 vs. negative control.

Both notebooks share the same structure: load CellProfiler soma/neurite features and BFP intensity, compute local neighbor count, call BFP-positive cells with a per-condition Gaussian mixture model, then sweep the densest-cell exclusion threshold and report effect size and significance (Area, Solidity, NeuriteLength, plus 3 additional neurite-contact metrics) at each one.
