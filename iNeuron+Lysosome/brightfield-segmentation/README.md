# Brightfield Neurite & Soma Segmentation  
Example Cellpose‑SAM + CellProfiler workflow
==================================================================

This repository contains:

* **Custom Cellpose‑SAM models** for neurite & soma segmentation of bright‑field iPSC‑derived neurons.  
* Example **Jupyter notebooks** that run inference on tiled images and save corresponding masks.  
* A template **CellProfiler v4.2.8** pipeline (`ipsc_neuron_tiles.cpproj`).  
* A **minimal example dataset** without downloading gigabytes of raw data.

---

## Quick installation

### Option 1 – recommended: create the Conda env from `environment.yml`

```bash
conda env create -f environment.yml
conda activate bf_seg
```

### Option 2 – manual install

```bash
# 1) create & activate a fresh env
conda create -n bf_seg python=3.10
conda activate bf_seg

# 2) upgrade pip (good hygiene)
python -m pip install --upgrade pip

# 3) install the core packages
python -m pip install     "cellpose[gui]" notebook matplotlib tqdm imageio fastremap fill-voids natsort
```

### (Optional, Recommended) GPU acceleration  
If you have an NVIDIA card, we recommend swapping the default CPU‑only Torch for a CUDA build  
(choose the CUDA tag that matches your driver; `cu126` is just an example):

```bash
pip uninstall torch
pip install torch torchvision --index-url https://download.pytorch.org/whl/cu126
```

---

## End‑to‑end example

```text
.
├── data/
│   └── example/
│        ├── raw/               ← input bright‑field tiles (mini subset)
│        └── interim/
│             ├── neurite/      ← where neurite masks will be saved
│             └── soma/         ← where soma masks will be saved
├── notebooks/
│   ├── cellpose_tile_neurite.ipynb
│   └── cellpose_tile_soma.ipynb
└── pipelines/
    └── ipsc_neuron_tiles.cpproj
```

1. **Run the notebooks**

   ```bash
   jupyter notebook notebooks/cellpose_tile_neurite.ipynb
   jupyter notebook notebooks/cellpose_tile_soma.ipynb
   ```

2. **Open the pipeline in CellProfiler 4.2.8**

   *File ► Open Project* → `pipelines/ipsc_neuron_tiles.cpproj`

   When prompted:

   | Dialog prompt   | What to pick                             |
   |-----------------|------------------------------------------|
   | **Default input folder**  | `data/example/raw`              |
   | **Default output folder** | `data/example/processed` (or any writable path) |

3. Click **Analyze Images**. All outputs export under your chosen output folder.

For further details, please refer to the workflow diagram  in [`reports/workflow_diagram.pdf`](reports/workflow_diagram.pdf).

---

## Citation

If you use this workflow in a paper, please cite **Lee et al., 2025 (in prep.)**  
and the original Cellpose & CellProfiler publications.

---

© 2025 The Welch Lab & contributors – released under the GNU GENERAL PUBLIC LICENSE.