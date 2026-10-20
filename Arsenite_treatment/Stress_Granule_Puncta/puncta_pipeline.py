"""Shared constants and functions for registering the 4-well SG-puncta plate
to Xenium and calling puncta per cell. Used by stress_granule_puncta_analysis.ipynb.

Coordinate systems:
  Xenium um    -- transcripts.x_location/y_location, boundaries.vertex_x/y
  Xenium px    -- Xenium um / PIXEL_SIZE_XEN (0.2125 um/px)
  SG mosaic px -- the stitched-plate pixel grid (39 rows x 18 cols of
                  2000x2000 FOVs, 10.75% overlap). Xenium px -> SG mosaic px
                  via the INVERSE of a well's affine matrix (matrix.csv maps
                  SG-mosaic-px -> Xenium-px, fit per well from manually
                  matched landmark pairs, keypoints.csv).
  FOV-local px -- pixel coords within one 2000x2000 FOV. SG mosaic px ->
                  FOV-local px by subtracting that FOV's origin in the mosaic.

Well membership is decided from each well's own keypoints.csv bounding box in
Xenium px -- the four boxes are cleanly separated.
"""
import glob
import os
import zipfile

import numpy as np
import pandas as pd
import tifffile
from skimage.draw import polygon as sk_polygon
from skimage.filters import gaussian, threshold_otsu
from skimage.measure import label, regionprops
from skimage.morphology import disk, white_tophat

# --------------------------------------------------------------------------
# paths -- replace with your own
# --------------------------------------------------------------------------
XENIUM_DIR = "path/to/xenium_output"
IMAGE_ROOT_MEDIAN_SUB = "path/to/median_subtracted_images"
AFFINE_ZIP_DIR = "path/to/affine_matrices"
AFFINE_EXTRACT_DIR = "path/to/affine_matrices_extracted"

PIXEL_SIZE_XEN = 0.2125  # um/px, Xenium
PUNCTA_TOKEN = "C4"      # stress-granule channel

FOV_SIZE = 2000
N_ROWS, N_COLS = 39, 18
OVERLAP = 0.1075
STEP = int(round(FOV_SIZE * (1 - OVERLAP)))  # 1785

# rows/cols are 1-indexed, inclusive (row 1 = FOVs 1-18)
WELL_GRID = {
    "top_left":     dict(rows=(2, 17),  cols=(2, 7),   treatment="arsenite"),
    "top_right":    dict(rows=(2, 17),  cols=(12, 17), treatment="arsenite"),
    "bottom_left":  dict(rows=(23, 38), cols=(2, 7),   treatment="untreated"),
    "bottom_right": dict(rows=(23, 38), cols=(12, 17), treatment="untreated"),
}
WELL_NAMES = list(WELL_GRID.keys())

TOPHAT_RADIUS = 6
SMOOTH_SIGMA = 1.0
MIN_SNR = 5.0
MIN_AREA = 5
MAX_AREA = 600


def fov_num(row, col):
    return (row - 1) * N_COLS + col


def fov_rowcol(fov):
    row, col0 = divmod(fov - 1, N_COLS)
    return row + 1, col0 + 1


def fov_list_for_well(well):
    g = WELL_GRID[well]
    r0, r1 = g["rows"]
    c0, c1 = g["cols"]
    return [fov_num(r, c) for r in range(r0, r1 + 1) for c in range(c0, c1 + 1)]


def fov_origin_in_mosaic(fov):
    """(y0, x0) of this FOV's top-left corner in the full stitched mosaic."""
    row, col = fov_rowcol(fov)
    return (row - 1) * STEP, (col - 1) * STEP


# --------------------------------------------------------------------------
# affines
# --------------------------------------------------------------------------
def _ensure_affines_extracted():
    os.makedirs(AFFINE_EXTRACT_DIR, exist_ok=True)
    for well in WELL_NAMES:
        d = os.path.join(AFFINE_EXTRACT_DIR, well)
        if os.path.exists(os.path.join(d, "matrix.csv")):
            continue
        zpath = os.path.join(AFFINE_ZIP_DIR, f"{well}_SGs.zip")
        with zipfile.ZipFile(zpath) as zf:
            zf.extractall(d)


def load_well_affine(well):
    """Return dict with matrix A (SG-mosaic-px -> Xenium-px), its inverse,
    the scale factor, derived SG pixel size, and the keypoint bbox in Xenium
    px (that well's footprint for cell assignment)."""
    _ensure_affines_extracted()
    d = os.path.join(AFFINE_EXTRACT_DIR, well)
    A = np.loadtxt(os.path.join(d, "matrix.csv"), delimiter=",")
    A_inv = np.linalg.inv(A)
    kp = pd.read_csv(os.path.join(d, "keypoints.csv"))
    scale = float(np.sqrt(abs(np.linalg.det(A[:2, :2]))))
    bbox = (
        float(kp["fixedX"].min()), float(kp["fixedX"].max()),
        float(kp["fixedY"].min()), float(kp["fixedY"].max()),
    )
    return {
        "A": A, "A_inv": A_inv, "keypoints": kp, "scale": scale,
        "sg_pixel_size_um": scale * PIXEL_SIZE_XEN,
        "bbox_xen_px": bbox,
    }


def xen_px_to_mosaic_px(xy, A_inv):
    """(N,2) Xenium px [x,y] -> (N,2) SG stitched-mosaic px [x,y]."""
    homo = np.hstack([xy, np.ones((len(xy), 1))])
    return homo.dot(A_inv.T)[:, :2]


# --------------------------------------------------------------------------
# cell / nucleus boundaries + well/FOV assignment
# --------------------------------------------------------------------------
def load_boundaries_px(kind):
    """kind: 'cell' or 'nucleus'. Returns df with x_px, y_px columns added."""
    fname = "cell_boundaries.csv.gz" if kind == "cell" else "nucleus_boundaries.csv.gz"
    df = pd.read_csv(os.path.join(XENIUM_DIR, fname))
    df["x_px"] = df["vertex_x"] / PIXEL_SIZE_XEN
    df["y_px"] = df["vertex_y"] / PIXEL_SIZE_XEN
    return df


def cell_centroids_px(cell_boundaries_px):
    g = cell_boundaries_px.groupby("cell_id")[["x_px", "y_px"]].mean()
    g.columns = ["cx_px", "cy_px"]
    return g


def assign_cells_to_wells(centroids, affines, pad_px=0.0):
    """centroids: index=cell_id, cols cx_px/cy_px. affines: well -> load_well_affine() result.
    Returns well -> set(cell_id)."""
    out = {}
    for well, aff in affines.items():
        x0, x1, y0, y1 = aff["bbox_xen_px"]
        m = ((centroids["cx_px"] >= x0 - pad_px) & (centroids["cx_px"] <= x1 + pad_px) &
             (centroids["cy_px"] >= y0 - pad_px) & (centroids["cy_px"] <= y1 + pad_px))
        out[well] = set(centroids.index[m])
    return out


def assign_cells_to_fovs(well, cell_ids, cell_boundaries_px, A_inv):
    """A cell is assigned to the single FOV whose 2000x2000 window fully
    contains its transformed boundary bbox. Cells straddling the border
    between two overlapping FOVs (not fully inside either) are dropped.

    Returns (assigned: dict cell_id -> fov, n_candidates: int).
    """
    sub = cell_boundaries_px[cell_boundaries_px["cell_id"].isin(cell_ids)]
    xy = sub[["x_px", "y_px"]].to_numpy()
    mxy = xen_px_to_mosaic_px(xy, A_inv)
    sub = sub.assign(mx=mxy[:, 0], my=mxy[:, 1])
    agg = sub.groupby("cell_id").agg(
        mx_min=("mx", "min"), mx_max=("mx", "max"),
        my_min=("my", "min"), my_max=("my", "max"),
    )

    assigned = {}
    for fov in fov_list_for_well(well):
        y0, x0 = fov_origin_in_mosaic(fov)
        m = ((agg["mx_min"] >= x0) & (agg["mx_max"] < x0 + FOV_SIZE) &
             (agg["my_min"] >= y0) & (agg["my_max"] < y0 + FOV_SIZE))
        for cid in agg.index[m]:
            assigned.setdefault(cid, fov)
    return assigned, len(agg)


# --------------------------------------------------------------------------
# puncta detection
# --------------------------------------------------------------------------
def load_fov_image_median_subtracted(fov):
    """Median-flat-field-subtracted C4 image for one FOV (background = the
    per-pixel median across all FOVs at this channel, subtracted per FOV,
    clipped >= 0). Still needs the tophat pass below before thresholding."""
    path = os.path.join(IMAGE_ROOT_MEDIAN_SUB, f"Field{fov:04d}", "C4_median_subtracted.npy")
    if not os.path.exists(path):
        return None
    return np.load(path)


def subtract_background(img, radius=TOPHAT_RADIUS, smooth_sigma=SMOOTH_SIGMA):
    x = gaussian(img, smooth_sigma, preserve_range=True) if smooth_sigma else img
    return white_tophat(x, footprint=disk(radius))


def bbox_from_polygons(*xy_arrays, pad=15, fov_size=FOV_SIZE):
    """Padded, clipped integer bbox (y0, y1, x0, x1) covering all given
    (x, y) arrays (e.g. cell + nucleus vertices)."""
    all_x = np.concatenate([a[:, 0] for a in xy_arrays if len(a)])
    all_y = np.concatenate([a[:, 1] for a in xy_arrays if len(a)])
    x0 = max(0, int(np.floor(all_x.min())) - pad)
    x1 = min(fov_size, int(np.ceil(all_x.max())) + pad)
    y0 = max(0, int(np.floor(all_y.min())) - pad)
    y1 = min(fov_size, int(np.ceil(all_y.max())) + pad)
    return y0, y1, x0, x1


def polygon_mask_crop(xs, ys, y0, x0, crop_shape):
    rr, cc = sk_polygon(ys - y0, xs - x0, shape=crop_shape)
    m = np.zeros(crop_shape, dtype=bool)
    m[rr, cc] = True
    return m


def detect_in_cell_crop(flat_crop, cell_mask, min_snr=MIN_SNR, min_area=MIN_AREA, max_area=MAX_AREA):
    """Per-cell puncta call: threshold = max(Otsu, median + min_snr*1.4826*MAD)
    over the cell's own pixels, then connected components kept if their area
    falls in [min_area, max_area] px. Returns (punctum_label_crop, threshold, frac_above)."""
    n_cell = int(cell_mask.sum())
    if n_cell == 0:
        return np.zeros(flat_crop.shape, np.int32), np.nan, 0.0

    vals = flat_crop[cell_mask]
    try:
        thr = float(threshold_otsu(vals))
    except (ValueError, RuntimeError):
        thr = float(vals.max()) + 1.0
    mad = float(np.median(np.abs(vals - np.median(vals))))
    floor = float(np.median(vals)) + min_snr * 1.4826 * mad
    thr = max(thr, floor)

    above = np.zeros(flat_crop.shape, dtype=bool)
    above[cell_mask] = flat_crop[cell_mask] > thr
    frac_above = above.sum() / n_cell

    lab = label(above)
    keep = np.zeros(flat_crop.shape, dtype=np.int32)
    for p in regionprops(lab):
        if min_area <= p.area <= max_area:
            keep[lab == p.label] = p.label
    return keep, thr, float(frac_above)
