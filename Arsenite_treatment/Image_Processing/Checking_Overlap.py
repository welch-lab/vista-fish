# Checking Overlapping % CQ1 Yokogawa Images. Skip if overlapping % is known

from os import path
import numpy as np
import pandas as pd
import m2stitch
import os
import cv2
import tifffile
import imageio
import matplotlib.pyplot as plt
from skimage.transform import resize
from skimage import exposure
from skimage.feature import match_template
import glob
import re

BASE = ("/path/to/Image")
SUFFIX = "T0001Z001C3.tif"

# Checking # of field of views
field_dirs = sorted(glob.glob(os.path.join(BASE, "Field*")),
                    key=lambda p: int(re.search(r"Field(\d+)", p).group(1)))
n_fields = len(field_dirs)
print(f"{n_fields} fields")
print("possible grids (rows x cols):",
      [(n_fields // a, a) for a in range(1, n_fields + 1) if n_fields % a == 0
       and 4 <= a <= n_fields // 2])

def load_field(idx):
    """idx is 0-based position in the sorted field list."""
    hits = glob.glob(os.path.join(field_dirs[idx], f"*{SUFFIX}"))
    if len(hits) != 1:
        raise RuntimeError(f"{field_dirs[idx]}: {len(hits)} matches")
    return tifffile.imread(hits[0]).squeeze().astype(np.float32)

A = load_field(0)
H, W = A.shape
print(f"tile {A.shape} {A.dtype}, range {A.min():.0f}-{A.max():.0f}")

def find_overlap(a, b, axis, template_px=24, max_frac=0.15, trim=150):
    if axis == 0:
        a, b = a.T, b.T           
    Hh, Ww = a.shape
    T = template_px
    S = int(Ww * max_frac)

    template = a[trim:Hh - trim, Ww - T:]
    search = b[:, :S]
    ncc = match_template(search, template)
    py, px = np.unravel_index(np.argmax(ncc), ncc.shape)
    return px + T, py - trim, float(ncc.max())

B = load_field(1)      
ov_px, perp_shift, score = find_overlap(A, B, axis=1)

print(f"horizontal overlap : {ov_px} px  ({100 * ov_px / W:.2f}% of {W})")
print(f"vertical drift     : {perp_shift:+d} px")
print(f"NCC peak           : {score:.3f}")
print()
print("NCC above ~0.5 means a real match. Below ~0.3, these two tiles are")
print("probably not neighbours -- the field numbering may not be raster,")
print("or field 2 sits in a different row.")

print("\ntemplate width -> overlap estimate")
for t in (12, 16, 24, 32, 48, 64):
    o, d, s = find_overlap(A, B, axis=1, template_px=t)
    print(f"  {t:3d} px -> {o:4d} px ({100*o/W:5.2f}%)  dy={d:+4d}  ncc={s:.3f}")

# checking multiple possible overlapping %
candidates = [a for a in range(4, n_fields // 2 + 1) if n_fields % a == 0]
print("stride  overlap%   dx    ncc")
best = None
for stride in candidates:
    if stride >= n_fields:
        continue
    try:
        D = load_field(stride)
    except Exception as e:
        print(f"{stride:6d}  load failed: {e}")
        continue
    o, d, s = find_overlap(A, D, axis=0)
    flag = ""
    if best is None or s > best[1]:
        best, flag = (stride, s), "  <-- best so far"
    print(f"{stride:6d}  {100*o/H:6.2f}%  {d:+4d}  {s:.3f}{flag}")

print(f"\nCOLS is most likely {best[0]}  (rows = {n_fields // best[0]})")
print("Sanity check: the vertical overlap % should be close to the")
print("horizontal one from cell 2. If it isn't, the stride is wrong.")

# Checking the overlap with an image
COLS = 18                 # number of columns
OV_H = ov_px                   # measured, not assumed
OV_V, _, _ = find_overlap(A, load_field(COLS), axis=0)

print(f"using COLS={COLS}, horizontal overlap {OV_H} px, "
      f"vertical {OV_V} px")

tl = A
tr = load_field(1)
bl = load_field(COLS)
br = load_field(COLS + 1)

step_x, step_y = W - OV_H, H - OV_V
mh, mw = step_y + H, step_x + W
mosaic = np.zeros((mh, mw), np.float32)

for img, (y, x) in [(tl, (0, 0)), (tr, (0, step_x)),
                    (bl, (step_y, 0)), (br, (step_y, step_x))]:
    mosaic[y:y+H, x:x+W] = img

lo, hi = np.percentile(mosaic[mosaic > 0], (1, 99.5))
plt.figure(figsize=(13, 13))
plt.imshow(mosaic, cmap="gray", vmin=lo, vmax=hi)
plt.axvline(step_x, color="r", lw=0.6)
plt.axvline(W, color="r", lw=0.6)
plt.axhline(step_y, color="r", lw=0.6)
plt.axhline(H, color="r", lw=0.6)
plt.title(f"2x2, overlap {100*OV_H/W:.2f}% h / {100*OV_V/H:.2f}% v "
          f"(red = overlap band edges)")
plt.axis("off")
plt.show()

# Zoomed in 2 x 2 junction
z = 250
seam = mosaic[step_y - z:step_y + z, step_x - z:step_x + z]
lo, hi = np.percentile(seam, (1, 99.5))
plt.figure(figsize=(9, 9))
plt.imshow(seam, cmap="gray", vmin=lo, vmax=hi, interpolation="nearest")
plt.axvline(z, color="r", lw=0.5, alpha=0.6)
plt.axhline(z, color="r", lw=0.5, alpha=0.6)
plt.title("centre of the 2x2 junction")
plt.axis("off")
plt.show()
