"""
Assign TrackMate lysosome tracks to Xenium-segmented cells, separately for
each of the 4 wells on the 20260615 4-well arsenite CQ1 plate, then combine.

Grid geometry:
  - Plate is tiled into 38 rows x 18 cols = 684 FOVs, numbered row-major
    starting at 1 (row 1 = FOVs 1-18).
  - FOV size 2000x2000 px, overlap 11% -> step between FOV origins = 1780 px.
  - FOV_num -> row = (FOV_num-1)//18 + 1, col = (FOV_num-1)%18 + 1
  - row,col  -> FOV_num = (row-1)*18 + col

Wells (1-indexed row/col, inclusive):
  top_left:     rows 2-16,  cols 2-7    (arsenite)
  top_right:    rows 2-16,  cols 12-17  (arsenite)
  bottom_left:  rows 23-37, cols 2-7    (untreated)
  bottom_right: rows 23-37, cols 12-17  (untreated)

Each well was registered independently to the Xenium image (own affine
matrix), producing 4 separate matrix.csv files under this project dir.
"""
import os
import time
import numpy as np
import pandas as pd
from skimage.draw import polygon

PROJECT_DIR = "path/to/project"
XENIUM_DIR = "path/to/xenium_output"
TRACKMATE_DIR = "path/to/trackmate_results"
OUT_DIR = "path/to/trackcell_output"

PIXEL_SIZE_XENIUM = 0.2125  # um/px
FOV_SIZE = 2000
N_COLS = 18
OVERLAP = 0.11
STEP = int(round(FOV_SIZE * (1 - OVERLAP)))  # 1780

WELLS = {
    "top_left": dict(rows=(2, 16), cols=(2, 7),
                      affine_dir="lyso_top_left", treatment="arsenite"),
    "top_right": dict(rows=(2, 16), cols=(12, 17),
                       affine_dir="top_right_alignment_files", treatment="arsenite"),
    "bottom_left": dict(rows=(23, 37), cols=(2, 7),
                         affine_dir="bottom_left", treatment="untreated"),
    "bottom_right": dict(rows=(23, 37), cols=(12, 17),
                          affine_dir="bottom_right_alignment_files", treatment="untreated"),
}


def fov_num(row, col):
    return (row - 1) * N_COLS + col


def fov_list_for_well(well):
    r0, r1 = well["rows"]
    c0, c1 = well["cols"]
    return [fov_num(r, c) for r in range(r0, r1 + 1) for c in range(c0, c1 + 1)]


def load_cell_boundaries():
    path = os.path.join(XENIUM_DIR, "cell_boundaries.csv.gz")
    df = pd.read_csv(path)
    return df


def transform_cell_boundaries(cell_boundaries_csv, affine_dir):
    matrix = np.loadtxt(os.path.join(PROJECT_DIR, affine_dir, "matrix.csv"), delimiter=",")
    matrix_inv = np.linalg.inv(matrix)

    df = cell_boundaries_csv.copy()
    df["vertex_x_coord"] = df["vertex_x"] / PIXEL_SIZE_XENIUM
    df["vertex_y_coord"] = df["vertex_y"] / PIXEL_SIZE_XENIUM

    homogeneous = np.hstack([
        df[["vertex_x_coord", "vertex_y_coord"]].values,
        np.ones((df.shape[0], 1)),
    ])
    transformed = homogeneous.dot(matrix_inv.T)
    df["x_transformed"] = np.round(transformed[:, 0]).astype(int)
    df["y_transformed"] = np.round(transformed[:, 1]).astype(int)
    return df


def process_fov(fov_number, cell_boundaries_csv, cell_boundaries_agg, well_name, treatment):
    row = (fov_number - 1) // N_COLS + 1
    col = (fov_number - 1) % N_COLS + 1
    x_sub = (col - 1) * STEP
    y_sub = (row - 1) * STEP

    field_str = f"{fov_number:04d}"
    base = f"Field{field_str}_stack"
    spots_path = os.path.join(TRACKMATE_DIR, f"Field{field_str}", f"{base}-spots.csv")
    tracks_path = os.path.join(TRACKMATE_DIR, f"Field{field_str}", f"{base}-tracks.csv")

    valid_ids = cell_boundaries_agg[
        (cell_boundaries_agg["x_transformed_min"] >= x_sub) &
        (cell_boundaries_agg["x_transformed_max"] < x_sub + FOV_SIZE) &
        (cell_boundaries_agg["y_transformed_min"] >= y_sub) &
        (cell_boundaries_agg["y_transformed_max"] < y_sub + FOV_SIZE)
    ].index

    if len(valid_ids) == 0:
        return None

    df_fov = cell_boundaries_csv[cell_boundaries_csv["cell_id"].isin(valid_ids)].copy()
    df_fov["x_transformed_subed"] = df_fov["x_transformed"] - x_sub
    df_fov["y_transformed_subed"] = df_fov["y_transformed"] - y_sub

    spots = pd.read_csv(spots_path, header=0, skiprows=[1, 2, 3], index_col=0)
    tracks = pd.read_csv(tracks_path, header=0, skiprows=[1, 2, 3], index_col=0)

    if spots.empty or tracks.empty:
        return None

    spots["POSITION_X_p"] = (spots["POSITION_X"] * 100).clip(upper=FOV_SIZE - 1)
    spots["POSITION_Y_p"] = (spots["POSITION_Y"] * 100).clip(upper=FOV_SIZE - 1)
    spots = spots.assign(
        col=spots["POSITION_X_p"].round().astype(int),
        row=spots["POSITION_Y_p"].round().astype(int),
    )

    height, width = (FOV_SIZE, FOV_SIZE)
    df_list = []
    for cell_id, grp in df_fov.groupby("cell_id"):
        xs = grp["x_transformed_subed"].to_numpy()
        ys = grp["y_transformed_subed"].to_numpy()
        cols_arr = np.round(xs).astype(int)
        rows_arr = np.round(ys).astype(int)
        rr, cc = polygon(rows_arr, cols_arr, shape=(height, width))
        df_list.append(pd.DataFrame({"cell_id": cell_id, "row": rr, "col": cc}))

    if not df_list:
        return None
    pixels_df = pd.concat(df_list, ignore_index=True)

    merged = spots.merge(pixels_df, on=["row", "col"], how="left")
    total_pos = spots.groupby("TRACK_ID").size().rename("total_positions")
    in_cell = (
        merged.dropna(subset=["cell_id"])
        .groupby(["TRACK_ID", "cell_id"])
        .size()
        .rename("positions_in_cell")
    )
    if in_cell.empty:
        return None

    result = in_cell.reset_index().merge(total_pos.reset_index(), on="TRACK_ID")
    result = result[
        (result["positions_in_cell"] > 0) &
        (result["positions_in_cell"] < result["total_positions"])
    ].reset_index(drop=True)
    if result.empty:
        return None

    merged_result = result.merge(tracks, how="left", on="TRACK_ID")
    merged_result = merged_result[merged_result["TRACK_DISPLACEMENT"] > 0.05]
    if merged_result.empty:
        return None

    merged_result["FOV"] = fov_number
    merged_result["well"] = well_name
    merged_result["treatment"] = treatment
    return merged_result


def run_well(well_name, well, cell_boundaries_csv):
    t0 = time.time()
    print(f"[{well_name}] transforming cell boundaries with {well['affine_dir']}/matrix.csv", flush=True)
    df = transform_cell_boundaries(cell_boundaries_csv, well["affine_dir"])
    agg = df.groupby("cell_id").agg({
        "x_transformed": ["min", "max"],
        "y_transformed": ["min", "max"],
    })
    agg.columns = ["_".join(c) for c in agg.columns]

    fovs = fov_list_for_well(well)
    print(f"[{well_name}] processing {len(fovs)} FOVs: {fovs[0]}..{fovs[-1]}", flush=True)

    out_well_dir = os.path.join(OUT_DIR, "Track_Cell_Result", well_name)
    os.makedirs(out_well_dir, exist_ok=True)

    results = []
    for i, n in enumerate(fovs):
        res = process_fov(n, df, agg, well_name, well["treatment"])
        if res is not None:
            out_path = os.path.join(out_well_dir, f"Field{n:04d}_track_cell.csv")
            res.to_csv(out_path, index=False)
            results.append(res)
        if (i + 1) % 15 == 0 or (i + 1) == len(fovs):
            print(f"[{well_name}] {i+1}/{len(fovs)} FOVs done ({time.time()-t0:.1f}s elapsed)", flush=True)

    if results:
        combined = pd.concat(results, ignore_index=True)
    else:
        combined = pd.DataFrame()
    combined_path = os.path.join(OUT_DIR, f"combined_{well_name}.csv")
    combined.to_csv(combined_path, index=False)
    print(f"[{well_name}] DONE in {time.time()-t0:.1f}s -> {len(combined)} track-cell rows -> {combined_path}", flush=True)
    return combined


def main():
    t0 = time.time()
    print("Loading shared cell_boundaries.csv.gz ...", flush=True)
    cell_boundaries_csv = load_cell_boundaries()
    print(f"Loaded {len(cell_boundaries_csv)} boundary vertices, "
          f"{cell_boundaries_csv['cell_id'].nunique()} unique cells "
          f"({time.time()-t0:.1f}s)", flush=True)

    all_combined = []
    for well_name, well in WELLS.items():
        combined = run_well(well_name, well, cell_boundaries_csv)
        all_combined.append(combined)

    master = pd.concat(all_combined, ignore_index=True)
    master_path = os.path.join(OUT_DIR, "combined_all_wells.csv")
    master.to_csv(master_path, index=False)
    print(f"\nALL DONE in {time.time()-t0:.1f}s. Master table: {len(master)} rows -> {master_path}", flush=True)


if __name__ == "__main__":
    main()
