#sgRNA analysis for outlier detection, sensitivity & specificity, and raw counts of well-expected vs unexpected 

import re
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from __future__ import annotations

# Loading files
DATA_DIR = Path('path/to/Xenium_file')
RUN_DIR = DATA_DIR / 'path/to/file'
CELL_FEATURE_MATRIX_H5_PATH = RUN_DIR / 'cell_feature_matrix.h5'
CELLS_PATH = RUN_DIR / 'cells.parquet'
CELL_FEATURE_MATRIX_H5_PATH, CELLS_PATH
adata = sc.read_10x_h5(CELL_FEATURE_MATRIX_H5_PATH)
adata

#gRNAs
sgrna_mask = adata.var['gene_ids'].astype(str).str.contains('-', regex=False)
sgrna_features = adata.var_names[sgrna_mask].tolist()
sgrna_gene_ids = adata.var.loc[sgrna_mask, 'gene_ids'].tolist()
len(sgrna_features), list(zip(sgrna_features[:10], sgrna_gene_ids[:10]))

#Cell IDs with gRNAs
adata_sgrna = adata[:, sgrna_features]
sgrna_matrix = adata_sgrna.X.tocsr()

detected_sgrna_features_per_cell = []
detected_sgrna_gene_ids_per_cell = []
detected_sgrna_counts_per_cell = []

for i in range(sgrna_matrix.shape[0]):
    start = sgrna_matrix.indptr[i]
    end = sgrna_matrix.indptr[i + 1]
    feature_idx = sgrna_matrix.indices[start:end]
    feature_counts = sgrna_matrix.data[start:end].tolist()
    detected_sgrna_features_per_cell.append(adata_sgrna.var_names[feature_idx].tolist())
    detected_sgrna_gene_ids_per_cell.append(adata_sgrna.var.iloc[feature_idx]['gene_ids'].tolist())
    detected_sgrna_counts_per_cell.append(feature_counts)

cell_sgrna_table = pd.DataFrame({
    'cell_id': adata_sgrna.obs_names,
    'detected_sgrna_features': detected_sgrna_features_per_cell,
    'detected_sgrna_gene_ids': detected_sgrna_gene_ids_per_cell,
    'detected_sgrna_counts': detected_sgrna_counts_per_cell,
    'n_sgrnas_detected': np.asarray((adata_sgrna.X > 0).sum(axis=1)).ravel(),
    'total_sgrna_counts': np.asarray(adata_sgrna.X.sum(axis=1)).ravel(),
})

cell_sgrna_table.head()

# Number of gRNA/cell and distinct sgRNA detected per cell
sg_per_cell_dist = cell_sgrna_table['n_sgrnas_detected'].value_counts().sort_index()
sg_per_cell_dist
plt.figure(figsize=(8, 4))
plt.bar(sg_per_cell_dist.index, sg_per_cell_dist.values)
plt.xlabel('Number of distinct sgRNAs detected in a cell')
plt.ylabel('Number of cells')
plt.title('Distinct sgRNAs detected per cell')
plt.xticks(sg_per_cell_dist.index)
plt.tight_layout()
plt.show()

#total sgRNA counts per cell
plt.figure(figsize=(8, 4))
plt.hist(cell_sgrna_table['total_sgrna_counts'], bins=50)
plt.xlabel('Total sgRNA counts per cell')
plt.ylabel('Number of cells')
plt.title('Total sgRNA counts per cell')
plt.tight_layout()
plt.show()

# Loading cells with spatial coordinates
cells = pd.read_parquet(CELLS_PATH)
cells[['cell_id', 'x_centroid', 'y_centroid']].head()

# combining cell IDs with spatial coordinates
cell_with_coords_table = cell_sgrna_table.merge(
    cells[['cell_id', 'x_centroid', 'y_centroid']],
    on='cell_id',
    how='left',
)

cell_with_coords_table.head()

# Coordinate splits used to infer the 4 wells from cell centroids
X_SPLIT = 5000
Y_SPLIT = 11000

cell_well_table = cell_with_coords_table.copy()
cell_well_table['well_col'] = np.where(cell_well_table['x_centroid'] < X_SPLIT, 'col_1', 'col_2')
cell_well_table['well_row'] = np.where(cell_well_table['y_centroid'] < Y_SPLIT, 'row_1', 'row_2')
cell_well_table['inferred_well'] = cell_well_table['well_row'] + '_' + cell_well_table['well_col']
cell_well_table[['cell_id', 'x_centroid', 'y_centroid', 'inferred_well']].head()

# Number of cells from each well
well_order = ['row_1_col_1', 'row_1_col_2', 'row_2_col_1', 'row_2_col_2']
cell_count_by_part = (
    cell_well_table['inferred_well']
    .value_counts()
    .reindex(well_order, fill_value=0)
    .rename_axis('part')
    .reset_index(name='n_cells')
)
cell_count_by_part

# Flatten the sparse cell x guide count matrix into one row per detection, with well labels attached
sgrna_counts = np.asarray(sgrna_matrix[cell_idx, sgrna_idx]).ravel()

cell_sgrna_long = pd.DataFrame({
    'cell_id': adata_sgrna.obs_names[cell_idx],
    'sgrna_feature': adata_sgrna.var_names[sgrna_idx],
    'sgrna_gene_id': adata_sgrna.var.iloc[sgrna_idx]['gene_ids'].to_numpy(),
    'sgrna_count_in_cell': sgrna_counts,
})

cell_sgrna_long['sgrna_gene_id_normalized'] = cell_sgrna_long['sgrna_gene_id'].map(
    lambda x: re.sub(r'(?<=[A-Za-z])-(?=\d)', '', str(x).strip())
)

cell_sgrna_long = cell_sgrna_long.merge(
    cell_well_table[['cell_id', 'inferred_well']],
    on='cell_id',
    how='left',
)

cell_sgrna_long.head()

#gRNAs and outliers from previous experiment
pool1_sgrnas_raw = ['ONECUT2-1', 'ONECUT2-2', 'BNC2-1', 'BNC2-2', 'EBF1-1', 'EBF1-2', 'ID4-1', 'ID4-2', 'EBF3-1', 'EBF3-2', 'ID2-1', 'ID2-2', 
                    'CARHSP1-1', 'CARHSP1-2', 'LHX9-1', 'LHX9-2', 'PROX1-1', 'PROX1-2', 'BCL6-1', 'BCL6-2', 'RUNX1-1', 'RUNX1-2', 'CREB5-1', 
                    'CREB5-2', 'ISL1-1', 'ISL1-2', 'SOX5-3', 'SOX5-4', 'SOX5-5', 'ONECUT1-3', 'ONECUT1-4', 'ONECUT1-5', 'NFIB-3', 'NFIB-4', 'NFIB-5', 
                    'JUNB-3', 'JUNB-4', 'JUNB-5', 'NFE2L2-3', 'NFE2L2-4', 'NFE2L2-5', 'KLF3-3', 'KLF3-4', 'KLF3-5', 'HMGB2-3', 'HMGB2-4', 'HMGB2-5', 'MEIS2-3', 
                    'MEIS2-4', 'MEIS2-5', 'AFF1-3', 'AFF1-4', 'AFF1-5', 'ZBTB7C-3', 'ZBTB7C-4', 'ZBTB7C-5', 'SIX1-3', 'SIX1-4', 'SIX1-5', 'Negative control-3', 
                    'Negative control-4', 'Negative control-5']
pool2_sgrnas_raw = ['ONECUT2-3', 'ONECUT2-4', 'ONECUT2-5', 'BNC2-3', 'BNC2-4', 'BNC2-5', 'EBF1-3', 'EBF1-4', 'EBF1-5', 'ID4-3', 'ID4-4', 'ID4-5', 
                    'EBF3-3', 'EBF3-4', 'EBF3-5', 'ID2-3', 'ID2-4', 'ID2-5', 'CARHSP1-3', 'CARHSP1-4', 'CARHSP1-5', 'LHX9-3', 'LHX9-4', 'LHX9-5', 'PROX1-3', 
                    'PROX1-4', 'PROX1-5', 'BCL6-3', 'BCL6-4', 'BCL6-5', 'RUNX1-3', 'RUNX1-4', 'RUNX1-5', 'CREB5-3', 'CREB5-4', 'CREB5-5', 'ISL1-3', 'ISL1-4', 
                    'ISL1-5', 'SOX5-1', 'SOX5-2', 'ONECUT1-1', 'ONECUT1-2', 'NFIB-1', 'NFIB-2', 'JUNB-1', 'JUNB-2', 'NFE2L2-1', 'NFE2L2-2', 'KLF3-1', 'KLF3-2', 'HMGB2-1', 
                    'HMGB2-2', 'MEIS2-1', 'MEIS2-2', 'AFF1-1', 'AFF1-2', 'ZBTB7C-1', 'ZBTB7C-2', 'SIX1-1', 'SIX1-2', 'Negative control-1', 'Negative control-2']
outlier_sgrnas_raw = ['SOX5-3', 'ID2-4', 'NFIB-2', 'ZBTB7C-1', 'ID2-2']

# Creating a table with cell info and gRNA
def normalize_sgrna_name(name):
    name = str(name).strip()
    negative_control_match = re.match(r'(?i)^negative[ _-]*(control|ctrl)[-_ ]*(\d+)$', name)
    if negative_control_match:
        return f"Negative_control-{negative_control_match.group(2)}"
    return re.sub(r'(?<=[A-Za-z])-(?=\d)', '', name)

cell_sgrna_long['sgrna_gene_id_normalized'] = cell_sgrna_long['sgrna_gene_id'].map(normalize_sgrna_name)

pool_annotation_table = pd.concat([
    pd.DataFrame({'sgrna_from_pool_file': pool1_sgrnas_raw, 'expected_pool': 'pool_1'}),
    pd.DataFrame({'sgrna_from_pool_file': pool2_sgrnas_raw, 'expected_pool': 'pool_2'}),
], ignore_index=True)
pool_annotation_table['sgrna_gene_id_normalized'] = pool_annotation_table['sgrna_from_pool_file'].map(normalize_sgrna_name)
pool_annotation_table = pool_annotation_table.drop_duplicates('sgrna_gene_id_normalized', keep='first')

outlier_annotation_table = pd.DataFrame({'requested_outlier': outlier_sgrnas_raw})
outlier_annotation_table['sgrna_gene_id_normalized'] = outlier_annotation_table['requested_outlier'].map(normalize_sgrna_name)
outlier_sgrnas = set(outlier_annotation_table['sgrna_gene_id_normalized'])

sgrna_pool_table = (
    cell_sgrna_long[['sgrna_gene_id', 'sgrna_gene_id_normalized']]
    .drop_duplicates()
    .sort_values('sgrna_gene_id')
    .reset_index(drop=True)
)
sgrna_pool_table = sgrna_pool_table.merge(
    pool_annotation_table[['sgrna_gene_id_normalized', 'expected_pool']],
    on='sgrna_gene_id_normalized',
    how='left',
)
sgrna_pool_table['expected_pool'] = sgrna_pool_table['expected_pool'].fillna('unassigned')
sgrna_pool_table['is_outlier'] = sgrna_pool_table['sgrna_gene_id_normalized'].isin(outlier_sgrnas)

requested_outlier_check = outlier_annotation_table.merge(
    sgrna_pool_table[['sgrna_gene_id', 'sgrna_gene_id_normalized', 'expected_pool', 'is_outlier']],
    on='sgrna_gene_id_normalized',
    how='left',
)

cell_sgrna_long = cell_sgrna_long.drop(columns=['expected_pool', 'is_outlier'], errors='ignore')
cell_sgrna_long = cell_sgrna_long.merge(
    sgrna_pool_table[['sgrna_gene_id', 'expected_pool', 'is_outlier']],
    on='sgrna_gene_id',
    how='left',
)
cell_sgrna_long['expected_pool'] = cell_sgrna_long['expected_pool'].fillna('unassigned')
cell_sgrna_long['is_outlier'] = cell_sgrna_long['is_outlier'].fillna(False)

cell_sgrna_long_no_outliers = cell_sgrna_long[~cell_sgrna_long['is_outlier']].copy()

display(requested_outlier_check)
display(sgrna_pool_table[sgrna_pool_table['sgrna_gene_id'].str.contains('NFIB|ID2|SOX5', regex=True)].sort_values('sgrna_gene_id'))
display(sgrna_pool_table[sgrna_pool_table['sgrna_gene_id'].str.contains('Negative', case=False, regex=True)].sort_values('sgrna_gene_id'))
display(sgrna_pool_table[sgrna_pool_table['is_outlier']].sort_values('sgrna_gene_id'))
sgrna_pool_table['is_outlier'].value_counts()


# Compare two wells guide-by-guide: cells detected in each, their percentage split, and log2 fold change.
def make_pairwise_well_summary(long_df, well_a, well_b):
    pair_df = long_df[long_df['inferred_well'].isin([well_a, well_b])].copy()

    summary = (
        pair_df
        .groupby(['sgrna_gene_id', 'inferred_well'])
        .agg(
            n_cells_detected=('cell_id', 'nunique'),
            total_counts=('sgrna_count_in_cell', 'sum'),
        )
        .reset_index()
    )

    wide = (
        summary
        .pivot(index='sgrna_gene_id', columns='inferred_well', values=['n_cells_detected', 'total_counts'])
        .fillna(0)
    )

    wide.columns = [f'{metric}_{well}' for metric, well in wide.columns]
    wide = wide.reset_index()

    for col in [
        f'n_cells_detected_{well_a}', f'n_cells_detected_{well_b}',
        f'total_counts_{well_a}', f'total_counts_{well_b}',
    ]:
        if col not in wide.columns:
            wide[col] = 0

    wide['total_cells_detected'] = wide[f'n_cells_detected_{well_a}'] + wide[f'n_cells_detected_{well_b}']
    wide[f'pct_{well_a}'] = np.where(
        wide['total_cells_detected'] > 0,
        100 * wide[f'n_cells_detected_{well_a}'] / wide['total_cells_detected'],
        0,
    )
    wide[f'pct_{well_b}'] = np.where(
        wide['total_cells_detected'] > 0,
        100 * wide[f'n_cells_detected_{well_b}'] / wide['total_cells_detected'],
        0,
    )
    wide[f'log2_fc_{well_a}_vs_{well_b}'] = np.log2(
        (wide[f'n_cells_detected_{well_a}'] + 1) / (wide[f'n_cells_detected_{well_b}'] + 1)
    )

    return wide.sort_values('total_cells_detected', ascending=False)

# Run all four well pairings, two left-right (across pools) and two top-bottom (across rows), with and without outlier guides.
row1_left_vs_right = make_pairwise_well_summary(cell_sgrna_long, 'row_1_col_1', 'row_1_col_2')
col1_top_vs_bottom = make_pairwise_well_summary(cell_sgrna_long, 'row_1_col_1', 'row_2_col_1')
row2_left_vs_right = make_pairwise_well_summary(cell_sgrna_long, 'row_2_col_1', 'row_2_col_2')
col2_top_vs_bottom = make_pairwise_well_summary(cell_sgrna_long, 'row_1_col_2', 'row_2_col_2')

row1_left_vs_right_no_outliers = make_pairwise_well_summary(cell_sgrna_long_no_outliers, 'row_1_col_1', 'row_1_col_2')
col1_top_vs_bottom_no_outliers = make_pairwise_well_summary(cell_sgrna_long_no_outliers, 'row_1_col_1', 'row_2_col_1')
row2_left_vs_right_no_outliers = make_pairwise_well_summary(cell_sgrna_long_no_outliers, 'row_2_col_1', 'row_2_col_2')
col2_top_vs_bottom_no_outliers = make_pairwise_well_summary(cell_sgrna_long_no_outliers, 'row_1_col_2', 'row_2_col_2')

row1_left_vs_right.head(20)

# Per-cell sensitivity/specificity of pool assignment, scored strictly (ambiguous cells dropped) and permissively.
well_pool_map = {
    'row_1_col_1': 'pool_1',
    'row_1_col_2': 'pool_2',
    'row_2_col_1': 'pool_1',
    'row_2_col_2': 'pool_2',
}
well_treatment_map = {
    'row_1_col_1': 'row_1',
    'row_1_col_2': 'row_1',
    'row_2_col_1': 'row_2',
    'row_2_col_2': 'row_2',
}
cell_pool_detection_table = cell_well_table[[
    'cell_id', 'inferred_well', 'x_centroid', 'y_centroid'
]].copy()
cell_pool_detection_table['expected_pool_of_well'] = cell_pool_detection_table['inferred_well'].map(well_pool_map)
cell_pool_detection_table['treatment_row'] = cell_pool_detection_table['inferred_well'].map(well_treatment_map)

valid_cell_summary = (
    cell_sgrna_long_no_outliers
    .groupby('cell_id')
    .agg(
        total_valid_sgrna_counts=('sgrna_count_in_cell', 'sum'),
        n_valid_sgrnas=('sgrna_gene_id_normalized', 'nunique'),
        has_pool1=('expected_pool', lambda x: (x == 'pool_1').any()),
        has_pool2=('expected_pool', lambda x: (x == 'pool_2').any()),
    )
    .reset_index()
)

cell_pool_detection_table = cell_pool_detection_table.merge(
    valid_cell_summary,
    on='cell_id',
    how='left',
)

cell_pool_detection_table['total_valid_sgrna_counts'] = cell_pool_detection_table['total_valid_sgrna_counts'].fillna(0)
cell_pool_detection_table['n_valid_sgrnas'] = cell_pool_detection_table['n_valid_sgrnas'].fillna(0)
cell_pool_detection_table['has_pool1'] = cell_pool_detection_table['has_pool1'].fillna(False)
cell_pool_detection_table['has_pool2'] = cell_pool_detection_table['has_pool2'].fillna(False)

cell_pool_detection_table['is_bfp_positive_proxy'] = cell_pool_detection_table['total_valid_sgrna_counts'] > 0
cell_pool_detection_table['is_ambiguous_pool_detection'] = (
    cell_pool_detection_table['has_pool1'] & cell_pool_detection_table['has_pool2']
)
cell_pool_detection_table['pool_detection_class'] = np.select(
    [
        cell_pool_detection_table['has_pool1'] & ~cell_pool_detection_table['has_pool2'],
        ~cell_pool_detection_table['has_pool1'] & cell_pool_detection_table['has_pool2'],
        cell_pool_detection_table['has_pool1'] & cell_pool_detection_table['has_pool2'],
    ],
    ['pool1_only', 'pool2_only', 'ambiguous_both'],
    default='no_pool_detected',
)


def compute_pool_metrics(cell_df, target_pool, treatment_row, mode='strict'):
    row_df = cell_df[
        (cell_df['treatment_row'] == treatment_row) &
        (cell_df['is_bfp_positive_proxy'])
    ].copy()

    n_bfp_proxy_before_filter = len(row_df)
    n_ambiguous_before_filter = int(row_df['is_ambiguous_pool_detection'].sum())

    if mode == 'strict':
        row_df = row_df[~row_df['is_ambiguous_pool_detection']].copy()

    detected_col = 'has_pool1' if target_pool == 'pool_1' else 'has_pool2'
    expected_positive = row_df['expected_pool_of_well'] == target_pool
    detected_positive = row_df[detected_col]

    tp = int((expected_positive & detected_positive).sum())
    fn = int((expected_positive & ~detected_positive).sum())
    fp = int((~expected_positive & detected_positive).sum())
    tn = int((~expected_positive & ~detected_positive).sum())

    sensitivity = tp / (tp + fn) if (tp + fn) > 0 else np.nan
    specificity = tn / (tn + fp) if (tn + fp) > 0 else np.nan

    return {
        'target_pool': target_pool,
        'treatment_row': treatment_row,
        'mode': mode,
        'n_bfp_proxy_before_filter': n_bfp_proxy_before_filter,
        'n_ambiguous_before_filter': n_ambiguous_before_filter,
        'n_cells_used': len(row_df),
        'TP': tp,
        'FN': fn,
        'FP': fp,
        'TN': tn,
        'sensitivity': sensitivity,
        'specificity': specificity,
    }


pool_detection_metrics = pd.DataFrame([
    compute_pool_metrics(cell_pool_detection_table, target_pool, treatment_row, mode)
    for target_pool in ['pool_1', 'pool_2']
    for treatment_row in ['row_1', 'row_2']
    for mode in ['strict', 'permissive']
])

pool1_detection_metrics = pool_detection_metrics[pool_detection_metrics['target_pool'] == 'pool_1'].copy()

display(cell_pool_detection_table.head())
display(pool1_detection_metrics)
display(pool_detection_metrics)

# Raw-count pool assignment accuracy: contingency table of guide pool x well pool, per treatment row and overall, with and without outlier guides
def add_well_pool_annotations(df):
    out = df.copy()
    out['well_pool'] = out['inferred_well'].map({
        'row_1_col_1': 'pool_1',
        'row_1_col_2': 'pool_2',
        'row_2_col_1': 'pool_1',
        'row_2_col_2': 'pool_2',
    })
    out['treatment_row'] = out['inferred_well'].map({
        'row_1_col_1': 'row_1',
        'row_1_col_2': 'row_1',
        'row_2_col_1': 'row_2',
        'row_2_col_2': 'row_2',
    })
    return out

def make_raw_count_metrics(df, group_label, filter_status):
    summary = (
        df.groupby(['expected_pool', 'well_pool'])
        .agg(raw_sgrna_counts=('sgrna_count_in_cell', 'sum'))
        .reset_index()
    )

    matrix = (
        summary.pivot(index='expected_pool', columns='well_pool', values='raw_sgrna_counts')
        .reindex(index=['pool_1', 'pool_2'], columns=['pool_1', 'pool_2'])
        .fillna(0)
    )

    p1_in_p1 = matrix.loc['pool_1', 'pool_1']
    p1_in_p2 = matrix.loc['pool_1', 'pool_2']
    p2_in_p1 = matrix.loc['pool_2', 'pool_1']
    p2_in_p2 = matrix.loc['pool_2', 'pool_2']

    expected_raw_counts = p1_in_p1 + p2_in_p2
    unexpected_raw_counts = p1_in_p2 + p2_in_p1
    total_raw_counts = expected_raw_counts + unexpected_raw_counts

    records = [
        {
            'filter_status': filter_status,
            'group': group_label,
            'positive_class': 'pool_1',
            'TP_like': p1_in_p1,
            'FN_like': p1_in_p2,
            'FP_like': p2_in_p1,
            'TN_like': p2_in_p2,
            'sensitivity': p1_in_p1 / (p1_in_p1 + p1_in_p2) if (p1_in_p1 + p1_in_p2) > 0 else np.nan,
            'specificity': p2_in_p2 / (p2_in_p1 + p2_in_p2) if (p2_in_p1 + p2_in_p2) > 0 else np.nan,
            'expected_raw_counts': expected_raw_counts,
            'unexpected_raw_counts': unexpected_raw_counts,
            'total_raw_counts': total_raw_counts,
            'on_target_fraction': expected_raw_counts / total_raw_counts if total_raw_counts > 0 else np.nan,
            'cross_pool_leakage_fraction': unexpected_raw_counts / total_raw_counts if total_raw_counts > 0 else np.nan,
        },
        {
            'filter_status': filter_status,
            'group': group_label,
            'positive_class': 'pool_2',
            'TP_like': p2_in_p2,
            'FN_like': p2_in_p1,
            'FP_like': p1_in_p2,
            'TN_like': p1_in_p1,
            'sensitivity': p2_in_p2 / (p2_in_p1 + p2_in_p2) if (p2_in_p1 + p2_in_p2) > 0 else np.nan,
            'specificity': p1_in_p1 / (p1_in_p1 + p1_in_p2) if (p1_in_p1 + p1_in_p2) > 0 else np.nan,
            'expected_raw_counts': expected_raw_counts,
            'unexpected_raw_counts': unexpected_raw_counts,
            'total_raw_counts': total_raw_counts,
            'on_target_fraction': expected_raw_counts / total_raw_counts if total_raw_counts > 0 else np.nan,
            'cross_pool_leakage_fraction': unexpected_raw_counts / total_raw_counts if total_raw_counts > 0 else np.nan,
        },
    ]

    return matrix, pd.DataFrame(records)

raw_count_inputs = {
    'with_outliers': add_well_pool_annotations(cell_sgrna_long),
    'without_outliers': add_well_pool_annotations(cell_sgrna_long_no_outliers),
}

raw_count_confusion_matrices = {}
raw_count_metrics_parts = []

for filter_status, df in raw_count_inputs.items():
    raw_count_confusion_matrices[filter_status] = {}
    for group_label, sub_df in [
        ('row_1', df[df['treatment_row'] == 'row_1']),
        ('row_2', df[df['treatment_row'] == 'row_2']),
        ('overall', df),
    ]:
        matrix, metrics = make_raw_count_metrics(sub_df, group_label, filter_status)
        raw_count_confusion_matrices[filter_status][group_label] = matrix
        raw_count_metrics_parts.append(metrics)

raw_count_metrics = pd.concat(raw_count_metrics_parts, ignore_index=True)

for filter_status in ['without_outliers', 'with_outliers']:
    print(f'Raw-count pool x well tables: {filter_status}')
    for group_label in ['row_1', 'row_2', 'overall']:
        print(group_label)
        display(raw_count_confusion_matrices[filter_status][group_label])

raw_count_metrics_summary = raw_count_metrics[[
    'filter_status',
    'group',
    'positive_class',
    'TP_like',
    'FN_like',
    'FP_like',
    'TN_like',
    'sensitivity',
    'specificity',
    'expected_raw_counts',
    'unexpected_raw_counts',
    'on_target_fraction',
    'cross_pool_leakage_fraction',
]].copy()

display(raw_count_metrics_summary)

plot_df = raw_count_metrics[(raw_count_metrics['group'] == 'overall')].copy()
plot_df['label'] = plot_df['filter_status'] + ' | ' + plot_df['positive_class']

fig, axes = plt.subplots(1, 2, figsize=(12, 4))

x = np.arange(len(plot_df))
width = 0.35
axes[0].bar(x - width / 2, plot_df['sensitivity'], width=width, label='Sensitivity', color='#1b9e77')
axes[0].bar(x + width / 2, plot_df['specificity'], width=width, label='Specificity', color='#d95f02')
axes[0].set_xticks(x)
axes[0].set_xticklabels(plot_df['label'], rotation=30, ha='right')
axes[0].set_ylim(0, 1.05)
axes[0].set_ylabel('Fraction')
axes[0].set_title('Overall raw-count sensitivity and specificity')
axes[0].legend()

x2 = np.arange(len(plot_df))
axes[1].bar(x2, plot_df['expected_raw_counts'], label='Expected raw counts', color='#4daf4a')
axes[1].bar(x2, plot_df['unexpected_raw_counts'], bottom=plot_df['expected_raw_counts'], label='Unexpected raw counts', color='#e41a1c')
axes[1].set_xticks(x2)
axes[1].set_xticklabels(plot_df['label'], rotation=30, ha='right')
axes[1].set_ylabel('Raw sgRNA transcript counts')
axes[1].set_title('Overall expected vs unexpected raw counts')
axes[1].legend()

plt.tight_layout()
plt.show()

row_plot_df = raw_count_metrics[(raw_count_metrics['positive_class'] == 'pool_1')].copy()
row_plot_df['label'] = row_plot_df['filter_status'] + ' | ' + row_plot_df['group']

plt.figure(figsize=(10, 4))
x3 = np.arange(len(row_plot_df))
plt.bar(x3 - width / 2, row_plot_df['sensitivity'], width=width, label='Pool 1 sensitivity', color='#377eb8')
plt.bar(x3 + width / 2, row_plot_df['specificity'], width=width, label='Pool 1 specificity', color='#984ea3')
plt.xticks(x3, row_plot_df['label'], rotation=30, ha='right')
plt.ylim(0, 1.05)
plt.ylabel('Fraction')
plt.title('Pool 1 raw-count sensitivity and specificity by treatment row')
plt.legend()
plt.tight_layout()
plt.show()

# Per-guide off-target table: expected vs wrong-pool raw counts and detected cells, normalized by cells per pool
WELL_TO_POOL = {
    "row_1_col_1": "pool_1",
    "row_1_col_2": "pool_2",
    "row_2_col_1": "pool_1",
    "row_2_col_2": "pool_2",
}

OTHER_POOL = {"pool_1": "pool_2", "pool_2": "pool_1"}

def build_per_guide_offtarget_table(
    cell_sgrna_long: pd.DataFrame,
    cell_well_table: pd.DataFrame,
    out_csv: str | Path | None = None,
) -> pd.DataFrame:
    detections = cell_sgrna_long.copy()
    detections["well_pool"] = detections["inferred_well"].map(WELL_TO_POOL)

    unassigned = sorted(
        detections.loc[~detections["expected_pool"].isin(OTHER_POOL), "sgrna_gene_id"].unique()
    )
    if unassigned:
        print(f"Dropping {len(unassigned)} guide(s) with no pool assignment: {unassigned}")
        detections = detections[detections["expected_pool"].isin(OTHER_POOL)]

    cells_per_pool = (
        cell_well_table.assign(well_pool=cell_well_table["inferred_well"].map(WELL_TO_POOL))
        .groupby("well_pool")["cell_id"]
        .nunique()
    )

    keys = ["sgrna_gene_id", "sgrna_gene_id_normalized", "expected_pool", "is_outlier"]
    wide = (
        detections.groupby(keys + ["well_pool"])
        .agg(raw_counts=("sgrna_count_in_cell", "sum"),
             detected_cells=("cell_id", "nunique"))
        .unstack("well_pool")
        .reindex(columns=pd.MultiIndex.from_product(
            [["raw_counts", "detected_cells"], ["pool_1", "pool_2"]]
        ))
        .fillna(0)
    )
    wide.columns = [f"{metric}_{pool}" for metric, pool in wide.columns]
    wide = wide.reset_index()
    expected_pool = wide["expected_pool"]
    unexpected_pool = expected_pool.map(OTHER_POOL)

    def by_pool(metric: str, pool: pd.Series) -> np.ndarray:
        return np.where(pool == "pool_1", wide[f"{metric}_pool_1"], wide[f"{metric}_pool_2"])

    wide["expected_raw_counts"] = by_pool("raw_counts", expected_pool)
    wide["unexpected_raw_counts"] = by_pool("raw_counts", unexpected_pool)
    wide["expected_detected_cells"] = by_pool("detected_cells", expected_pool)
    wide["unexpected_detected_cells"] = by_pool("detected_cells", unexpected_pool)

    wide["total_raw_counts"] = wide["expected_raw_counts"] + wide["unexpected_raw_counts"]
    wide["unexpected_count_fraction"] = (
        wide["unexpected_raw_counts"] / wide["total_raw_counts"].replace(0, np.nan)
    )

    wide["total_cells_in_expected_wells"] = expected_pool.map(cells_per_pool)
    wide["total_cells_in_unexpected_wells"] = unexpected_pool.map(cells_per_pool)

    wide["expected_rate_per_cell"] = (
        wide["expected_raw_counts"] / wide["total_cells_in_expected_wells"]
    )
    wide["unexpected_rate_per_cell"] = (
        wide["unexpected_raw_counts"] / wide["total_cells_in_unexpected_wells"]
    )
    wide["unexpected_detected_cell_fraction"] = (
        wide["unexpected_detected_cells"] / wide["total_cells_in_unexpected_wells"]
    )

    wide["is_negative_control"] = (
        wide["sgrna_gene_id_normalized"]
        .astype(str)
        .str.match(r"(?i)^negative[ _-]*control", na=False)
    )

    table = wide.sort_values(
        ["unexpected_rate_per_cell", "unexpected_raw_counts"], ascending=False
    ).reset_index(drop=True)

    if out_csv is not None:
        table.to_csv(out_csv, index=False)
        print(f"Saved per-guide off-target table: {out_csv}")
    return table

SUMMARY_COLUMNS = [
    "sgrna_gene_id",
    "expected_pool",
    "is_outlier",
    "is_negative_control",
    "expected_raw_counts",
    "unexpected_raw_counts",
    "total_raw_counts",
    "unexpected_count_fraction",
    "expected_detected_cells",
    "unexpected_detected_cells",
    "expected_rate_per_cell",
    "unexpected_rate_per_cell",
    "unexpected_detected_cell_fraction",
]

# Generating as an sgRNA off-target count CSV file
per_sgrna_leakage_metrics = build_per_guide_offtarget_table(
    cell_sgrna_long, cell_well_table, out_csv="per_sgrna_offtarget_counts.csv"
)
per_sgrna_leakage_metrics[SUMMARY_COLUMNS].head(10)

# Generating overall sensitivity & specificity, expected vs unexpected counts and off-target counts per guide
def plot_combined_overall_on_target_off_target_summary(
    raw_count_metrics,
    per_sgrna_leakage_metrics,
    guide_filter_status='with_outliers',
    guide_lo_max=None,      
    guide_hi_min=None,      
    guide_width_ratios=(3, 1),
    guide_label_size=4,
):
    overall_df = (
        raw_count_metrics.loc[
            (raw_count_metrics['group'] == 'overall')
            & (raw_count_metrics['positive_class'].isin(['pool_1', 'pool_2']))
            & (raw_count_metrics['filter_status'].isin([
                'with_outliers', 'without_outliers',
            ]))
        ]
        .copy()
    )
    filter_order = ['with_outliers', 'without_outliers']
    filter_labels = {
        'with_outliers': 'With outliers',
        'without_outliers': 'Without outliers',
    }
    pool_labels = {'pool_1': 'Pool 1', 'pool_2': 'Pool 2'}

    overall_df['filter_status'] = pd.Categorical(
        overall_df['filter_status'], categories=filter_order, ordered=True
    )
    overall_df['positive_class'] = pd.Categorical(
        overall_df['positive_class'], categories=['pool_1', 'pool_2'], ordered=True
    )
    overall_df = overall_df.sort_values(
        ['positive_class', 'filter_status']
    ).reset_index(drop=True)
    if len(overall_df) != 4:
        raise RuntimeError(
            'Expected overall rows for pool_1/pool_2 under both outlier filter statuses.'
        )

    guide_df = per_sgrna_leakage_metrics.copy()
    if guide_filter_status == 'without_outliers':
        guide_df = guide_df.loc[~guide_df['is_outlier']].copy()
    guide_df = guide_df.sort_values(
        ['unexpected_raw_counts', 'sgrna_gene_id'], ascending=[True, True]
    ).copy()
    guide_df['plot_label'] = np.where(
        guide_df['is_outlier'],
        guide_df['sgrna_gene_id'] + ' *',
        guide_df['sgrna_gene_id'],
    )

    fig = plt.figure(figsize=(16, 8))
    outer = fig.add_gridspec(1, 3, width_ratios=[1.25, 1.25, 1.65], wspace=0.42)
    sens_ax = fig.add_subplot(outer[0])
    count_ax = fig.add_subplot(outer[1])
    guide_gs = outer[2].subgridspec(1, 2, width_ratios=list(guide_width_ratios), wspace=0.05)
    guide_lo = fig.add_subplot(guide_gs[0])
    guide_hi = fig.add_subplot(guide_gs[1], sharey=guide_lo)

    sens_df = overall_df.loc[
        overall_df['positive_class'].astype(str).eq('pool_1')
    ].copy()
    sens_df = sens_df.sort_values('filter_status').reset_index(drop=True)
    sens_plot_labels = [
        filter_labels[str(row.filter_status)] for row in sens_df.itertuples()
    ]
    well_count_rows = []
    for row in sens_df.itertuples():
        well_count_rows.extend([
            {'filter_status': row.filter_status, 'well_pool': 'pool_1',
             'expected_raw_counts': row.TP_like, 'unexpected_raw_counts': row.FP_like},
            {'filter_status': row.filter_status, 'well_pool': 'pool_2',
             'expected_raw_counts': row.TN_like, 'unexpected_raw_counts': row.FN_like},
        ])
    well_count_df = pd.DataFrame(well_count_rows)
    well_count_df['filter_status'] = pd.Categorical(
        well_count_df['filter_status'], categories=filter_order, ordered=True
    )
    well_count_df['well_pool'] = pd.Categorical(
        well_count_df['well_pool'], categories=['pool_1', 'pool_2'], ordered=True
    )
    well_count_df = well_count_df.sort_values(
        ['well_pool', 'filter_status']
    ).reset_index(drop=True)
    count_plot_labels = [
        f'{pool_labels[str(row.well_pool)]} wells\n{filter_labels[str(row.filter_status)]}'
        for row in well_count_df.itertuples()
    ]

    x_sens = np.arange(len(sens_df))
    x_count = np.arange(len(well_count_df))
    width = 0.34

    sens_bars = [
        sens_ax.bar(x_sens - width / 2, sens_df['sensitivity'], width=width,
                    label='Sensitivity', color='#1b9e77', edgecolor='none'),
        sens_ax.bar(x_sens + width / 2, sens_df['specificity'], width=width,
                    label='Specificity', color='#d95f02', edgecolor='none'),
    ]
    sens_ax.set_xticks(x_sens)
    sens_ax.set_xticklabels(sens_plot_labels, rotation=20, ha='right')
    sens_ax.set_ylim(0, 1.05)
    sens_ax.set_ylabel('Fraction')
    sens_ax.set_title('Overall raw-count\nsensitivity and specificity')
    sens_ax.legend(loc='upper left', bbox_to_anchor=(0.1, 1.0), fontsize=8, frameon=False)
    for bars in sens_bars:
        sens_ax.bar_label(bars, labels=[f'{b.get_height():.1%}' for b in bars],
                          padding=3, fontsize=8, color='black')

    expected_bars = count_ax.bar(
        x_count, well_count_df['expected_raw_counts'], label='Expected raw counts',
        color='#4daf4a', edgecolor='none', width=0.56)
    unexpected_bars = count_ax.bar(
        x_count, well_count_df['unexpected_raw_counts'],
        bottom=well_count_df['expected_raw_counts'],
        label='Unexpected raw counts', color='#e41a1c', edgecolor='none', width=0.56)
    count_ax.set_xticks(x_count)
    count_ax.set_xticklabels(count_plot_labels, rotation=25, ha='right')
    count_ax.set_ylabel('Raw sgRNA transcript counts')
    count_ax.set_title('Well-specific expected vs\nunexpected raw counts')
    count_ax.legend(fontsize=8, frameon=False)
    count_ax.bar_label(expected_bars,
                       labels=[f'{int(v):,}' for v in well_count_df['expected_raw_counts']],
                       label_type='center', fontsize=8, color='black', fontweight='bold')
    count_ax.bar_label(unexpected_bars,
                       labels=[f'{int(v):,}' for v in well_count_df['unexpected_raw_counts']],
                       label_type='center', fontsize=8, color='black', fontweight='bold')

    guide_values = guide_df['unexpected_raw_counts'].to_numpy(dtype=float)
    labels = guide_df['plot_label'].to_numpy()
    y = np.arange(len(guide_values))

    if guide_lo_max is None:
        bulk = guide_values[~guide_df['is_outlier'].to_numpy(dtype=bool)]
        bulk = bulk[np.isfinite(bulk)]
        guide_lo_max = float(bulk.max()) * 1.15 if len(bulk) and bulk.max() > 0 \
            else float(np.nanmax(guide_values))
    above = guide_values[guide_values > guide_lo_max]
    if guide_hi_min is None:
        guide_hi_min = float(above.min()) * 0.90 if len(above) else guide_lo_max * 1.5
    guide_hi_max = float(np.nanmax(guide_values)) * 1.05 if len(above) else guide_hi_min * 1.1

    for ax in (guide_lo, guide_hi):
        ax.barh(y, guide_values, color='#9c755f', edgecolor='none', height=0.55)

    guide_lo.set_xlim(0, guide_lo_max)
    guide_hi.set_xlim(guide_hi_min, guide_hi_max)
    guide_lo.set_yticks(y)
    guide_lo.set_yticklabels(labels)
    guide_lo.set_ylim(-0.8, len(y) - 0.2)
    guide_lo.tick_params(axis='y', labelsize=guide_label_size, pad=1)
    guide_lo.tick_params(axis='x', labelsize=8)
    guide_hi.tick_params(axis='y', left=False, labelleft=False)
    guide_hi.tick_params(axis='x', labelsize=8)

    for yi, val in zip(y, guide_values):
        if val >= guide_hi_min:
            guide_hi.text(val, yi, f'  {int(val):,}', va='center', ha='left',
                          fontsize=7, fontweight='bold', color='black', clip_on=False)

    guide_lo.spines['right'].set_visible(False)
    guide_hi.spines['left'].set_visible(False)
    guide_lo.spines['top'].set_visible(False)
    guide_hi.spines['top'].set_visible(False)
    d = 0.5
    break_kw = dict(marker=[(-1, -d), (1, d)], markersize=10, linestyle='none',
                    color='k', mec='k', mew=1, clip_on=False)
    guide_lo.plot([1, 1], [0, 1], transform=guide_lo.transAxes, **break_kw)
    guide_hi.plot([0, 0], [0, 1], transform=guide_hi.transAxes, **break_kw)

    guide_lo.set_xlabel('Off-target raw counts')
    guide_lo.xaxis.set_label_coords(0.5 + (1 - guide_width_ratios[0] / sum(guide_width_ratios)), -0.06)
    guide_lo.set_title('Off-target counts per guide\n'
                       f'({filter_labels.get(guide_filter_status, guide_filter_status)})',
                       x=0.5 + (1 - guide_width_ratios[0] / sum(guide_width_ratios)))
    if guide_filter_status == 'with_outliers':
        guide_hi.text(0.98, 0.01, '* outlier', transform=guide_hi.transAxes,
                      ha='right', va='bottom', fontsize=7)

    for axis in (sens_ax, count_ax):
        axis.spines['top'].set_visible(False)
        axis.spines['right'].set_visible(False)

    fig.suptitle('Overall gRNAs pool assignment summary',
                 fontsize=14, fontweight='bold', y=1.03)
    return fig, (sens_ax, count_ax, guide_lo, guide_hi)


combined_summary_fig, combined_summary_axes = plot_combined_overall_on_target_off_target_summary(
    raw_count_metrics=raw_count_metrics,
    per_sgrna_leakage_metrics=per_sgrna_leakage_metrics,
    guide_filter_status='with_outliers',
)

combined_summary_svg_path = '/path/file_name.svg'
combined_summary_fig.savefig(combined_summary_svg_path, format='svg', bbox_inches='tight')
print(f'Saved SVG: {combined_summary_svg_path}')
plt.show()
