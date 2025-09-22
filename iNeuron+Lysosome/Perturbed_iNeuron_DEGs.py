print(adata.var['feature_types'].unique())
probe_vars = [v for v in adata_clean.var_names if re.search(r"-\d+$", v)]
print(f"Probe‐style features in adata.var_names: {probe_vars}")

probes_to_remove = [
    'AFF1-1', 'AFF1-2', 'AFF1-3', 'AFF1-4',
    'BCL6-1', 'BCL6-2', 'BCL6-3', 'BCL6-4', 'BCL6-5',
    'BNC2-1', 'BNC2-2', 'BNC2-3', 'BNC2-4', 'BNC2-5',
    'CARHSP1-1', 'CARHSP1-2', 'CARHSP1-3', 'CARHSP1-4', 'CARHSP1-5',
    'CREB5-1', 'CREB5-2', 'CREB5-3', 'CREB5-4', 'CREB5-5',
    'EBF1-1', 'EBF1-2', 'EBF1-3', 'EBF1-4', 'EBF1-5',
    'EBF3-1', 'EBF3-2', 'EBF3-3', 'EBF3-4', 'EBF3-5',
    'HMGB2-1', 'HMGB2-2', 'HMGB2-3', 'HMGB2-4', 'HMGB2-5',
    'ID2-1', 'ID2-2', 'ID2-3', 'ID2-4', 'ID2-5',
    'ID4-1', 'ID4-2', 'ID4-3', 'ID4-4', 'ID4-5',
    'ISL1-1', 'ISL1-2', 'ISL1-3', 'ISL1-4', 'ISL1-5',
    'JUNB-1', 'JUNB-2', 'JUNB-3', 'JUNB-4', 'JUNB-5',
    'KLF3-1', 'KLF3-2', 'KLF3-3', 'KLF3-4', 'KLF3-5',
    'LHX9-1', 'LHX9-3', 'LHX9-4', 'LHX9-5',
    'MEIS2-1', 'MEIS2-2', 'MEIS2-3', 'MEIS2-4', 'MEIS2-5',
    'NFE2L2-1', 'NFE2L2-2', 'NFE2L2-3', 'NFE2L2-4', 'NFE2L2-5',
    'NFIB-1', 'NFIB-2', 'NFIB-3', 'NFIB-4', 'NFIB-5',
    'Negative_control-1', 'Negative_control-2', 'Negative_control-3', 'Negative_control-4', 'Negative_control-5',
    'ONECUT1-1', 'ONECUT1-2', 'ONECUT1-3', 'ONECUT1-4', 'ONECUT1-5',
    'ONECUT2-1', 'ONECUT2-2', 'ONECUT2-3', 'ONECUT2-4', 'ONECUT2-5',
    'PROX1-1', 'PROX1-2', 'PROX1-3', 'PROX1-5',
    'RUNX1-1', 'RUNX1-2', 'RUNX1-3', 'RUNX1-4', 'RUNX1-5',
    'SIX1-1', 'SIX1-2', 'SIX1-3', 'SIX1-4', 'SIX1-5',
    'SOX5-1', 'SOX5-2', 'SOX5-3', 'SOX5-4', 'SOX5-5',
    'ZBTB7C-1', 'ZBTB7C-2', 'ZBTB7C-3', 'ZBTB7C-4', 'ZBTB7C-5'
]
to_remove = [p for p in probes_to_remove if p in adata.var_names]
keep_vars = [v for v in adata.var_names if v not in to_remove]
adata_clean = adata[:, keep_vars].copy()

# DEG
keep_genes = [g for g in adata_clean.var_names if not re.search(r'-\d+$', g)]
adata = adata_clean[:, keep_genes].copy()

exclude   = {'total_counts', 'Negative_control'}
pert_cols = [
    c for c in adata.obs.columns
    if c not in exclude and adata.obs[c].dtype.kind in ('i','b','u')
]
df_label = adata.obs[pert_cols + ['Negative_control']]
group_labels = (
    df_label
      .idxmax(axis=1)
      .replace({'Negative_control': 'Control'})
)

X = adata.X
if sparse.issparse(X):
    X = X.toarray()
else:
    X = np.array(X)
genes = adata.var_names

lfc = pd.DataFrame(index=pert_cols, columns=genes, dtype=float)
ctrl_mask = group_labels == 'Control'
for pert in pert_cols:
    pert_mask = group_labels == pert
    if pert_mask.sum() < 10:
        continue
    mu_p = X[pert_mask, :].mean(axis=0)
    mu_c = X[ctrl_mask, :].mean(axis=0)
    lfc.loc[pert] = mu_p - mu_c

all_groups = pert_cols + ['Control']
pvals_omni = pd.Series(index=genes, dtype=float)

for i, gene in enumerate(genes):
    expr_lists = [X[group_labels == grp, i] for grp in all_groups]
    if any(arr.size == 0 for arr in expr_lists):
        pvals_omni[gene] = np.nan
    else:
        _, p = kruskal(*expr_lists, nan_policy='omit')
        pvals_omni[gene] = p

pvals_omni.dropna(inplace=True)

alpha     = 0.05
p_flat    = pvals_omni.values

_, bonf_flat, _, _ = multipletests(p_flat, alpha=alpha, method='bonferroni')
bonf_omni = pd.Series(bonf_flat, index=pvals_omni.index)

_, fdr_flat, _, _ = multipletests(p_flat, alpha=alpha, method='fdr_bh')
fdr_omni = pd.Series(fdr_flat, index=pvals_omni.index)

sig_raw  = pvals_omni.index[pvals_omni  <= alpha]
sig_bonf = bonf_omni.index[bonf_omni    <= alpha]
sig_fdr  = fdr_omni.index[fdr_omni      <= alpha]

print(f"Genes significant (raw  p ≤ {alpha}): {len(sig_raw)}")
print(f"Genes significant (Bonf  p ≤ {alpha}): {len(sig_bonf)}")
print(f"Genes significant (FDR   p ≤ {alpha}): {len(sig_fdr)}")

sig_genes_fdr = list(sig_fdr)
sig_genes_bonf = list(sig_bonf)
sig_genes_raw = list(sig_raw)

with open('sig_DEGs_FDR.txt', 'w') as f:
    for gene in sig_genes_fdr:
        f.write(f"{gene}\n")

with open('sig_DEGs_bonferroni.txt', 'w') as f:
    for gene in sig_genes_bonf:
        f.write(f"{gene}\n")

df_sig = pd.DataFrame({
    'p_raw':  pvals_omni.loc[sig_genes_fdr],
    'p_bonf': bonf_omni.loc[sig_genes_fdr],
    'p_fdr':  fdr_omni.loc[sig_genes_fdr],
})
df_sig.to_csv('sig_DEGs_FDR_with_pvalues.csv')

sig_genes = list(sig_fdr)
plot_df   = lfc[sig_genes].dropna(axis=1, how='all')

plot_df_z = (
    plot_df
    .subtract(plot_df.mean(axis=0), axis=1)
    .div(plot_df.std(axis=0),   axis=1)
)

sns.set(context='notebook', style='white')
cg = sns.clustermap(
    plot_df_z,
    cmap='vlag',
    center=0,
    row_cluster=True,
    col_cluster=True,
    figsize=(36, 10),
    xticklabels=1,
    yticklabels=1,
    dendrogram_ratio=(.1, .2),
    cbar_pos=(0.02, .8, .05, .18)
)

cg.ax_heatmap.set_xlabel('Genes')
cg.ax_heatmap.set_ylabel('Perturbation target')
plt.suptitle('DEG heatmap (FDR-significant genes, omnibus KW test)', y=1.02)
plt.show()
