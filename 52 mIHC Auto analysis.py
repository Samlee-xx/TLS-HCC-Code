#%%
from sklearn.metrics import roc_auc_score, accuracy_score, recall_score, confusion_matrix
import pandas as pd
import itertools
from sklearn.ensemble import RandomForestClassifier
from sklearn.preprocessing import LabelEncoder
import os
import napari
import tifffile
import numpy as np
import imageio
from qtpy.QtWidgets import QPushButton
import warnings
from alpineer import io_utils
from skimage import io
from ark.segmentation import marker_quantification, segmentation_utils
from ark.utils import deepcell_service_utils, example_dataset, plot_utils
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
from scipy.ndimage import gaussian_filter1d
from skimage.filters import threshold_multiotsu
import numpy as np  #1.22.4
import matplotlib.pyplot as plt  #3.5.1
import seaborn as sns  #0.11.2
from shapely.geometry import MultiPoint, Point, Polygon  #1.8.2
from scipy.spatial import Voronoi #1.7.3
import shapely.geometry #1.8.2
import pandas as pd  #1.4.2
from sklearn.preprocessing import LabelEncoder  #1.0.2
import os
import pandas as pd
from alpineer import io_utils
from ark.utils.plot_utils import cohort_cluster_plot, color_segmentation_by_stat
import natsort as ns
import ark.settings as settings
import colorcet as cc
import colorcet
from matplotlib import colormaps
from sklearn.cluster import DBSCAN
from sklearn.neighbors import NearestNeighbors
from scipy.stats import rankdata
import os, numpy as np, matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter1d
from scipy.signal import find_peaks
from skimage.filters import threshold_multiotsu
from scipy.stats import rankdata

base_dir    = r"F:\Analysis\HEB\01 QUPATH OUT"
result_path = r"F:\Analysis\HEB\02 Segmentation"
cell_table_csv = os.path.join(result_path, r"cell_table\cell_table_arcsinh_transformed.csv")

# change the marker names
MARKERS = {
    "CD3":   "Cy7",
    "CD20":  "SpGold",
    "CD45":  "SpGreen",
    "CD66b": "Cy5",
}
dist_dir = os.path.join(result_path, "threshold_plots_manual")
os.makedirs(dist_dir, exist_ok=True)

CELL_SIZE_MIN = 30
CELL_SIZE_MAX = 500


df = pd.read_csv(cell_table_csv)
df = df[(df['cell_size'] > CELL_SIZE_MIN) & (df['cell_size'] < CELL_SIZE_MAX)].copy()

def parse_sample_region(fov_name):
    parts = str(fov_name).split('_')
    if len(parts) >= 2:
        sample = parts[0]
        region = parts[1]
    else:
        sample = 'unknown'
        region = str(fov_name)
    return sample, region

df[["sample_id", "region_id"]] = df["fov"].apply(lambda x: pd.Series(parse_sample_region(x)))
df["is_tls"] = df["region_id"].str.startswith("FOV")
df["is_bg"]  = df["region_id"].str.startswith("B")

tls_df = df[df["is_tls"]].copy()
bg_df  = df[df["is_bg"]].copy()

def plot_fov_vs_bg(values_fov, values_bg, marker_name, fov_id, save_path,
                   bins=120, smooth_sigma=2.0):
    vals_fov = np.asarray(values_fov)
    vals_fov = vals_fov[np.isfinite(vals_fov)]
    vals_bg  = np.asarray(values_bg)
    vals_bg  = vals_bg[np.isfinite(vals_bg)]

    merged = np.concatenate([vals_fov, vals_bg])
    if merged.size == 0:
        return

    lo, hi = np.percentile(merged, [0.5, 99.5])
    vals_fov = np.clip(vals_fov, lo, hi)
    vals_bg  = np.clip(vals_bg, lo, hi)

    hist_f, edges = np.histogram(vals_fov, bins=bins, density=True)
    xs = 0.5*(edges[:-1] + edges[1:])
    kde_f = gaussian_filter1d(hist_f.astype(float), sigma=smooth_sigma)

    hist_b, _ = np.histogram(vals_bg, bins=edges, density=True)
    kde_b = gaussian_filter1d(hist_b.astype(float), sigma=smooth_sigma)

    plt.figure(figsize=(7.2, 4.2))
    plt.hist(vals_bg, bins=edges, density=True, alpha=0.40, label="Background", edgecolor='none', color="#4C78A8")
    plt.plot(xs, kde_b, lw=1.5, color="#4C78A8", label="BG KDE")
    plt.hist(vals_fov, bins=edges, density=True, alpha=0.35, label="FOV", edgecolor='none', color="#F58518")
    plt.plot(xs, kde_f, lw=1.5, color="#F58518", label="FOV KDE")

    plt.title(f"{marker_name} | {fov_id}")
    plt.xlabel("Expression (arcsinh)")
    plt.ylabel("Density")
    plt.legend()
    plt.tight_layout()
    os.makedirs(os.path.dirname(save_path), exist_ok=True)
    plt.savefig(save_path, dpi=300)
    plt.close()

rows = []
for fov_id, fov_df in tls_df.groupby("fov"):
    sample_id = fov_df["sample_id"].iloc[0]
    sample_bg_df = bg_df[bg_df["sample_id"] == sample_id]

    for mk, col in MARKERS.items():
        if col not in fov_df.columns:
            continue
        bg_vals = sample_bg_df[col].values if col in sample_bg_df.columns else np.array([])
        out_png = os.path.join(dist_dir, f"{fov_id}_{mk}.png")

        plot_fov_vs_bg(fov_df[col].values, bg_vals, mk, fov_id, out_png)

        rows.append({
            "fov": fov_id,
            "sample_id": sample_id,
            "marker": mk,
            "column": col,
            "threshold": np.nan
        })

template = pd.DataFrame(rows, columns=["fov", "sample_id", "marker", "column", "threshold"])
template["threshold"] = pd.to_numeric(template["threshold"], errors="coerce")

template_long = os.path.join(result_path, "manual_thresholds_template_long.csv")
template.to_csv(template_long, index=False)

template_wide = template.pivot(index="fov", columns="marker", values="threshold")
template_wide = template_wide.reindex(columns=list(MARKERS.keys()))
template_wide_path = os.path.join(result_path, "manual_thresholds_template_wide.csv")
template_wide.to_csv(template_wide_path)

# %% 
import os
import numpy as np
import pandas as pd
from scipy.stats import rankdata
from sklearn.cluster import DBSCAN

result_path = r"F:\Analysis\HEB\02 Segmentation"
cell_table_csv = os.path.join(result_path, r"cell_table\cell_table_arcsinh_transformed.csv")
manual_thresholds_wide_csv = os.path.join(result_path, "manual_thresholds_template_wide.csv")

MARKERS = {
    "CD3":   "Cy7",
    "CD20":  "SpGold",
    "CD45":  "SpGreen",
    "CD66b": "Cy5",
}

seg_dir    = os.path.join(result_path, "deepcell_output")
plot_dir   = os.path.join(result_path, "seg_plot_celltype")
dbscan_dir = os.path.join(result_path, "seg_plot_dbscan")
for d in [plot_dir, dbscan_dir]:
    os.makedirs(d, exist_ok=True)

CELL_SIZE_MIN = 30
CELL_SIZE_MAX = 500
DBSCAN_EPS    = 90
DBSCAN_MIN_SAMPLES = 50
MULTI_MARGIN = 0.00

def load_manual_thresholds_wide(path):
    wide = pd.read_csv(path)
    wide["fov"] = wide["fov"].astype(str)
    for mk in MARKERS.keys():
        if mk not in wide.columns:
            wide[mk] = np.nan
    rows = []
    for _, r in wide.iterrows():
        for mk in MARKERS:
            rows.append({"fov": r["fov"], "marker": mk, "threshold": r[mk]})
    long_df = pd.DataFrame(rows).dropna(subset=["threshold"])
    thr_map = {(row["fov"], row["marker"]): float(row["threshold"])
               for _, row in long_df.iterrows()}
    return wide, long_df, thr_map
wide_thr, long_thr, thr_map = load_manual_thresholds_wide(manual_thresholds_wide_csv)

cell_seg = pd.read_csv(cell_table_csv)
cell_seg = cell_seg[
    (cell_seg['cell_size'] > CELL_SIZE_MIN) &
    (cell_seg['cell_size'] < CELL_SIZE_MAX)
].copy()
cell_seg["fov"] = cell_seg["fov"].astype(str)

def parse_sample_region(fov_name):
    parts = str(fov_name).split('_')
    if len(parts) >= 2:
        return parts[0], parts[1]
    else:
        return 'unknown', fov_name

cell_seg[["sample_id", "region_id"]] = cell_seg["fov"].apply(
    lambda x: pd.Series(parse_sample_region(x))
)
cell_seg["is_tls"] = cell_seg["region_id"].str.startswith("FOV")
allowed_fovs = set(wide_thr["fov"].astype(str))
df = cell_seg[cell_seg["is_tls"] & cell_seg["fov"].isin(allowed_fovs)].copy()
custom_cmap = pd.DataFrame({
    'cell_meta_cluster': ['CD20+', 'CD3+', 'CD66B+', 'Non-immune', 'Others'],
    'color': ['#1f77b4', '#ff7f0e', '#d62728', '#aec7e8', 'gray']
})
def compute_percentile(series, mask):
    vals = series[mask]
    if len(vals) == 0:
        return pd.Series([], dtype=float)
    pct = rankdata(vals, method='average') / len(vals)
    return pd.Series(pct, index=vals.index)


def _safe_concat(dfs, **kw):
    dfs = [d for d in dfs if hasattr(d, "empty") and not d.empty]
    return pd.concat(dfs, **kw) if len(dfs) > 0 else pd.DataFrame()

all_thr, all_counts, all_ctprop, all_rows, all_cluster_counts = [], [], [], [], []

for fov_id, df_fov in df.groupby("fov"):
    thresholds = {mk: thr_map[(fov_id, mk)] for mk in MARKERS}
    thr_records = []
    for mk, col in MARKERS.items():
        pos_rate = float((df_fov[col] > thresholds[mk]).mean())
        thr_records.append({
            "fov": fov_id,
            "marker": mk,
            "threshold": thresholds[mk],
            "method": "manual",
            "pos_rate": pos_rate
        })
    cd45_pos = df_fov[MARKERS["CD45"]] > thresholds["CD45"]
    df_immune = df_fov[cd45_pos].copy()
    cd3_pos   = df_immune[MARKERS["CD3"]]   > thresholds["CD3"]
    cd20_pos  = df_immune[MARKERS["CD20"]]  > thresholds["CD20"]
    cd66b_pos = df_immune[MARKERS["CD66b"]] > thresholds["CD66b"]
    df_immune["celltype"] = "Unassigned"
    only_cd3   =  cd3_pos & ~cd20_pos & ~cd66b_pos
    only_cd20  = ~cd3_pos &  cd20_pos & ~cd66b_pos
    only_cd66b = ~cd3_pos & ~cd20_pos &  cd66b_pos
    all_neg    = ~cd3_pos & ~cd20_pos & ~cd66b_pos
    df_immune.loc[only_cd3,   "celltype"] = "CD3+"
    df_immune.loc[only_cd20,  "celltype"] = "CD20+"
    df_immune.loc[only_cd66b, "celltype"] = "CD66B+"
    df_immune.loc[all_neg,    "celltype"] = "Others"
    df_immune["CD3_percentile"]   = compute_percentile(df_immune[MARKERS["CD3"]], cd3_pos)
    df_immune["CD20_percentile"]  = compute_percentile(df_immune[MARKERS["CD20"]], cd20_pos)
    df_immune["CD66B_percentile"] = compute_percentile(df_immune[MARKERS["CD66b"]], cd66b_pos)

    multi_pos = (cd3_pos.astype(int) + cd20_pos.astype(int) + cd66b_pos.astype(int)) > 1
    multi_idx = df_immune[(df_immune["celltype"] == "Unassigned") & multi_pos].index

    for idx in multi_idx:
        percs = {
            "CD3+": float(df_immune.at[idx, "CD3_percentile"]) if pd.notna(df_immune.at[idx, "CD3_percentile"]) else 0,
            "CD20+":float(df_immune.at[idx, "CD20_percentile"])if pd.notna(df_immune.at[idx, "CD20_percentile"]) else 0,
            "CD66B+":float(df_immune.at[idx,"CD66B_percentile"])if pd.notna(df_immune.at[idx,"CD66B_percentile"])else 0,
        }
        top2 = sorted(percs.items(), key=lambda kv: kv[1], reverse=True)[:2]
        if top2[0][1] - top2[1][1] >= MULTI_MARGIN:
            df_immune.at[idx, "celltype"] = top2[0][0]
        else:
            df_immune.at[idx, "celltype"] = "Unassigned"
    df_out = df_fov.copy()
    df_out["celltype"] = "Non-immune"
    df_out.loc[df_immune.index, "celltype"] = df_immune["celltype"]
    df_out["cluster_label"] = -2
    immune_idx = df_immune.index

    if len(immune_idx) >= DBSCAN_MIN_SAMPLES:
        X = df_out.loc[immune_idx, ["centroid-0","centroid-1"]].values
        db = DBSCAN(eps=DBSCAN_EPS, min_samples=DBSCAN_MIN_SAMPLES).fit(X)
        df_out.loc[immune_idx, "cluster_label"] = db.labels_
    else:
        print(f"[INFO] fov {fov_id}: Skipp DBSCAN")

    valid_mask = df_out["cluster_label"] >= 0
    cluster_counts = df_out.loc[valid_mask, ["fov","cluster_label"]] \
                           .groupby(["fov","cluster_label"]) \
                           .size().reset_index(name="cell_count")
    all_cluster_counts.append(cluster_counts)

    cell_cluster = df_out[["fov", "label", "celltype"]].rename(columns={"celltype": "cell_meta_cluster"})
    cohort_cluster_plot(
        fovs=[fov_id], seg_dir=seg_dir, save_dir=plot_dir,
        cell_data=cell_cluster, erode=True,
        fov_col=settings.FOV_ID, label_col=settings.CELL_LABEL,
        cluster_col=settings.CELL_TYPE, 
        seg_suffix="_whole_cell.tiff", cmap=custom_cmap,
        display_fig=False, fig_file_type="png", dpi=300
    )

    df_dbscan_plot = df_out[df_out["cluster_label"] >= -1].copy()
    cell_cluster2 = df_dbscan_plot[["fov", "label", "cluster_label"]].rename(columns={"cluster_label": "cell_meta_cluster"})
    cohort_cluster_plot(
        fovs=[fov_id], seg_dir=seg_dir, save_dir=dbscan_dir,
        cell_data=cell_cluster2, erode=True,
        fov_col=settings.FOV_ID, label_col=settings.CELL_LABEL,
        cluster_col="cell_meta_cluster",  
        seg_suffix="_whole_cell.tiff",
        display_fig=False, fig_file_type="png", dpi=300
    )

    df_valid = df_out[df_out["cluster_label"] > -1]
    if df_valid.empty:
        ct_prop = pd.DataFrame()
    else:
        ct = pd.crosstab(df_valid["cluster_label"], df_valid["celltype"])
        ct_prop = ct.div(ct.sum(axis=1), axis=0) * 100

    all_thr.append(pd.DataFrame(thr_records))
    all_counts.append(df_out["celltype"].value_counts().to_frame(name="count").assign(fov=fov_id))
    if not ct_prop.empty:
        all_ctprop.append(ct_prop.assign(fov=fov_id))
    all_rows.append(df_out)

thr_summary  = _safe_concat(all_thr,  ignore_index=True)
counts_sum   = _safe_concat(all_counts)
ctprop_sum   = _safe_concat(all_ctprop)
df_labeled   = _safe_concat(all_rows, ignore_index=True)

cluster_counts_sum = _safe_concat(all_cluster_counts, ignore_index=True)

thr_summary.to_csv(os.path.join(result_path, "thresholds_per_fov.csv"), index=False)
counts_sum.to_csv(os.path.join(result_path, "celltype_counts_per_fov.csv"))
df_labeled.to_csv(os.path.join(result_path, "cells_with_types_and_clusters.csv"), index=False)

if not ctprop_sum.empty:
    ctprop_sum.to_csv(os.path.join(result_path, "cluster_celltype_prop_per_fov.csv"), index=False)

if not cluster_counts_sum.empty:
    cluster_counts_sum.to_csv(os.path.join(result_path, "cluster_cell_counts_per_fov.csv"), index=False)
    cluster_counts_wide = cluster_counts_sum.pivot_table(
        index="fov", columns="cluster_label",
        values="cell_count", aggfunc="sum", fill_value=0
    )
    cluster_counts_wide.to_csv(os.path.join(result_path, "cluster_cell_counts_per_fov_wide.csv"))

# %%
