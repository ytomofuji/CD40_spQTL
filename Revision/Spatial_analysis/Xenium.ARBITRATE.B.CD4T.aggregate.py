#!/usr/bin/env python3
## ---------------------------
##
## Script name: Xenium.ARBITRATE.B.CD4T.aggregate.py
##
## Purpose of script: Define B/CD4T aggregates based on local cell density
##
##
## Date Created: 2026-11-13
##
## ---------------------------
##
## Notes:
##   1. for each cell, count cells within 50 um (cell itself included) and calculate fraction of B + CD4T cells
##   2. candidate cells: >= 10 cells within 50 um and B + CD4T fraction >= 0.3
##   3. DBSCAN on candidate cells (eps = 40 um, min_samples = 10)
##   4. aggregate: DBSCAN cluster with >= 30 candidate cells and >= 50 B or CD4T cells
##   all other cells: background
##   samples with > 5,000 cells and >= 100 B + CD4T cells were processed
##   usage: python Xenium.ARBITRATE.B.CD4T.aggregate.py [work_dir]
##
## ---------------------------

#cd /path/to/work_dir
#sbatch: --mem=64G --time=03:00:00


import sys
import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.dataset as ds
import pyarrow.compute as pc
from scipy.spatial import cKDTree
from sklearn.cluster import DBSCAN

work_dir = "/path/to/work_dir"
if len(sys.argv) > 1:
	work_dir = sys.argv[1]

name = "ARBITRATE.pre"
celltype_file = f"{work_dir}/{name}.celltype.tsv.gz"
polygon_file = "/path/to/allcells_qc_harmumapclust_lineage[metadata-tiles_coarsetypes_finetypes].parquet"

lymphoid_types = ["B", "CD4T"]
radius = 50
min_neighbors = 10
frac_thr = 0.3
eps = 40
min_samples = 10
min_cells = 30
min_lymphoid_cells = 50
min_cells_per_sample = 5000
min_lymphoid_per_sample = 100


def find_aggregates(df, frac, lym_mask):
	aggregate_id = np.full(len(df), -1, dtype=int)
	cand = (df["n_total_within_50um"].to_numpy() >= min_neighbors) & (np.nan_to_num(frac, nan=-1) >= frac_thr)
	rows = []
	if cand.sum() < min_cells:
		return pd.DataFrame(rows), aggregate_id, cand
	ci = np.where(cand)[0]
	labels = DBSCAN(eps=eps, min_samples=min_samples).fit_predict(df[["x", "y"]].to_numpy(np.float32)[ci])
	valid = []
	for lbl in np.unique(labels):
		if lbl < 0:
			continue
		m = labels == lbl
		if m.sum() < min_cells:
			continue
		n_lym = int(lym_mask[ci[m]].sum())
		if n_lym < min_lymphoid_cells:
			continue
		valid.append((lbl, int(m.sum()), n_lym))
	#aggregate id: ordered by size
	valid.sort(key=lambda v: -v[1])
	for new_id, (lbl, size, n_lym) in enumerate(valid):
		idx = ci[labels == lbl]
		g = df.iloc[idx]
		rows.append({"aggregate_id": new_id, "x": float(g["x"].mean()), "y": float(g["y"].mean()), "n_cells": size,
			"n_B_CD4T": n_lym, "n_B": int((g["celltype"] == "B").sum()), "n_CD4T": int((g["celltype"] == "CD4T").sum())})
		aggregate_id[idx] = new_id
	return pd.DataFrame(rows), aggregate_id, cand


#cell types and coordinates
meta = pd.read_csv(celltype_file, sep="\t", usecols=["cell_id", "sid", "celltype"])
geo = ds.dataset(polygon_file).to_table(columns=["cell_id", "x", "y"],
	filter=pc.is_in(pc.field("cell_id"), value_set=pa.array(meta["cell_id"].tolist()))).to_pandas()
cells = meta.merge(geo, on="cell_id", how="inner")
assert len(cells) == len(meta)

out_list, agg_list = [], []
for sid, sub in cells.groupby("sid", sort=True):
	df = sub.reset_index(drop=True).copy()
	coords = df[["x", "y"]].to_numpy(np.float32)
	lym_mask = df["celltype"].isin(lymphoid_types).to_numpy()

	#number of all cells and B + CD4T cells within 50 um
	df["n_total_within_50um"] = cKDTree(coords).query_ball_point(coords, r=radius, return_length=True)
	df["n_B_CD4T_within_50um"] = 0
	for t in lymphoid_types:
		m = (df["celltype"] == t).to_numpy()
		if m.any():
			df["n_B_CD4T_within_50um"] += cKDTree(coords[m]).query_ball_point(coords, r=radius, return_length=True)
	frac = df["n_B_CD4T_within_50um"].to_numpy() / np.clip(df["n_total_within_50um"].to_numpy(), 1, None)
	frac = np.where(df["n_total_within_50um"].to_numpy() < min_neighbors, np.nan, frac)
	df["frac_B_CD4T_within_50um"] = frac

	n_agg = 0
	if (len(df) > min_cells_per_sample) and (lym_mask.sum() >= min_lymphoid_per_sample):
		aggs, aggregate_id, cand = find_aggregates(df, frac, lym_mask)
		df["candidate"] = cand
		df["aggregate_id"] = aggregate_id
		df["region"] = np.where(aggregate_id >= 0, "aggregate", "background")
		n_agg = len(aggs)
		if len(aggs):
			aggs.insert(0, "sid", sid)
			agg_list.append(aggs)
	else:
		df["candidate"] = np.nan
		df["aggregate_id"] = np.nan
		df["region"] = np.nan
	out_list.append(df)
	print(f"{sid}: {len(df)} cells, {n_agg} aggregates", flush=True)

out = pd.concat(out_list, ignore_index=True)
out[["x", "y"]] = out[["x", "y"]].round(2)
out.drop(columns=["celltype"]).to_csv(f"{work_dir}/{name}.B.CD4T.aggregate.percell.tsv.gz", sep="\t", index=False, na_rep="NA")
pd.concat(agg_list, ignore_index=True).to_csv(f"{work_dir}/{name}.B.CD4T.aggregate.list.tsv", sep="\t", index=False)

summary = out.groupby("sid").apply(lambda g: pd.Series({"n_cells": len(g),
	"n_aggregates": g["aggregate_id"].dropna().pipe(lambda s: s[s >= 0].nunique()),
	"n_cells_in_aggregates": int((g["region"] == "aggregate").sum()),
	"n_B_in_aggregates": int(((g["region"] == "aggregate") & (g["celltype"] == "B")).sum()),
	"n_B_background": int(((g["region"] == "background") & (g["celltype"] == "B")).sum())}))
print(summary.to_string())
summary.to_csv(f"{work_dir}/{name}.B.CD4T.aggregate.per.sample.tsv", sep="\t")
