#!/usr/bin/env python3
## ---------------------------
##
## Script name: Xenium.ARBITRATE.B.CD4T.distance.py
##
## Purpose of script: Calculate distance from each B cell to the nearest CD4T cell (cell polygons)
##
##
## Date Created: 2026-11-13
##
## ---------------------------
##
## Notes:
##   distance: polygon-to-polygon distance (shapely) to the nearest CD4T cell in the same sample
##   count: number of CD4T cells whose polygon is within r um of the B cell centroid
##   usage: python Xenium.ARBITRATE.B.CD4T.distance.py [work_dir]
##
## ---------------------------

#cd /path/to/work_dir
#sbatch: --mem=64G --time=06:00:00


import sys
import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.dataset as ds
import pyarrow.compute as pc
import shapely
from shapely.strtree import STRtree

work_dir = "/path/to/work_dir"
if len(sys.argv) > 1:
	work_dir = sys.argv[1]

name = "ARBITRATE.pre"
celltype_file = f"{work_dir}/{name}.celltype.tsv.gz"
score_file = f"{work_dir}/{name}.B.cell.CRAFT.score.tsv.gz"
polygon_file = "/path/to/allcells_qc_harmumapclust_lineage[metadata-tiles_coarsetypes_finetypes].parquet"
radii = [10, 50]


#B cells (after sample QC) and CD4T cells
meta = pd.read_csv(celltype_file, sep="\t", usecols=["cell_id", "sid", "celltype"])
b_cells = pd.read_csv(score_file, sep="\t", usecols=["cell_id"])["cell_id"]
meta = meta[(meta["cell_id"].isin(b_cells)) | (meta["celltype"] == "CD4T")].copy()
meta["isB"] = meta["cell_id"].isin(b_cells)

#cell centroids and polygons
tab = ds.dataset(polygon_file).to_table(columns=["cell_id", "x", "y", "polygon_wkt"],
	filter=pc.is_in(pc.field("cell_id"), value_set=pa.array(meta["cell_id"].tolist())))
df = meta.merge(tab.to_pandas(), on="cell_id", how="inner")
assert len(df) == len(meta)

res = []
for sid, sub in df.groupby("sid"):
	b = sub[sub["isB"]].reset_index(drop=True)
	t = sub[(~sub["isB"]) & (sub["celltype"] == "CD4T")].reset_index(drop=True)
	if len(b) == 0:
		continue
	out = b[["cell_id", "sid"]].copy()
	b_poly = shapely.from_wkt(b["polygon_wkt"].to_numpy())
	b_pts = shapely.points(b[["x", "y"]].to_numpy())
	t_poly = shapely.from_wkt(t["polygon_wkt"].to_numpy())
	tree = STRtree(t_poly)

	#distance to the nearest CD4T cell
	(bi, ti), d = tree.query_nearest(b_poly, return_distance=True, all_matches=False)
	dist = np.full(len(b), np.nan)
	dist[bi] = d
	out["dist_CD4T_um"] = dist

	#number of CD4T cells within r um of the B cell centroid
	pi, pt = tree.query(b_pts, predicate="dwithin", distance=max(radii))
	dd = shapely.distance(b_pts[pi], t_poly[pt])
	for r in radii:
		out[f"n_CD4T_within_{r}um"] = np.bincount(pi[dd <= r], minlength=len(b))
	res.append(out)
	print(f"{sid}: {len(b)} B cells, {len(t)} CD4T cells", flush=True)

pd.concat(res, ignore_index=True).to_csv(f"{work_dir}/{name}.B.CD4T.distance.tsv.gz", sep="\t", index=False)
