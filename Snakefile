
import os

#####################################
### Load Configurations
#####################################

configfile: "config.yaml"

name_of_project       = config["name_of_project"]
metadata_file         = config["metadata_file"]
HALO_data_file        = config["HALO_data_file"]
output_directory      = config["output_directory"]
pca_dims              = int(config["pca_dims"])
clustering_resolution = float(config["clustering_resolution"])

# Absolute path to the pipeline folder (where this Snakefile lives).
# Passed into run: blocks so app.py is always found, whether running inside
# Docker (/pipeline) or locally on macOS.
PIPELINE_DIR = os.path.dirname(os.path.abspath(workflow.snakefile))

######################################
### Prevent thread over-subscription in numpy/scipy backends
######################################
os.environ["OMP_NUM_THREADS"]        = "1"
os.environ["MKL_NUM_THREADS"]        = "1"
os.environ["OPENBLAS_NUM_THREADS"]   = "1"
os.environ["BLAS_NUM_THREADS"]       = "1"
os.environ["VECLIB_MAXIMUM_THREADS"] = "1"
os.environ["NUMEXPR_NUM_THREADS"]    = "1"

######################################
### Auto-generate metadata_markers.csv if absent.
### Runs at Snakefile parse time — works whether the pipeline is launched
### via RUN_COMET_PIPELINE.BAT on Windows or `snakemake` directly on macOS.
######################################
if not os.path.exists(metadata_file):
	try:
		import sys
		if PIPELINE_DIR not in sys.path:
			sys.path.insert(0, PIPELINE_DIR)
		from generate_metadata import generate_from_csv
		print(
			"\n[COMET] '{}' not found — auto-generating from '{}' ...".format(
				metadata_file, HALO_data_file
			)
		)
		generate_from_csv(HALO_data_file, metadata_file)
		print(
			"[COMET] Done. Please open '{}' and confirm the "
			"type/localization columns before proceeding.\n".format(metadata_file)
		)
	except ImportError:
		raise FileNotFoundError(
			"'{}' not found and generate_metadata.py is not available "
			"in the pipeline directory. Either create metadata_markers.csv "
			"manually or add generate_metadata.py to the pipeline folder.".format(
				metadata_file
			)
		)


######################################

rule all:
	input:
		f"{output_directory}/{name_of_project}_annotated_seurat.rds",
		f"{output_directory}/images_final/umap_cell_types.png"


rule process_halo:
	input:
		csv  = HALO_data_file,
		meta = metadata_file
	output:
		h5ad = f"{output_directory}/{name_of_project}.h5ad"
	run:
		import os
		import re
		import unicodedata

		import numpy as np
		import pandas as pd
		import scanpy as sc
		import anndata as ad
		import spatialdata as sd
		import geopandas as gpd
		import matplotlib
		matplotlib.use("Agg")
		import matplotlib.pyplot as plt
		import harpy as hp
		from shapely.geometry import Point
		from spatialdata.models import ShapesModel, TableModel
		from spatialdata.transformations import Identity

		sc.settings.autoshow = False
		os.makedirs(output_directory, exist_ok=True)
		os.makedirs("{}/images".format(output_directory), exist_ok=True)

		# ── 1. Read the CSV robustly ──────────────────────────────────────────
		# HALO exports come in two layouts depending on locale / version:
		#
		#  Standard: each row is a normal CSV row, possibly with quoted fields.
		#    Image Location,"Analysis Region","Object Id","XMin",...
		#    \\path\file.tiff,"field of view","0","456,789",...
		#
		#  Wrapped:  each entire row is itself enclosed in one outer pair of
		#    quotes, with all internal quotes doubled as "".  Example:
		#    "Image Location,""Analysis Region"",""Object Id"",""XMin"",..."
		#    "\\path\file.tiff,""field of view"",""0"",""456,789"",..."
		#
		# The wrapped format is produced by certain HALO versions on
		# European-locale Windows.  Pandas cannot read it directly because the
		# outer quotes make every row look like a single field.
		# We detect and unescape the wrapped format before handing off to pandas.
		import io as _io

		def _try_standard(path, enc):
			"""Standard CSV: try comma / semicolon / tab as separator."""
			for _s in (",", ";", "\t"):
				try:
					_c = pd.read_csv(path, sep=_s, decimal=",", encoding=enc, index_col=False)
					if _c.shape[1] >= 5 and not any("¬" in str(col) for col in _c.columns):
						return _c
				except Exception:
					pass
			return None

		def _try_wrapped(path, enc):
			"""
			Wrapped HALO format: strip the outer per-row quotes and unescape
			doubled internal quotes, then parse the resulting standard CSV.

			Raw:   "Image Location,""Analysis Region"",""XMin"",..."
			After: Image Location,"Analysis Region","XMin",...
			"""
			try:
				with open(path, "r", encoding=enc, errors="replace", newline="") as _fh:
					raw_lines = _fh.read().splitlines()
				# Detect: first non-empty line starts with " and contains "" inside
				_first = next((l for l in raw_lines if l), "")
				if not (_first.startswith('"') and '""' in _first):
					return None
				# Unescape each line: strip one outer " from start and end, then "" -> "
				unescaped = []
				for _line in raw_lines:
					if _line.startswith('"') and _line.endswith('"'):
						_line = _line[1:-1].replace('""', '"')
					unescaped.append(_line)
				_content = "\n".join(unescaped)
				_c = pd.read_csv(_io.StringIO(_content), decimal=",", index_col=False)
				if _c.shape[1] >= 5 and not any("¬" in str(col) for col in _c.columns):
					return _c
			except Exception:
				pass
			return None

		df = None
		for _enc in ("utf-8-sig", "cp1252", "latin-1"):
			df = _try_standard(input.csv, _enc)
			if df is None:
				df = _try_wrapped(input.csv, _enc)
			if df is not None:
				break

		if df is None:
			raise ValueError(
				"Could not read '{}' cleanly with any known encoding or format. "
				"Check that the file is a valid HALO or Horizon CSV export. "
				"The file should be either a standard CSV or the HALO "
				"whole-row-quoted format.".format(input.csv)
			)

		# ── 2. Normalise column names ─────────────────────────────────────────
		# NFKC collapses µ (U+00B5) → μ (U+03BC) and ² (U+00B2) → 2, giving
		# consistent column names regardless of which software wrote the file.
		# Replacing ? and α → a unifies αSMA/?SMA with the "aSMA" convention
		# used in metadata_markers.csv.
		def _norm(c):
			return (
				unicodedata.normalize("NFKC", str(c))
				.replace("?", "a")
				.replace("\u03b1", "a")   # α
				.strip()
			)
		df.columns = [_norm(c) for c in df.columns]

		# ── 3. Load & validate metadata_markers.csv ───────────────────────────
		markers_df = pd.read_csv(input.meta, index_col=False, skiprows=[0, 1])
		markers_df.columns = [c.strip() for c in markers_df.columns]
		for _col in markers_df.select_dtypes(include="object").columns:
			markers_df[_col] = markers_df[_col].str.strip()

		bad_types = markers_df.type[~markers_df.type.isin(["Protein", "Transcript", "Other"])]
		if len(bad_types):
			raise ValueError(
				"Unknown value(s) in 'type' column of metadata_markers.csv: {}. "
				"Allowed values: Protein, Transcript, Other.".format(sorted(bad_types.unique()))
			)
		bad_locs = markers_df.localization[
			~markers_df.localization.isin(["Nucleus", "Cytoplasm", "NA"])
		]
		if len(bad_locs):
			raise ValueError(
				"Unknown value(s) in 'localization' column of metadata_markers.csv: {}. "
				"Allowed values: Nucleus, Cytoplasm, NA.".format(sorted(bad_locs.unique()))
			)

		# ── 4. Detect input format ────────────────────────────────────────────
		_halo_cols    = {"Object Id", "XMin", "XMax", "YMin", "YMax"}
		_horizon_cols = {"Annotation Group", "Annotation Index"}

		if _halo_cols.issubset(df.columns):
			input_format = "HALO"
			for _col in ["XMin", "XMax", "YMin", "YMax", "Cell Area (μm2)"]:
				if _col in df.columns:
					df[_col] = pd.to_numeric(
						df[_col].astype(str).str.replace(",", ".", regex=False),
						errors="coerce",
					)
			df["x"]               = (df["XMin"] + df["XMax"]) / 2
			df["y"]               = (df["YMin"] + df["YMax"]) / 2
			df["Object_Id_Clean"] = df["Object Id"].astype(int) + 1
			df["cell_area_um2"]   = df["Cell Area (μm2)"]

		elif _horizon_cols.issubset(df.columns) and any(c.startswith("Nuclei/") for c in df.columns):
			input_format = "HORIZON"
			df["x"]               = pd.to_numeric(df["Nuclei/X Position in μm"], errors="coerce")
			df["y"]               = pd.to_numeric(df["Nuclei/Y Position in μm"], errors="coerce")
			# Annotation Index resets per ROI — use global row numbers for a unique id
			df["Object_Id_Clean"] = np.arange(1, len(df) + 1)
			df["cell_area_um2"]   = pd.to_numeric(df["Nuclei/Area in μm2"], errors="coerce")

		else:
			raise ValueError(
				"'{}' does not look like a HALO or Horizon export. "
				"Expected 'Object Id'/'XMin'/'XMax'/'YMin'/'YMax' (HALO) "
				"or 'Annotation Group'/'Nuclei/...' (Horizon).".format(input.csv)
			)

		# Coerce shared spatial columns to float (handles any stragglers)
		for _col in ["x", "y", "cell_area_um2"]:
			df[_col] = pd.to_numeric(
				df[_col].astype(str).str.replace(",", ".", regex=False),
				errors="coerce",
			)

		pixel_size          = 0.5    # µm per pixel
		df["cell_area_px"]  = df["cell_area_um2"] / (pixel_size ** 2)
		df["radius_px"]     = np.sqrt(df["cell_area_px"] / np.pi)

		# ── 5. Match metadata rows to data columns ────────────────────────────
		def find_marker_column(columns, channel, col_type, localization, fmt):
			"""
			Return the export column that matches one metadata_markers.csv row,
			or None if no match is found.

			Handles:
			HALO plain:   'CD3 Cytoplasm Intensity'
			HALO cycled:  '1 | CD3_500x Cytoplasm Intensity'
			Horizon:      'Nuclei/Mean Intensity (CD3_500x - TRITC Protein Autofluo)'
			"""
			esc      = re.escape(channel)
			dilution = r"(?:_\d+(?:\s?\d+)?x)?"   # optional _500x / _20 000x

			if fmt == "HALO":
				if col_type in ("Protein", "Other"):
					suffix = "Cytoplasm Intensity" if localization == "Cytoplasm" else "Nucleus Intensity"
					pat = re.compile(
						r"^(?:\d+\s*\|\s*)?" + esc + dilution + r"\s+" + suffix + r"$",
						re.IGNORECASE,
					)
				else:  # Transcript
					pat = re.compile(
						r"^(?:\d+\s*\|\s*)?T\d+\s+" + esc + dilution + r"\s+Copies$",
						re.IGNORECASE,
					)
			else:  # HORIZON
				kind = "RNA" if col_type == "Transcript" else "Protein"
				pat = re.compile(
					r"Nuclei/Mean Intensity\s*\(" + esc + dilution + r"\s*-\s*\w+\s+" + kind + r"\s+Autofluo\)$",
					re.IGNORECASE,
				)

			# If a marker appears in multiple cycles, keep the last (highest cycle)
			matches = sorted(c for c in columns if pat.search(c))
			return matches[-1] if matches else None

		matched_rows    = []
		matched_columns = []
		for _, row in markers_df.iterrows():
			col = find_marker_column(
				df.columns, row["channel"], row["type"], row["localization"], input_format
			)
			if col is None:
				print(
					"WARNING: no column found in '{}' for marker '{}' "
					"(type={}, localization={}) — skipping.".format(
						input.csv, row["channel"], row["type"], row["localization"]
					)
				)
				continue
			matched_rows.append(row)
			matched_columns.append(col)

		if not matched_rows:
			raise ValueError(
				"None of the markers in '{}' matched any column in '{}' "
				"(detected format: {}). Check that channel names in "
				"metadata_markers.csv match the export.".format(
					input.meta, input.csv, input_format
				)
			)

		matched_meta = pd.DataFrame(matched_rows).reset_index(drop=True)
		clean_names  = matched_meta["channel"].tolist()

		# ── 6. Build AnnData ──────────────────────────────────────────────────
		# Ensure expression columns are numeric before casting to float32
		for _col in matched_columns:
			df[_col] = pd.to_numeric(
				df[_col].astype(str).str.replace(",", ".", regex=False),
				errors="coerce",
			)

		adata = ad.AnnData(
			X   = df[matched_columns].values.astype(np.float32),
			obs = df[["Object_Id_Clean"]].copy(),
			var = matched_meta.set_index(pd.Index(clean_names)),
		)
		adata.layers["raw"]        = adata.X.copy()
		adata.obs["region"]        = pd.Categorical(["halo_labels"] * len(adata))
		adata.obs["instance_id"]   = df["Object_Id_Clean"].astype(int).values
		adata.obs_names            = df["Object_Id_Clean"].astype(str).values
		adata.obs["cell_area_um2"] = df["cell_area_um2"].values
		adata.obsm["spatial"]      = df[["x", "y"]].values.astype(np.float32)
		adata.var["feature_type"]  = adata.var["type"]

		# ── 7. Build SpatialData ──────────────────────────────────────────────
		sdata    = sd.SpatialData()
		geometry = [Point(xy) for xy in zip(df["x"], df["y"])]
		gdf      = gpd.GeoDataFrame(df, geometry=geometry)
		gdf["geometry"] = gdf.buffer(df["radius_px"], cap_style=1)
		gdf.index       = df["Object_Id_Clean"].astype(str).values

		cells_shapes = ShapesModel.parse(gdf, transformations={"global": Identity()})
		sdata.shapes["halo_cells"] = cells_shapes

		target_shape = (44643, 44643)
		hp.im.rasterize(
			sdata        = sdata,
			shapes_layer = "halo_cells",
			output_layer = "halo_labels",
			out_shape    = target_shape,
			chunks       = 2048,
			overwrite    = True,
		)

		sdata.table = TableModel.parse(
			adata,
			region       = "halo_labels",
			region_key   = "region",
			instance_key = "instance_id",
		)

		# ── 8. Normalise & scale ──────────────────────────────────────────────
		print("Normalizing")
		is_protein = adata.var["feature_type"] == "Protein"
		is_rna     = adata.var["feature_type"] == "Transcript"

		# .copy() is essential — without it, 'normalized' and 'raw' share memory
		# and updating one silently overwrites the other.
		adata.layers["normalized"] = adata.layers["raw"].copy()
		adata.layers["normalized"][:, is_protein] = np.arcsinh(
			adata.layers["raw"][:, is_protein] / 50.0
		)
		adata.layers["normalized"][:, is_rna] = np.log1p(
			adata.layers["raw"][:, is_rna]
		)
		adata.X = adata.layers["normalized"].copy()

		if not isinstance(adata.X, np.ndarray):
			adata.X = adata.X.toarray()
		adata.X = np.ascontiguousarray(adata.X, dtype=np.float32)

		print("Scaling")
		sc.pp.scale(adata, max_value=10)

		# ── 9. PCA / UMAP / Leiden (protein markers only) ────────────────────
		clustering_markers = [
			m for m in sdata.table.var_names
			if sdata.table.var.loc[m, "feature_type"] == "Protein"
		]
		adata_subset = sdata.table[:, clustering_markers].copy()

		print("Performing PCA")
		n_comps = min(10, adata_subset.shape[1])
		sc.tl.pca(adata_subset, n_comps=n_comps)
		sc.pl.pca_variance_ratio(adata_subset, n_pcs=n_comps, show=False)
		plt.savefig("{}/images/pca_elbowPlot.png".format(output_directory), bbox_inches="tight")
		plt.close()

		print("Finding Neighbors")
		sc.pp.neighbors(adata_subset, n_pcs=pca_dims)
		print("Performing UMAP")
		sc.tl.umap(adata_subset)

		sdata.table.obsm["X_pca"]         = adata_subset.obsm["X_pca"]
		sdata.table.obsm["X_umap"]        = adata_subset.obsm["X_umap"]
		sdata.table.uns["neighbors"]       = adata_subset.uns["neighbors"]
		sdata.table.obsp["distances"]      = adata_subset.obsp["distances"]
		sdata.table.obsp["connectivities"] = adata_subset.obsp["connectivities"]

		# Store PCA loadings for all vars (zeros for non-protein features)
		full_pcs    = np.zeros((sdata.table.shape[1], adata_subset.varm["PCs"].shape[1]))
		var_indices = [sdata.table.var_names.get_loc(m) for m in clustering_markers]
		full_pcs[var_indices, :] = adata_subset.varm["PCs"]
		sdata.table.varm["PCs"] = full_pcs

		print("Performing Leiden Clustering")
		sc.tl.leiden(
			sdata.table,
			resolution = clustering_resolution,
			key_added  = "leiden_clusters",
		)

		# ── 10. Diagnostic plots ──────────────────────────────────────────────
		sc.pl.umap(sdata.table, color="cell_area_um2", show=False)
		plt.savefig("{}/images/umap_cell_areas.png".format(output_directory), bbox_inches="tight")
		plt.close()

		sc.pl.umap(sdata.table, color="leiden_clusters", show=False)
		plt.savefig(
			"{}/images/umap_leiden_clustering_res{}.png".format(output_directory, clustering_resolution),
			bbox_inches="tight",
		)
		plt.close()

		for _var in sdata.table.var_names:
			sc.pl.umap(sdata.table, color=_var, show=False)
			plt.savefig("{}/images/umap_{}.png".format(output_directory, _var), bbox_inches="tight")
			plt.close()

		sc.pl.dotplot(
			sdata.table, var_names=sdata.table.var_names,
			groupby="leiden_clusters", show=False
		)
		plt.savefig("{}/images/dotplot_allChannels.png".format(output_directory), bbox_inches="tight")
		plt.close()

		protein_markers = [
			m for m in sdata.table.var_names
			if sdata.table.var.loc[m, "feature_type"] == "Protein"
		]
		if protein_markers:
			sc.pl.dotplot(
				sdata.table, var_names=protein_markers,
				groupby="leiden_clusters", show=False
			)
			plt.savefig("{}/images/dotplot_proteinChannels.png".format(output_directory), bbox_inches="tight")
			plt.close()

		rna_markers = [
			m for m in sdata.table.var_names
			if sdata.table.var.loc[m, "feature_type"] == "Transcript"
		]
		if rna_markers:
			sc.pl.dotplot(
				sdata.table, var_names=rna_markers,
				groupby="leiden_clusters", show=False
			)
			plt.savefig("{}/images/dotplot_rnaChannels.png".format(output_directory), bbox_inches="tight")
			plt.close()

		print("Performing Neighborhood Enrichment Analysis")
		hp.tb.nhood_enrichment(
			sdata,
			labels_layer    = "halo_labels",
			table_layer     = "table",
			output_layer    = "table_score_genes_enrichment",
			celltype_column = "leiden_clusters",
		)
		hp.pl.nhood_enrichment(
			sdata,
			table_layer     = "table_score_genes_enrichment",
			celltype_column = "leiden_clusters",
			output          = "{}/images/nhood_enrichment.png".format(output_directory),
		)

		# ── 11. Save ──────────────────────────────────────────────────────────
		print("Saving Output")
		adata.obs["x"]         = df["x"].values.astype(np.float32)
		adata.obs["y"]         = df["y"].values.astype(np.float32)
		adata.obs["radius_px"] = df["radius_px"].values.astype(np.float32)

		sdata.table.copy().write_h5ad(output.h5ad)
		print("Finished Successfully!")


rule annotate_clusters:
	input:
		h5ad = f"{output_directory}/{name_of_project}.h5ad"
	output:
		annotated_h5ad = f"{output_directory}/{name_of_project}_annotated.h5ad"
	run:
		import os
		import subprocess
		import time

		print("\n" + "=" * 60)
		print("PIPELINE PAUSED: LAUNCHING INTERACTIVE ANNOTATION PORTAL")
		print("Open http://localhost:8501 in your browser.")
		print("Complete annotations and click 'Finalize' to resume.")
		print("=" * 60 + "\n")

		# Locate app.py relative to this Snakefile so the path works both
		# inside Docker (/pipeline/app.py) and locally on macOS.
		app_path = os.path.join(PIPELINE_DIR, "app.py")
		if not os.path.exists(app_path):
			raise FileNotFoundError(
				"app.py not found at '{}'. "
				"Ensure app.py is in the same folder as this Snakefile.".format(app_path)
			)

		cmd = [
			"streamlit", "run", app_path,
			"--server.address=0.0.0.0",
			"--server.port=8501",
			"--server.headless=true",
			"--",
			"--input",   input.h5ad,
			"--output",  output.annotated_h5ad,
			"--img_dir", "{}/images".format(output_directory),
		]

		process = subprocess.Popen(cmd)

		# Block Snakemake until the researcher saves annotations
		while not os.path.exists(output.annotated_h5ad):
			time.sleep(2)

		print("\nAnnotations saved. Resuming pipeline ...")
		process.terminate()
		process.wait()


rule post_annotation_viz:
	input:
		h5ad = f"{output_directory}/{name_of_project}_annotated.h5ad"
	output:
		umap_img  = f"{output_directory}/images_final/umap_cell_types.png",
		dot_img   = f"{output_directory}/images_final/dotplot_final_proteins.png",
		nhood_img = f"{output_directory}/images_final/nhood_enrichment_final.png"
	run:
		import os
		import matplotlib
		matplotlib.use("Agg")
		import matplotlib.pyplot as plt
		import scanpy as sc
		import harpy as hp
		import spatialdata as sd
		import geopandas as gpd
		from shapely.geometry import Point
		from spatialdata.models import ShapesModel, TableModel
		from spatialdata.transformations import Identity

		os.makedirs(os.path.dirname(output.umap_img), exist_ok=True)
		adata = sc.read_h5ad(input.h5ad)

		# UMAP coloured by annotated cell type
		sc.pl.umap(adata, color="cell_type", show=False)
		plt.savefig(output.umap_img, bbox_inches="tight")
		plt.close()

		# Dotplot of protein markers grouped by cell type
		protein_markers = [
			m for m in adata.var_names
			if adata.var.loc[m, "feature_type"] == "Protein"
		]
		if protein_markers:
			sc.pl.dotplot(
				adata, var_names=protein_markers,
				groupby="cell_type", show=False
			)
			plt.savefig(output.dot_img, bbox_inches="tight")
			plt.close()

		# Neighbourhood enrichment with biological cell-type names.
		# Re-construct the SpatialData object from coordinates saved in obs.
		sdata    = sd.SpatialData()
		geometry = [Point(xy) for xy in zip(adata.obs["x"], adata.obs["y"])]
		gdf      = gpd.GeoDataFrame(adata.obs, geometry=geometry)
		gdf["geometry"] = gdf.buffer(adata.obs["radius_px"], cap_style=1)
		gdf.index       = adata.obs["instance_id"].astype(str).values

		cells_shapes = ShapesModel.parse(gdf, transformations={"global": Identity()})
		sdata.shapes["halo_cells"] = cells_shapes

		target_shape = (44643, 44643)
		hp.im.rasterize(
			sdata        = sdata,
			shapes_layer = "halo_cells",
			output_layer = "halo_labels",
			out_shape    = target_shape,
			chunks       = 2048,
			overwrite    = True,
		)

		if "spatialdata_attrs" in adata.uns:
			del adata.uns["spatialdata_attrs"]

		sdata.table = TableModel.parse(
			adata,
			region       = "halo_labels",
			region_key   = "region",
			instance_key = "instance_id",
		)

		hp.tb.nhood_enrichment(
			sdata,
			labels_layer    = "halo_labels",
			table_layer     = "table",
			celltype_column = "cell_type",
			output_layer    = "table_annotated_enrichment",
		)
		hp.pl.nhood_enrichment(
			sdata,
			table_layer     = "table_annotated_enrichment",
			celltype_column = "cell_type",
			output          = output.nhood_img,
		)

		print("Final annotated visualizations generated successfully!")


rule convert_to_seurat:
	input:
		h5ad = f"{output_directory}/{name_of_project}_annotated.h5ad"
	output:
		rds = f"{output_directory}/{name_of_project}_annotated_seurat.rds"
	run:
		import os
		import subprocess
		import scipy.sparse
		import pandas as pd
		import scanpy as sc

		print("Converting AnnData to Seurat object ...")
		adata = sc.read_h5ad(input.h5ad)

		tmp_prefix = "{}/tmp_transfer".format(output_directory)
		expr_path  = "{}_expr.csv".format(tmp_prefix)
		meta_path  = "{}_metadata.csv".format(tmp_prefix)

		if scipy.sparse.issparse(adata.X):
			expr_df = pd.DataFrame(
				adata.X.toarray(), index=adata.obs_names, columns=adata.var_names
			)
		else:
			expr_df = pd.DataFrame(
				adata.X, index=adata.obs_names, columns=adata.var_names
			)

		expr_df.T.to_csv(expr_path)
		adata.obs.to_csv(meta_path)

		r_script = """
library(Seurat)

counts   <- read.csv('{expr}', row.names=1, check.names=FALSE)
metadata <- read.csv('{meta}', row.names=1, check.names=FALSE)

counts_sparse <- as(as.matrix(counts), "dgCMatrix")
seurat_obj    <- CreateSeuratObject(counts = counts_sparse, meta.data = metadata)

if ('cell_type' %in% colnames(seurat_obj@meta.data)) {{
    Idents(seurat_obj) <- 'cell_type'
}}

saveRDS(seurat_obj, file = '{rds}')
""".format(expr=expr_path, meta=meta_path, rds=output.rds)

		r_script_path = "{}_generator.R".format(tmp_prefix)
		with open(r_script_path, "w") as _fh:
			_fh.write(r_script)

		subprocess.run(["Rscript", r_script_path], check=True)

		for _f in [expr_path, meta_path, r_script_path]:
			if os.path.exists(_f):
				os.remove(_f)

		print("Seurat conversion complete.")
