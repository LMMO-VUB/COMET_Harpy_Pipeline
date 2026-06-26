
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


######################################
######################################
######################################
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["BLAS_NUM_THREADS"] = "1"
os.environ["VECLIB_MAXIMUM_THREADS"] = "1"
os.environ["NUMEXPR_NUM_THREADS"] = "1"
######################################

rule all:
	input:
		f"{output_directory}/{name_of_project}_annotated_seurat.rds",
		f"{output_directory}/images_final/umap_cell_types.png"

rule process_halo:
	input:
		csv = HALO_data_file,
		meta = metadata_file
	output:
		h5ad = f"{output_directory}/{name_of_project}.h5ad",
	run:
		import os
		import pandas as pd
		import numpy as np
		import scanpy as sc
		import anndata as ad
		import spatialdata as sd
		import squidpy as sq
		import geopandas as gpd
		import matplotlib
		matplotlib.use('Agg')
		import matplotlib.pyplot as plt
	 	import harpy as hp
		import re
		from shapely.geometry import Point
		from spatialdata.models import ShapesModel, TableModel
		from spatialdata.transformations import Identity

		sc.settings.autoshow = False

		os.makedirs(f"{output_directory}", exist_ok=True)
		os.makedirs(f"{output_directory}/images", exist_ok=True)

		# Load raw HALO data
		df = pd.read_csv(input.csv, decimal=',')
		markers_df = pd.read_csv(input.meta, index_col=False, skiprows=[0,1])

		########## Double-check the markers_df is formatted correctly
		### Error checkin
		if not all(markers_df.type.isin(['Protein', 'Transcript', 'Other'])):
			var = markers_df.type[markers_df.type.isin(['Protein', 'Transcript', 'Other'])]
			raise ValueError(f"Unkown value in 'type' column of metadata: {var}")
		if not all(markers_df.localization.isin(["Nucleus", "Cytoplasm", "NA"])):
			var = markers_df.localization[markers_df.localization.isin(["Nucleus", "Cytoplasm", "NA"])]
			raise ValueError(f"Unkown value in 'type' column of metadata: {var}")

		# Build the feature set based on the CSV instructions
		# Example: HALO column might be "CD3 Cytoplasm Intensity"


		###### Create some new columns
		df.columns = [c.replace('?', 'alpha_').strip() for c in df.columns]
		# Calculate Centroids
		df['x'] = (df['XMin'] + df['XMax']) / 2
		df['y'] = (df['YMin'] + df['YMax']) / 2
		df['Object_Id_Clean'] = df['Object Id'].astype(int) + 1
		# Convert Area (µm²) to Area (pixels²) then calculate Radius in pixels
		pixel_size = 0.5 
		df['cell_area_px'] = df['Cell Area (µm²)'] / (pixel_size**2)
		df['radius_px'] = np.sqrt(df['cell_area_px'] / np.pi)


		def clean_name(name):
			### A Function which cleans up the name of the channel intensities because they seem
			### to come out of HALO inconsistent and with a bunch of illegal characters
			name = name.replace('?', 'a').replace('α', 'a')

			cycle_match = re.search(r'(\d+)\s*\|', name)
			cycle_prefix = f"{cycle_match.group(1)}_" if cycle_match else ""
			clean = name.split('|')[-1].strip()
			clean = clean.replace('T2 ', 'T2_').replace('T3 ', 'T3_')

			# 5. Remove all the various intensity/area suffixes
			suffixes = [
				' Nucleus Intensity', ' Cell Intensity', ' Avg Intensity', 
				' Intensity', ' Positive', ' Positive Nucleus', 
				' Positive Cytoplasm', ' Cytoplasm Intensity', ' Copies'
			]
			for s in suffixes:
				clean = clean.replace(s, '')
			clean = re.sub(r'_\d+x', '', clean)
			clean = re.sub(r'_\d+\s\d+x', '', clean) # handles the space in 20 000x
			clean = re.sub(r'[^a-zA-Z0-9_]', '_', clean).strip('_')
			clean = re.sub(r'\s+', '_', clean).strip('_')

			return f"{cycle_prefix}{clean}"


		### Selects exactly which intensity columns we are interested in
		channel_columns = {}
		halo_cols=[]
		print(df.columns)
		for _, row in markers_df.iterrows():
			if (row['type']=='Protein') or (row['type']=='Other'):
				col_name = f"{row['channel']} {row['localization']} Intensity"
			elif row['type']=='Transcript':
				col_name = f"{row['channel']} Copies"

			### Extract the column of interest
			for halo_col in df.columns:
				if col_name in halo_col:
					print(col_name)
					print(halo_col)
					channel_columns[row['channel']] = df[halo_col]
					halo_cols.append(halo_col)
		clean_marker_names = [clean_name(c) for c in halo_cols]

		### Make Anndata
		adata = ad.AnnData(
			X=df[halo_cols].values.astype(np.float32),
			obs=df[['Object_Id_Clean']].copy(),
			var=markers_df.set_index(pd.Index(clean_marker_names))
		)
		adata.layers["raw"] = adata.X.copy()
		adata.obs['region'] = 'halo_labels'
		adata.obs['region'] = adata.obs['region'].astype('category')
		adata.obs['instance_id'] = df['Object_Id_Clean'].astype(int).values
		adata.obs_names = df['Object_Id_Clean'].astype(str).values
		adata.obs['cell_area_um2'] = df['Cell Area (µm²)'].values
		"Cell Area (µm²)"
		spatial_coords = df[['x', 'y']].values.astype(np.float32)
		adata.obsm['spatial'] = spatial_coords
		adata.var['feature_type'] = adata.var['type']
	        
		

		####################
		###### Make Sdata
		sdata = sd.SpatialData()
		
		# Shapes Slot
		geometry = [Point(xy) for xy in zip(df['x'], df['y'])]
		gdf = gpd.GeoDataFrame(df, geometry=geometry)
		# Buffer using the PIXEL radius to ensure alignment with TIFF
		gdf['geometry'] = gdf.buffer(df['radius_px'], cap_style=1)
		gdf.index = df['Object_Id_Clean'].astype(str).values
		# Parse into SpatialData ShapesModel
		transform = Identity()
		cells_shapes = ShapesModel.parse(
			gdf, 
			transformations={'global': transform}
		)
		sdata.shapes['halo_cells'] = cells_shapes        
	        
		# Labels Slot
		target_shape = (44643, 44643)
		hp.im.rasterize(
			sdata=sdata,
			shapes_layer='halo_cells',
			output_layer='halo_labels',
			out_shape=target_shape,
			chunks=2048, # Chunking improves performance for large COMET images
			overwrite=True
		)

		# Table Slot
		sdata.table = TableModel.parse(
			adata,
			region="halo_labels",
			region_key="region",
			instance_key="instance_id"
		)

		####### Sdata Made
		#######################




		########################
		######## Data Processing
		# Normalization
		print("Normalizing")
		is_protein = adata.var['feature_type'] == 'Protein'
		is_rna = adata.var['feature_type'] == 'Transcript'

		adata.layers['normalized'] = adata.layers["raw"]
		adata.layers['normalized'][:,is_protein] = np.arcsinh(adata.layers['raw'][:, is_protein] / 50.0)
		adata.layers['normalized'][:,is_rna] = np.log1p(adata.layers['raw'][:, is_rna])
		adata.X = adata.layers['normalized'].copy()

		if not isinstance(adata.X, np.ndarray):
			adata.X = adata.X.toarray()  # Ensure dense matrix
		adata.X = np.ascontiguousarray(adata.X, dtype=np.float32)

		print("Scaling")
		sc.pp.scale(adata, max_value=10)

		clustering_markers = [m for m in sdata.table.var_names if (sdata.table.var.loc[m, 'feature_type']=='Protein')]
		adata_subset = sdata.table[:, clustering_markers].copy()

		# Dimensional Reduction
		print("Performing PCA")
		sc.tl.pca(adata_subset, n_comps= min(10, adata_subset.shape[1]) )
		sc.pl.pca_variance_ratio(adata_subset, n_pcs=min(10, adata_subset.shape[1]), show=False)
		plt.savefig(f"{output_directory}/images/pca_elbowPlot.png", bbox_inches='tight')
		plt.close()
		print("Finding Neighbors")
		sc.pp.neighbors(adata_subset, n_pcs=pca_dims)
		print("Performing UMAP")
		sc.tl.umap(adata_subset)

		sdata.table.obsm['X_pca'] = adata_subset.obsm['X_pca']
		sdata.table.obsm['X_umap'] = adata_subset.obsm['X_umap']
		sdata.table.uns['neighbors'] = adata_subset.uns['neighbors']
		sdata.table.obsp['distances'] = adata_subset.obsp['distances']
		sdata.table.obsp['connectivities'] = adata_subset.obsp['connectivities']
		full_pcs = np.zeros((sdata.table.shape[1], adata_subset.varm['PCs'].shape[1]))
		var_indices = [sdata.table.var_names.get_loc(m) for m in clustering_markers]
		full_pcs[var_indices, :] = adata_subset.varm['PCs']
		sdata.table.varm['PCs'] = full_pcs

		# Clustering
		print("Performing Leiden Clustering")
		sc.tl.leiden(sdata.table, resolution=clustering_resolution, key_added='leiden_clusters')

		sc.pl.umap(sdata.table, color='cell_area_um2', show=False)
		plt.savefig(f"{output_directory}/images/umap_cell_areas.png", bbox_inches='tight')
		plt.close()
		sc.pl.umap(sdata.table, color='leiden_clusters', show=False)
		plt.savefig(f"{output_directory}/images/umap_leiden_clustering_res{clustering_resolution}.png", bbox_inches='tight')
		plt.close()

		for var in sdata.table.var_names:
			sc.pl.umap(sdata.table, color=var, show=False)
			plt.savefig(f"{output_directory}/images/umap_{var}.png", bbox_inches='tight')
			plt.close()

		####### Make lots of DotPlots
		sc.pl.dotplot(sdata.table, var_names=sdata.table.var_names, 
	    		groupby='leiden_clusters', show=False)
		plt.savefig(f"{output_directory}/images/dotplot_allChannels.png", bbox_inches='tight')
		plt.close()

		protein_markers = [m for m in sdata.table.var_names if (sdata.table.var.loc[m, 'feature_type']=='Protein')]
		if len(protein_markers)>0:
			sc.pl.dotplot(sdata.table, var_names=protein_markers, 
					groupby='leiden_clusters', show=False)
			plt.savefig(f"{output_directory}/images/dotplot_proteinChannels.png", bbox_inches='tight')
			plt.close()

		rna_markers = [m for m in sdata.table.var_names if (sdata.table.var.loc[m, 'feature_type']=='Transcript')]
		if len(rna_markers)>0:
			sc.pl.dotplot(sdata.table, var_names=rna_markers, 
					groupby='leiden_clusters', show=False)
			plt.savefig(f"{output_directory}/images/dotplot_rnaChannels.png", bbox_inches='tight')
			plt.close()

		print("Performing Neighborhood Enrichment Analysis")
		### Neighborhood Enrichment
		hp.tb.nhood_enrichment(sdata, 
			labels_layer='halo_labels', 
			table_layer='table', 
			output_layer='table_score_genes_enrichment', 
			celltype_column='leiden_clusters'
		)
		hp.pl.nhood_enrichment(
			sdata,
			table_layer="table_score_genes_enrichment",
			celltype_column='leiden_clusters',
			output=f'{output_directory}/images/nhood_enrichment.png'
		)
		#adata.write_h5ad(f"{output.directory}/{name_of_project}.h5ad")

		print("Saving Output")
		adata.obs['x'] = df['x'].values.astype(np.float32)
		adata.obs['y'] = df['y'].values.astype(np.float32)
		adata.obs['radius_px'] = df['radius_px'].values.astype(np.float32)

		final_adata = sdata.table.copy() 
		final_adata.write_h5ad(output.h5ad)

		print("Finished Successfully!")


rule annotate_clusters:
	input:
		h5ad = f"{output_directory}/{name_of_project}.h5ad"
	output:
		annotated_h5ad = f"{output_directory}/{name_of_project}_annotated.h5ad"
	run:
		import subprocess
		import sys
		import os
		import time

		print("\n" + "="*60)
		print("PIPELINE PAUSED: LAUNCHING INTERACTIVE ANNOTATION PORTAL")
		print("Opening your browser. Complete annotations and click 'Save' to resume.")
		print("="*60 + "\n")

		# Launch streamlit as a background process
		# Passing the target file paths so the app knows where to save
		cmd = [
			"streamlit", "run", "app.py", "--", 
			"--input", input.h5ad, 
			"--output", output.annotated_h5ad,
			"--img_dir", f"{output_directory}/images"
		]
        
		process = subprocess.Popen(cmd)

		# Loop and block Snakemake until the researcher creates the final file
		while not os.path.exists(output.annotated_h5ad):
			time.sleep(2)

		# Once the file is detected, terminate the Streamlit server cleanly
		print("\nAnnotations saved successfully! Resuming Snakemake pipeline...")
		process.terminate()
		process.wait()



rule post_annotation_viz:
	input:
		h5ad = f"{output_directory}/{name_of_project}_annotated.h5ad"
	output:
		umap_img = f"{output_directory}/images_final/umap_cell_types.png",
		dot_img = f"{output_directory}/images_final/dotplot_final_proteins.png",
		nhood_img = f"{output_directory}/images_final/nhood_enrichment_final.png"
	run:
		import matplotlib
		matplotlib.use('Agg')
		import scanpy as sc
		import harpy as hp
		import matplotlib.pyplot as plt
		import spatialdata as sd
		import geopandas as gpd
		
		# MAKE SURE THIS LINE IS PRESENT INSIDE THIS RUN BLOCK:
		from shapely.geometry import Point 
		
		from spatialdata.models import ShapesModel, TableModel
		from spatialdata.transformations import Identity
		import os

		os.makedirs(os.path.dirname(output.umap_img), exist_ok=True)
		adata = sc.read_h5ad(input.h5ad)

		# Remake UMAPs with Annotated cell types
		sc.pl.umap(adata, color='cell_type', show=False)
		plt.savefig(output.umap_img, bbox_inches='tight')
		plt.close()

		# Remake Marker DotPlots
		protein_markers = [m for m in adata.var_names if adata.var.loc[m, 'feature_type']=='Protein']
		if protein_markers:
			sc.pl.dotplot(adata, var_names=protein_markers, groupby='cell_type', show=False)
			plt.savefig(output.dot_img, bbox_inches='tight')
			plt.close()

		# 3. Neighborhood Enrichment with Biological Names
		# We need to reload the SpatialData object to link the table back
		# Note: Ensure sdata was saved or re-construct for spatial analysis
		# Assuming 'sdata' structure is accessible or re-loaded
		sdata = sd.SpatialData() 

		geometry = [Point(xy) for xy in zip(adata.obs['x'], adata.obs['y'])]
		gdf = gpd.GeoDataFrame(adata.obs, geometry=geometry)
		gdf['geometry'] = gdf.buffer(adata.obs['radius_px'], cap_style=1)
		gdf.index = adata.obs['instance_id'].astype(str).values
		
		# Map shapes and link your annotated data frame
		cells_shapes = ShapesModel.parse(gdf, transformations={'global': Identity()})
		sdata.shapes['halo_cells'] = cells_shapes 
		
		# Rasterize into the tracking target label name expected by harpy
		target_shape = (44643, 44643)
		hp.im.rasterize(
			sdata=sdata,
			shapes_layer='halo_cells',
			output_layer='halo_labels',
			out_shape=target_shape,
			chunks=2048,
			overwrite=True
		)
		
		if 'spatialdata_attrs' in adata.uns:
			del adata.uns['spatialdata_attrs']

		sdata.table = TableModel.parse(
			adata,
			region="halo_labels",
			region_key="region",
			instance_key="instance_id"
		)

		hp.tb.nhood_enrichment(
			sdata, 
			labels_layer='halo_labels',    
			table_layer='table', 
			celltype_column='cell_type',
			output_layer='table_annotated_enrichment'
		)
		
		hp.pl.nhood_enrichment(
			sdata,
			table_layer="table_annotated_enrichment",
			celltype_column='cell_type',
			output=output.nhood_img
		)


		print("Final annotated visualizations generated successfully!")



rule convert_to_seurat:
    input:
        h5ad = f"{output_directory}/{name_of_project}_annotated.h5ad"
    output:
        rds = f"{output_directory}/{name_of_project}_annotated_seurat.rds"
    run:
        import scanpy as sc
        import pandas as pd
        import subprocess
        import os
        import scipy.sparse

        print("Converting anndata to Seurat object")
        
        adata = sc.read_h5ad(input.h5ad)
        
        tmp_prefix = f"{output_directory}/tmp_transfer"
        expr_path = f"{tmp_prefix}_expr.csv"
        meta_path = f"{tmp_prefix}_metadata.csv"

        if scipy.sparse.issparse(adata.X):
            expr_df = pd.DataFrame(adata.X.toarray(), index=adata.obs_names, columns=adata.var_names)
        else:
            expr_df = pd.DataFrame(adata.X, index=adata.obs_names, columns=adata.var_names)
            
        # Temp files
        expr_df.T.to_csv(expr_path)
        adata.obs.to_csv(meta_path)
        
        r_payload = f"""
        library(Seurat)
        
        # Load matrices using native csv parsers
        counts <- read.csv('{expr_path}', row.names=1, check.names=FALSE)
        metadata <- read.csv('{meta_path}', row.names=1, check.names=FALSE)
        
        # Safely convert to a sparse matrix format inside R to optimize memory footprint
        counts_sparse <- as(as.matrix(counts), "dgCMatrix")
        
        # Construct the final Seurat Object
        seurat_obj <- CreateSeuratObject(counts = counts_sparse, meta.data = metadata)
        
        if ('cell_type' %in% colnames(seurat_obj@meta.data)) {{
            Idents(seurat_obj) <- 'cell_type'
        }}
        
        saveRDS(seurat_obj, file = '{output.rds}')
        """
        
        r_script_path = f"{tmp_prefix}_generator.R"
        with open(r_script_path, "w") as f:
            f.write(r_payload)
        subprocess.run(["Rscript", r_script_path], check=True)
        
        for f in [expr_path, meta_path, r_script_path]:
            if os.path.exists(f):
                os.remove(f)
                
        print("Seurat RDS File Conversion Success")
