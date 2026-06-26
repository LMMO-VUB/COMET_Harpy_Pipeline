import streamlit as st
import scanpy as sc
import pandas as pd
import argparse
import os

st.set_page_config(layout="wide", page_title="COMET Cluster Annotator")

# Handle arguments passed from Snakemake
parser = argparse.ArgumentParser()
parser.add_argument("--input", required=True)
parser.add_argument("--output", required=True)
parser.add_argument("--img_dir", required=True)
args, unknown = parser.parse_known_args()

st.title("Multi-Omics Cluster Annotation Portal")
st.markdown("Inspect your diagnostic plots on the right to accurately name the clusters on the left.")
st.divider()

# Load Data securely
@st.cache_resource
def load_spatial_data(path):
    return sc.read_h5ad(path)

adata = load_spatial_data(args.input)
clusters = sorted(adata.obs['leiden_clusters'].unique(), key=int)

# Use session state to maintain the base configuration safely
if "base_table" not in st.session_state:
    st.session_state.base_table = pd.DataFrame({
        "Cluster ID": clusters,
        "Assigned Cell Type": ["" for _ in clusters]
    })

# Main Structural Split: Form Sheet vs Image Dashboard
col_left, col_right = st.columns([3, 7], gap="large")

with col_left:
    st.subheader("Assign Cell Types")
    st.caption("Type labels below. Press Enter to commit, use Down Arrow to navigate, and Enter to edit the next row.")
    
    # Passing the dataframe directly without an on_change callback prevents reruns mid-typing
    # Key parameter ensures structural persistence across navigation events
    edited_output = st.data_editor(
        st.session_state.base_table, 
        disabled=["Cluster ID"], 
        hide_index=True,
        width='stretch',
        num_rows="fixed",
        key="static_editor"
    )

    st.markdown("---")
    if st.button("Finalize & Resume Pipeline", type="primary", width='stretch'):
        # Map labels from the live editor output state
        mapping = dict(zip(edited_output["Cluster ID"].astype(str), edited_output["Assigned Cell Type"]))
        
        if any(str(val).strip() == "" for val in mapping.values()):
            st.error("Please provide an annotation for all clusters before saving.")
        else:
            with st.spinner("Applying annotations and saving dataset..."):
                adata.obs['cell_type'] = adata.obs['leiden_clusters'].map(mapping)
                adata.obs['cell_type'] = adata.obs['cell_type'].astype('category')
                
                adata.write_h5ad(args.output)
                st.success("Success! This window will close and your pipeline will resume momentarily.")
                st.rerun(scope="app")

with col_right:
    st.subheader("Spatial & Diagnostic Plot Dashboard")
    
    img_col1, img_col2 = st.columns(2)
    
    # Permanent Anchor: Left Side Image
    with img_col1:
        st.markdown("**Global Cluster Map**")
        
        if os.path.exists(args.img_dir):
            all_files = os.listdir(args.img_dir)
            umap_file = [f for f in all_files if f.startswith("umap_leiden_clustering")]
            
            if umap_file:
                anchor_path = os.path.join(args.img_dir, umap_file[0])
                st.image(anchor_path, caption=f"Reference: {umap_file[0]}", width='stretch')
            else:
                st.warning("Could not find a file matching 'umap_leiden_clustering_res*.png' in the images directory.")
        else:
            st.error("Images folder directory path not found.")
            
    # Interactive Visualizer: Right Side Image
    with img_col2:
        st.markdown("**Expression Profile Explorer**")
        
        if os.path.exists(args.img_dir):
            available_plots = [
                f for f in os.listdir(args.img_dir) 
                if f.endswith(('.png', '.jpg', '.jpeg')) and not f.startswith("umap_leiden_clustering")
            ]
            
            if available_plots:
                selected_plot = st.selectbox(
                    "Switch marker heatmap / dotplot view:", 
                    options=sorted(available_plots),
                    index=0
                )
                
                img_path = os.path.join(args.img_dir, selected_plot)
                st.image(img_path, caption=f"Viewing: {selected_plot}", width='stretch')
            else:
                st.info("No alternative channel marker plots found.")