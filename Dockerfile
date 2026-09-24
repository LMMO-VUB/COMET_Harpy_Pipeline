# Dockerfile for the COMET Spatial Pipeline (image: spatial_pipeline:v1).
# Built via `docker compose build` / `docker compose up --build`, driven by
# docker-compose.yml in this same folder. This file must live in the repo
# (not just as a locally-saved image/tar) so the image can be reproduced from
# source alone.
#
# Use a stable Ubuntu base to support both Python and R cleanly
FROM ubuntu:22.04

# Prevent interactive prompts during package installations
ENV DEBIAN_FRONTEND=noninteractive

# Install system dependencies required for spatial computing, graphics, and R compiling
RUN apt-get update && apt-get install -y --no-install-recommends \
    wget \
    curl \
    git \
    ca-certificates \
    build-essential \
    libgl1-mesa-glx \
    libglib2.0-0 \
    gfortran \
    liblapack-dev \
    libopenblas-dev \
    libxml2-dev \
    libssl-dev \
    libcurl4-openssl-dev \
    && rm -rf /var/lib/apt/lists/*

# Install Miniforge3 (which natively includes Mamba) to manage our unified bioinformatics stack
RUN wget https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh -O miniforge.sh \
    && bash miniforge.sh -b -p /opt/conda \
    && rm miniforge.sh

# Put conda/mamba in the system path environment
ENV PATH=/opt/conda/bin:$PATH

# Configure channels for reproducible bioinformatics resolving
RUN conda config --add channels defaults \
    && conda config --add channels bioconda \
    && conda config --add channels conda-forge \
    && conda config --set channel_priority strict

# Install exact Python and spatial math library foundations from your Mac history.
#
# scanpy/umap-learn/scikit-learn/igraph, and especially pynndescent, numba, and
# leidenalg, are pinned explicitly here (not just left as transitive deps of
# squidpy) because a version mismatch in exactly these three packages between
# this image and the Mac's sproteo_fresh env was confirmed (2026-09-24) to
# produce genuinely different Leiden cluster counts at the same
# clustering_resolution, and a visibly different UMAP -- not just cosmetic
# layout noise. pynndescent builds the approximate neighbor graph Leiden
# clusters on, and leidenalg is the clustering implementation itself, so a
# minor-version drift in either is enough to change results; numba is what
# both are JIT-compiled through. Versions below match the Mac exactly as of
# that date -- if you ever update the Mac env, update these to match (and
# vice versa), don't let them drift independently again.
RUN mamba install -y \
    python=3.10 \
    numpy=1.26 \
    scipy=1.11 \
    scikit-image \
    squidpy=1.6.5 \
    scanpy=1.11.5 \
    umap-learn=0.5.12 \
    pynndescent=0.5.13 \
    scikit-learn=1.7.2 \
    numba=0.65.1 \
    leidenalg=0.11.0 \
    python-igraph=1.0.0 \
    snakemake \
    && mamba clean --all -y

# Install R and Seurat directly via Mamba to guarantee binary library alignments.
#
# NOTE on version pin: this was briefly bumped to r-seurat=5.1.0 during a code
# review, on the theory that convert_to_seurat's generated R script -- which
# calls SetAssayData(..., layer = "data", ...), a SeuratObject-5-only argument
# name (renamed from "slot") -- would fail against Seurat 4.x with "unused
# argument (layer = ...)". That was reverted: confirmed runs on the Windows
# machine, using an image built from this same 4.4.0 pin (loaded via
# spatial_pipeline_v1.tar), have completed convert_to_seurat successfully. So
# either SeuratObject added `layer` as a pre-5.0 forward-compat alias, or
# something else about that theory was wrong -- either way, real observed
# behavior beats the theory. Left at 4.4.0. If you ever do bump this, retest
# the convert_to_seurat rule specifically (it's the last rule in the DAG, so a
# failure there is easy to miss until the very end of a run) before relying on
# the new version. Note envs/harpy_r_env.yaml (the macOS reference env) still
# lists r-seurat=5.1.0, so the two platforms are running different Seurat
# majors -- worth reconciling deliberately rather than by accident.
RUN mamba install -y \
    r-base=4.3 \
    r-seurat=4.4.0 \
    && mamba clean --all -y

# Install exact PyPI spatial and dashboard packages from your Mac environment
RUN pip install --no-cache-dir \
    anndata==0.11.4 \
    streamlit==1.58.0 \
    harpy-analysis==0.3.0

# harpy-analysis only formally requires leidenalg>=0.9.1 (which 0.11.0 already
# satisfies), but pip's resolver still re-pinned leidenalg down to 0.10.2 while
# satisfying something in harpy's much larger transitive dependency tree
# (spatialdata-plot, flowsom, scanpy, etc. all ship their own constraints).
# Confirmed via `docker run --rm spatial_pipeline:v1 python -c "import
# leidenalg; print(leidenalg.version)"` printing 0.10.2 right after a build
# that explicitly mamba-installed 0.11.0 above -- pip silently overrode it.
#
# leidenalg's own PyPI metadata declares python-igraph as a dependency, so
# pulling in a replacement leidenalg drags a replacement python-igraph along
# with it too -- confirmed the hard way: forcing ONLY leidenalg back via
# `pip install --no-deps --force-reinstall leidenalg==0.11.0` still left
# `igraph.__version__` reporting "0.11.9" instead of the conda-forge-correct
# "1.0.0", because the pip install of harpy-analysis above had already
# replaced python-igraph before this step ever ran, and --no-deps on
# leidenalg alone doesn't touch a package that's already sitting there.
# A prior attempt at fixing this also tried `pip install --no-deps
# --force-reinstall leidenalg==0.11.0 python-igraph==1.0.0` together, which
# didn't work either: PyPI's "python-igraph" package at tag 1.0.0 itself
# reports `igraph.__version__ == "0.11.9"` at runtime -- a real
# inconsistency between how conda-forge and PyPI version this package. So
# there's no pip incantation that gets python-igraph back to a correctly
# self-reporting 1.0.0.
#
# The reliable fix is to stop fighting pip's resolver with more pip and
# instead re-run mamba -- the same tool and channel that installed a
# confirmed-correct build of all of these the first time (the "Install
# exact Python and spatial math library foundations" step above).
#
# A plain `mamba install --force-reinstall leidenalg=0.11.0
# python-igraph=1.0.0` here was ALSO tried and ALSO didn't work (confirmed
# on Windows: still reported 0.10.2 / 0.11.9 afterward, identical to what
# pip had left behind). The likely cause: pip overwrote the on-disk files
# for these packages without updating conda's own installed-package
# records, so conda's records still say "leidenalg 0.11.0 is already
# installed here" and --force-reinstall's usual linking shortcuts can skip
# actually re-copying files it believes are already correct, leaving pip's
# files in place underneath. So this step first deletes the on-disk
# package directories directly (bypassing both pip's and conda's
# bookkeeping, which by this point disagree with each other and with
# reality) before asking mamba for a genuinely fresh install into empty
# space, and then asserts the result at build time instead of trusting it
# -- if this ever regresses again, the build fails immediately with the
# actual versions printed, rather than silently shipping a bad image that
# only gets caught hours later during a pipeline run.
RUN SITE=$(python -c "import sysconfig; print(sysconfig.get_paths()['purelib'])") \
    && (pip uninstall -y leidenalg python-igraph igraph 2>/dev/null || true) \
    && rm -rf "$SITE"/leidenalg* "$SITE"/igraph* "$SITE"/python_igraph* \
    && mamba install -y --force-reinstall \
        scanpy=1.11.5 \
        umap-learn=0.5.12 \
        pynndescent=0.5.13 \
        scikit-learn=1.7.2 \
        numba=0.65.1 \
        leidenalg=0.11.0 \
        python-igraph=1.0.0 \
    && mamba clean --all -y \
    && python -c "\
import leidenalg, igraph; \
print('VERSION CHECK: leidenalg=' + leidenalg.version + ' igraph=' + igraph.__version__); \
assert leidenalg.version == '0.11.0', 'leidenalg is ' + leidenalg.version + ', expected 0.11.0'; \
assert igraph.__version__ == '1.0.0', 'igraph is ' + igraph.__version__ + ', expected 1.0.0'"

# Expose the standard port used by your Streamlit interactive cell-type assignment interface
EXPOSE 8501

# Set up a working directory mapping inside the container
WORKDIR /data

# Default action when launching the container is a bash terminal prompt
CMD ["/bin/bash"]