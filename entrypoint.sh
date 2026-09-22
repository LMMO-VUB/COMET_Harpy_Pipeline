#!/usr/bin/env bash
# entrypoint.sh  –  COMET Spatial Pipeline (Docker)
#
# Mount layout expected by docker-compose.yml:
#   /pipeline  →  snake_project directory (Snakefile, scripts, envs …) [read-only]
#   /data      →  user's run folder (config.yaml, input CSV, outputs)  [read-write]

set -euo pipefail

echo "[entrypoint] Pipeline dir : /pipeline"
echo "[entrypoint] Data dir     : /data"
echo "[entrypoint] Config       : /data/config.yaml"

# ── Seed metadata_markers.csv into the run folder if not already there.
#    The user may have placed a custom one; if so, respect it.
if [ ! -f "/data/metadata_markers.csv" ]; then
    if [ -f "/pipeline/metadata_markers.csv" ]; then
        echo "[entrypoint] Copying metadata_markers.csv from pipeline directory…"
        cp /pipeline/metadata_markers.csv /data/metadata_markers.csv
    else
        echo "[entrypoint] WARNING: metadata_markers.csv not found in /pipeline or /data." \
             "The pipeline will attempt to auto-generate it."
    fi
fi

# ── Run Snakemake.
#    --directory /data  makes all relative paths in the config (HALO_data_file,
#    output_directory, metadata_file) resolve under /data, which is writable.
#    --snakefile points to the actual Snakefile in the read-only pipeline mount.
echo "[entrypoint] Starting Snakemake…"
snakemake \
    --snakefile  /pipeline/Snakefile \
    --configfile /data/config.yaml \
    --directory  /data \
    --cores      all \
    --rerun-incomplete
