#!/usr/bin/env Rscript
#
# merge_tma_seurat.R
#
# Combine several per-TMA Seurat .rds outputs from the COMET Spatial
# Pipeline (one produced per run of rule convert_to_seurat, e.g.
# "<project>_results/<project>_annotated_seurat.rds") into a single Seurat
# object for joint downstream analysis.
#
# WHY SPLIT BY TMA IN THE FIRST PLACE:
#   Each HALO/Horizon export has its own local pixel/micron coordinate
#   space, and the pipeline numbers cells sequentially within one export
#   (Object_Id_Clean). Feeding six TMAs through as one combined CSV would
#   silently collide cell IDs and mix unrelated spatial coordinate systems
#   together (bad for rasterization and neighborhood enrichment, which
#   assume one coherent tissue section). Running each TMA through the
#   pipeline separately, then merging the resulting Seurat objects
#   afterward, is the right order of operations.
#
# WHAT THIS SCRIPT DOES:
#   - Loads each TMA's Seurat .rds
#   - Tags every cell with which TMA it came from (a new "TMA" metadata
#     column), and prefixes cell IDs with that label so cells from
#     different TMAs can never collide once merged
#   - Merges all objects into one via Seurat::merge()
#   - Drops the per-TMA "pca"/"umap" reductions that came along for the
#     ride, since each was computed independently on that TMA's cells alone
#     -- they do NOT share a coordinate space, so plotting them together
#     post-merge would be meaningless. Prints how to recompute a single,
#     shared PCA/UMAP on the combined object if you want one.
#
# IMPORTANT CAVEAT -- READ THIS BEFORE TRUSTING A MERGED "cell_type":
#   Leiden cluster numbers (and skip_annotation_portal's cluster-ID-as-
#   cell_type auto-assignment) are NOT comparable across separate pipeline
#   runs. Cluster "3" in TMA1 and cluster "3" in TMA3 can be two completely
#   different cell populations -- each TMA's clustering is computed
#   independently, so the numbers just happen to repeat. If you want a
#   single, meaningfully comparable "cell_type" across all 6 TMAs, you
#   must go through the manual annotation portal for every TMA (i.e. do
#   NOT set skip_annotation_portal: true for this workflow) and use the
#   SAME cell-type names every time (e.g. always "CD8 T cell", never
#   sometimes "CD8+ T cell" / "CD8 Tcell") -- the portal's job is exactly
#   to translate each TMA's arbitrary cluster numbers into one shared,
#   human-chosen vocabulary. This script does not and cannot fix
#   inconsistent naming for you.
#
# USAGE:
#   Rscript merge_tma_seurat.R --out combined_seurat.rds \
#     TMA1=/path/to/TMA1_results/TMA1_annotated_seurat.rds \
#     TMA2=/path/to/TMA2_results/TMA2_annotated_seurat.rds \
#     TMA3=/path/to/TMA3_results/TMA3_annotated_seurat.rds \
#     TMA4=/path/to/TMA4_results/TMA4_annotated_seurat.rds \
#     TMA5=/path/to/TMA5_results/TMA5_annotated_seurat.rds \
#     TMA6=/path/to/TMA6_results/TMA6_annotated_seurat.rds
#
#   Works with any number >= 2 of <label>=<path> pairs, not just 6 --
#   label is whatever short, unique name you want each TMA tagged with
#   (shows up in the new "TMA" metadata column and as a cell-id prefix).

suppressPackageStartupMessages(library(Seurat))

args <- commandArgs(trailingOnly = TRUE)
out_path <- "combined_seurat.rds"
pairs <- list()

i <- 1
while (i <= length(args)) {
  a <- args[i]
  if (a == "--out") {
    if (i == length(args)) stop("--out needs a value", call. = FALSE)
    out_path <- args[i + 1]
    i <- i + 2
  } else if (grepl("=", a, fixed = TRUE)) {
    kv <- strsplit(a, "=", fixed = TRUE)[[1]]
    label <- kv[1]
    path  <- paste(kv[-1], collapse = "=")
    if (label %in% names(pairs)) {
      stop("Duplicate TMA label: '", label, "'. Labels must be unique.", call. = FALSE)
    }
    pairs[[label]] <- path
    i <- i + 1
  } else {
    stop(
      "Unrecognized argument: '", a, "'. Expected --out <path> or ",
      "<label>=<path_to_seurat.rds>.", call. = FALSE
    )
  }
}

if (length(pairs) < 2) {
  stop(
    "Need at least two <label>=<path_to_seurat.rds> pairs to merge.\n",
    "Example:\n",
    "  Rscript merge_tma_seurat.R --out combined.rds \\\n",
    "    TMA1=.../TMA1_results/TMA1_annotated_seurat.rds \\\n",
    "    TMA2=.../TMA2_results/TMA2_annotated_seurat.rds",
    call. = FALSE
  )
}

labels <- names(pairs)

cat("Loading", length(pairs), "Seurat object(s):\n")
objs <- list()
for (label in labels) {
  path <- pairs[[label]]
  if (!file.exists(path)) stop("File not found for '", label, "': ", path, call. = FALSE)
  cat(sprintf("  %-10s <- %s\n", label, path))
  obj <- readRDS(path)
  obj$TMA <- label  # record which TMA each cell came from, for grouping/splitting later
  objs[[label]] <- obj
}

cat("\nMerging (prefixing cell IDs with each TMA label to avoid collisions)...\n")
combined <- merge(
  x            = objs[[1]],
  y            = objs[-1],
  add.cell.ids = labels,
  project      = "COMET_combined"
)

cat(sprintf(
  "\nCombined object: %d cells across %d TMA(s), assays: %s\n",
  ncol(combined), length(labels), paste(Assays(combined), collapse = ", ")
))

if ("cell_type" %in% colnames(combined@meta.data)) {
  cat("\ncell_type counts per TMA (sanity-check your naming consistency here --\n")
  cat("the same biological label should show up with similar proportions across\n")
  cat("TMAs from the same tissue/condition; wildly different label SETS per TMA\n")
  cat("usually means the portal's naming wasn't kept consistent):\n")
  print(table(combined$cell_type, combined$TMA))
  Idents(combined) <- "cell_type"
} else {
  cat(
    "\nNote: no 'cell_type' column found. If skip_annotation_portal was used,\n",
    "this will instead be whatever raw Leiden cluster IDs each TMA produced --\n",
    "see the caveat at the top of this script before treating those as\n",
    "comparable across TMAs.\n",
    sep = ""
  )
}

# Each input object's "pca"/"umap" reduction was computed independently, on
# that TMA's cells alone -- they are different coordinate spaces that just
# happen to share axis names. Keeping them around post-merge invites
# plotting cells from different TMAs on axes that don't mean the same thing
# for each of them, so drop them here rather than leave that trap in place.
dropped <- character(0)
for (reduction in c("pca", "umap")) {
  if (reduction %in% Reductions(combined)) {
    combined[[reduction]] <- NULL
    dropped <- c(dropped, reduction)
  }
}
if (length(dropped) > 0) {
  cat(sprintf(
    "\nDropped per-TMA reduction(s) not shared across TMAs: %s\n",
    paste(dropped, collapse = ", ")
  ))
}

saveRDS(combined, file = out_path)
cat(sprintf("\nSaved combined Seurat object -> %s\n", out_path))

cat(
  "\nThis object now has one shared 'TMA' column and (if you used the\n",
  "annotation portal consistently) one shared 'cell_type' column, but no\n",
  "shared PCA/UMAP embedding yet. For a single joint embedding across all\n",
  "TMAs -- e.g. for one combined UMAP figure, or before integration/batch\n",
  "correction -- run something like (adjust npcs/dims to match pca_dims in\n",
  "your config.yaml):\n",
  "    DefaultAssay(combined) <- \"Protein\"\n",
  "    combined <- ScaleData(combined)\n",
  "    combined <- RunPCA(combined, npcs = 5)\n",
  "    combined <- RunUMAP(combined, dims = 1:5)\n",
  "(The Protein assay's 'data' layer was already set by the pipeline to its\n",
  "normalized values, so you likely do NOT want to call NormalizeData()\n",
  "again here -- that would double-normalize.)\n",
  sep = ""
)
