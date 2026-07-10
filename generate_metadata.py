#!/usr/bin/env python3
"""
generate_metadata.py

Auto-generates metadata_markers.csv by inspecting the column headers of a
HALO (or Horizon) per-cell export CSV, for use when no metadata_markers.csv
already exists in a researcher's run folder.

Called automatically from the Snakefile before rule `process_halo` runs, but
can also be run by hand:

    python generate_metadata.py --input data/halo_exports/my_export.csv \
                                 --output metadata_markers.csv
"""
import argparse
import csv
import re

# ---------------------------------------------------------------------------
# Known marker -> (type, localization) lookup.
#
# HALO/Horizon exports always contain BOTH a "Nucleus Intensity" and a
# "Cytoplasm Intensity" column for every marker, regardless of that marker's
# actual biology -- so the correct readout channel (e.g. FOXP3 is a nuclear
# transcription factor, CD3 is a surface/cytoplasmic marker) can NOT be
# determined from the header text alone. This lookup captures the lab's
# already-validated panel (mirrored from the existing metadata_markers.csv)
# so those assignments are reused automatically. Extend it as the panel
# grows. Anything not listed here gets a Protein/Cytoplasm default and is
# flagged in the console output for manual review.
# ---------------------------------------------------------------------------
KNOWN_MARKERS = {
    "DAPI":        ("Other",   "Cytoplasm"),
    "CD3":         ("Protein", "Cytoplasm"),
    "CD4":         ("Protein", "Cytoplasm"),
    "CD8":         ("Protein", "Cytoplasm"),
    "FOXP3":       ("Protein", "Nucleus"),
    "CD20":        ("Protein", "Nucleus"),
    "CD56":        ("Protein", "Nucleus"),
    "CD11C":       ("Protein", "Nucleus"),
    "CD68":        ("Protein", "Nucleus"),
    "ASMA":        ("Protein", "Nucleus"),
    "PD-L1":       ("Protein", "Cytoplasm"),
    "CD45":        ("Protein", "Cytoplasm"),
    "PD-1":        ("Protein", "Cytoplasm"),
    "CYTOKERATIN": ("Protein", "Cytoplasm"),
    "PANCK":       ("Protein", "Cytoplasm"),
    "KI-67":       ("Protein", "Cytoplasm"),
    "KI67":        ("Protein", "Cytoplasm"),
}

DEFAULT_PROTEIN_LOCALIZATION = "Cytoplasm"

# HALO/Horizon exports are frequently produced on Windows machines with
# regional (non-UTF-8) settings. Forcing a UTF-8 read of such a file is what
# produces the classic "¬µm¬≤" mojibake instead of "µm²". Try encodings in
# this order and sanity-check the result.
CANDIDATE_ENCODINGS = ["utf-8-sig", "cp1252", "latin-1"]


def get_delimiter(sample_line):
    tabs   = sample_line.count("\t")
    semis  = sample_line.count(";")
    commas = sample_line.count(",")
    if tabs >= semis and tabs >= commas:
        return "\t"
    if semis > commas:
        return ";"
    return ","


def read_header(path):
    last_err = None
    for enc in CANDIDATE_ENCODINGS:
        try:
            with open(path, "r", encoding=enc, newline="") as fh:
                first_line = fh.readline()
                
                clean_line = first_line.replace('"', '')
                
                delimiter = get_delimiter(clean_line)
                # Parse the quote-free string directly split by the delimiter
                header = [h.strip() for h in clean_line.split(delimiter)]
                # -----------------------------------------------------------------

            if any("\ufffd" in h or "¬" in h for h in header):
                last_err = f"encoding '{enc}' produced mojibake"
                continue
            return header, enc
        except (UnicodeDecodeError, StopIteration) as e:
            last_err = e
            continue
    raise ValueError(f"Could not read a clean header row from {path}: {last_err}")

# ---------------------------------------------------------------------------
# Column-name parsers for the export schemas we know about.
# ---------------------------------------------------------------------------

# HALO per-cell FL export, e.g.:
#   "1 | CD3_500x Nucleus Intensity"
#   "1 | CD3_500x Cytoplasm Intensity"
#   "1 | T2 AXL Copies"
HALO_INTENSITY_RE = re.compile(
    r"^(?:\d+\s*\|\s*)?(?P<marker>.+?)\s+"
    r"(Nucleus Intensity|Cytoplasm Intensity|Cell Intensity|Avg Intensity)$"
)
HALO_TRANSCRIPT_RE = re.compile(
    r"^(?:\d+\s*\|\s*)?T\d+\s+(?P<marker>.+?)\s+"
    r"(Copies|Area.*|Classification|Cell Intensity|Avg Intensity)$"
)

# Horizon per-nucleus export, e.g.:
#   "Nuclei/Mean Intensity (CD3_500x - TRITC Protein Autofluo)"
HORIZON_RE = re.compile(
    r"Nuclei/Mean Intensity \((?P<marker>.+?)\s*-\s*\w+ (?:Protein|RNA) Autofluo\)"
)

DILUTION_SUFFIX_RE = re.compile(r"_\d+(?:\s?\d+)?x$", re.IGNORECASE)


def strip_dilution(marker_raw):
    marker = marker_raw.replace("?", "a").replace("\u03b1", "a").strip()
    marker = DILUTION_SUFFIX_RE.sub("", marker)
    return marker.strip()


def detect_markers(header):
    """Returns (markers: {canonical_name: is_transcript}, detected_format)."""
    markers = {}
    fmt = None

    for col in header:
        m = HALO_TRANSCRIPT_RE.match(col)
        if m:
            fmt = fmt or "HALO"
            markers[strip_dilution(m.group("marker"))] = True
            continue
        m = HALO_INTENSITY_RE.match(col)
        if m:
            fmt = fmt or "HALO"
            markers.setdefault(strip_dilution(m.group("marker")), False)
            continue
        m = HORIZON_RE.search(col)
        if m:
            fmt = fmt or "HORIZON"
            markers.setdefault(strip_dilution(m.group("marker")), False)
            continue

    return markers, fmt


def classify(marker, is_transcript):
    key = marker.upper()
    if key in KNOWN_MARKERS:
        return KNOWN_MARKERS[key]
    if is_transcript:
        return ("Transcript", "NA")
    if key == "DAPI":
        return ("Other", "Cytoplasm")
    return ("Protein", DEFAULT_PROTEIN_LOCALIZATION)


def generate(input_csv, output_csv):
    input_csv = str(input_csv)
    output_csv = str(output_csv)

    header, enc = read_header(input_csv)
    markers, fmt = detect_markers(header)

    if not markers:
        raise ValueError(
            f"Could not identify any marker channels in the header of '{input_csv}'. "
            "This usually means the file is neither a recognized HALO nor Horizon "
            "export. Please create metadata_markers.csv by hand instead (see the "
            "template shipped with the pipeline)."
        )

    if fmt == "HORIZON":
        print(
            "[COMET] Detected Horizon export. The pipeline supports Horizon "
            "format natively. Localizations for markers not in the known-marker "
            "lookup will be defaulted to Nucleus (Horizon reports nuclei "
            "intensities). Please review the generated metadata_markers.csv "
            "before running the pipeline."
        )

    rows = []
    unknown = []
    for marker in sorted(markers, key=str.lower):
        mtype, localization = classify(marker, markers[marker])
        if marker.upper() not in KNOWN_MARKERS and mtype != "Other":
            unknown.append(marker)
        rows.append((marker, mtype, localization))

    with open(output_csv, "w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh)
        writer.writerow([
            "### Metadata for the snakemake automatic pipeline. Insert channel "
            "(exact name of protein/transcript)", "", ""
        ])
        writer.writerow([
            '### type ("Protein" or "Transcript")   & localization '
            '(For Proteins: "Cytoplasm" or "Nucleus")', "", ""
        ])
        writer.writerow(["channel", "type", "localization"])
        for marker, mtype, localization in rows:
            writer.writerow([marker, mtype, localization])

    print(f"Detected input format: {fmt} (read using {enc} encoding)")
    print(f"Wrote {len(rows)} marker rows to {output_csv}")
    if unknown:
        print(
            "These markers were not in the known-marker lookup and were defaulted "
            f"to Protein/{DEFAULT_PROTEIN_LOCALIZATION} -- please double-check them "
            f"(and edit {output_csv} if needed) before running the pipeline:"
        )
        for m in unknown:
            print(f"  - {m}")

    return rows, fmt


# Alias expected by the Snakefile's parse-time import
generate_from_csv = generate


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--input", required=True, help="Path to a HALO/Horizon export CSV")
    p.add_argument("--output", required=True, help="Path to write metadata_markers.csv")
    args = p.parse_args()
    generate(args.input, args.output)


if __name__ == "__main__":
    main()
