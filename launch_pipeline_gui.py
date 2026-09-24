#!/usr/bin/env python3
"""
launch_pipeline_gui.py

Windows GUI launcher for the COMET Spatial Pipeline.
Double-click RUN_COMET_PIPELINE.BAT to open this window.

This file lives in the root of the pipeline git repository,
alongside docker-compose.yml and the Snakefile.
"""

import os
import re
import sys
import time
import shutil
import platform
import subprocess
import threading
import webbrowser

try:
    import tkinter as tk
    from tkinter import ttk, filedialog, messagebox, scrolledtext
except ImportError:
    sys.exit(
        "ERROR: tkinter is not available.\n"
        "On Windows, reinstall Python from https://www.python.org/ and make sure\n"
        "'tcl/tk and IDLE' is checked during installation."
    )

# ---------------------------------------------------------------------------
# Path resolution
# ---------------------------------------------------------------------------
SCRIPT_DIR   = os.path.dirname(os.path.abspath(__file__))
COMPOSE_FILE = os.path.join(SCRIPT_DIR, "docker-compose.yml")

_DOCKER_DESKTOP_CANDIDATES = [
    os.path.join(
        os.environ.get("ProgramFiles", r"C:\Program Files"),
        "Docker", "Docker", "Docker Desktop.exe",
    ),
    os.path.join(
        os.environ.get("LOCALAPPDATA", ""),
        "Programs", "Docker", "Docker", "Docker Desktop.exe",
    ),
]

DOCKER_READY_TIMEOUT = 120

# ---------------------------------------------------------------------------
# Theme palette  (Option B — Clean Light + Blue Accent)
# ---------------------------------------------------------------------------
T_BG     = "#f5f7fa"
T_FG     = "#1e2937"
T_ACCENT = "#1a56db"
T_BORDER = "#c5d3e8"
T_LF_LBL = "#2c3e6b"
T_HINT   = "#6b7a99"
T_RED    = "#d93025"
T_TEAL   = "#0d9488"
T_HDR    = "#1a3a5c"


# ---------------------------------------------------------------------------
# UNC → drive-letter conversion for Docker volume mounts
# ---------------------------------------------------------------------------

def _find_drive_for_share(net_use_out: str, unc_share: str) -> "str | None":
    """
    Scan 'net use' output for a line that maps a drive letter to unc_share.
    Returns the drive letter (e.g. 'Z') or None.
    """
    target = unc_share.rstrip("\\").lower()
    for line in net_use_out.splitlines():
        cols = line.split()
        drives = [c for c in cols if len(c) == 2 and c[1] == ":" and c[0].isalpha()]
        uncs   = [c for c in cols if c.startswith("\\\\")]
        if drives and uncs and uncs[0].rstrip("\\").lower() == target:
            return drives[0][0].upper()
    return None


def _free_drive_letter(net_use_out: str) -> str:
    """Return a drive letter not currently in use (checked via net use + fsutil)."""
    used: set[str] = set()
    for line in net_use_out.splitlines():
        for col in line.split():
            if len(col) == 2 and col[1] == ":" and col[0].isalpha():
                used.add(col[0].upper())
    r = subprocess.run(["fsutil", "fsinfo", "drives"],
                       capture_output=True, text=True)
    for token in r.stdout.split():
        if len(token) >= 2 and token[1] == ":" and token[0].isalpha():
            used.add(token[0].upper())
    for letter in "ZYXWVUTSRQPONMLKJIHGFE":
        if letter not in used:
            return letter
    raise RuntimeError("No free drive letters available to map the network share.")


def _is_network_drive(drive_letter: str) -> bool:
    """Return True if a drive letter (e.g. 'Z:') is a mapped network drive."""
    r = subprocess.run(
        ["net", "use", f"{drive_letter}:"],
        capture_output=True, text=True,
    )
    return r.returncode == 0


def _resolve_docker_path(path: str) -> str:
    """
    Docker Desktop on Windows cannot mount UNC paths (\\\\server\\share\\...).
    If path is a UNC path, this function maps the share to a drive letter
    (reusing an existing mapping where possible) and returns the equivalent
    drive-letter path.  Non-UNC paths are returned unchanged.
    """
    if not path.startswith("\\\\"):
        return path

    # Parse  \\server\share\optional\sub\path
    bare  = path[2:]
    parts = bare.split("\\", 2)
    if len(parts) < 2:
        raise RuntimeError(f"Cannot parse UNC path: {path}")

    unc_share = f"\\\\{parts[0]}\\{parts[1]}"
    remainder = parts[2] if len(parts) > 2 else ""

    # Check existing mappings
    net_out = subprocess.run(
        ["net", "use"], capture_output=True, text=True
    ).stdout

    drive = _find_drive_for_share(net_out, unc_share)

    if drive is None:
        # Attempt to map the share to a free drive letter
        drive = _free_drive_letter(net_out)
        r = subprocess.run(
            ["net", "use", f"{drive}:", unc_share],
            capture_output=True, text=True,
        )
        if r.returncode != 0:
            err = (r.stderr or r.stdout).strip()
            raise RuntimeError(
                f"Could not map {unc_share} to {drive}: — {err}\n\n"
                "Please map the network share to a drive letter manually in\n"
                "File Explorer (right-click → Map network drive) and then\n"
                "select the mapped path in the Output Folder field."
            )

    mapped = f"{drive}:\\"
    if remainder:
        mapped = os.path.join(mapped, remainder)
    return mapped


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _safe_project_name(csv_path: str) -> str:
    stem = os.path.splitext(os.path.basename(csv_path))[0]
    stem = re.sub(r"[^A-Za-z0-9_\-]", "_", stem)
    stem = re.sub(r"_+", "_", stem).strip("_")
    return stem or "project"


def _write_config(
    run_dir: str,
    project_name: str,
    csv_filename: str,
    pca_dims: int,
    resolution: float,
    port: int,
) -> str:
    output_dir = f"{project_name}_results"
    content = (
        "# COMET Spatial Pipeline Configuration\n"
        "# Auto-generated by launch_pipeline_gui.py\n\n"
        f'name_of_project: "{project_name}"\n\n'
        f'HALO_data_file: "{csv_filename}"\n'
        'metadata_file: "metadata_markers.csv"\n\n'
        f'output_directory: "{output_dir}"\n\n'
        f"pca_dims: {pca_dims}\n"
        f"clustering_resolution: {resolution}\n"
        f"streamlit_port: {port}\n"
    )
    config_path = os.path.join(run_dir, "config.yaml")
    with open(config_path, "w", encoding="utf-8") as fh:
        fh.write(content)
    return config_path


def _popen_kwargs() -> dict:
    kw: dict = {}
    if platform.system() == "Windows":
        si = subprocess.STARTUPINFO()
        si.dwFlags |= subprocess.STARTF_USESHOWWINDOW
        si.wShowWindow = subprocess.SW_HIDE
        kw["startupinfo"] = si
    return kw


def _docker_is_running() -> bool:
    try:
        r = subprocess.run(
            ["docker", "info"], capture_output=True, timeout=6, **_popen_kwargs()
        )
        return r.returncode == 0
    except Exception:
        return False


def _find_docker_desktop() -> "str | None":
    for path in _DOCKER_DESKTOP_CANDIDATES:
        if path and os.path.isfile(path):
            return path
    return None


# ---------------------------------------------------------------------------
# Main GUI
# ---------------------------------------------------------------------------

class CometLauncherApp(tk.Tk):
    def __init__(self) -> None:
        super().__init__()
        self.title("COMET Pipeline Launcher")
        self.resizable(True, True)
        self.minsize(860, 760)

        self._process: "subprocess.Popen | None" = None
        self._running = False
        self._browser_opened = False
        self._active_port = 8501

        self._csv_var  = tk.StringVar()
        self._name_var = tk.StringVar()
        self._out_var  = tk.StringVar()

        self._apply_theme()
        self._build_ui()
        self._center_window()

    # ── Theme ────────────────────────────────────────────────────────────────

    def _apply_theme(self) -> None:
        self.configure(bg=T_BG)
        s = ttk.Style(self)
        s.theme_use("clam")

        s.configure("TFrame",   background=T_BG)
        s.configure("TLabel",   background=T_BG, foreground=T_FG,
                    font=("Helvetica", 11))
        s.configure("TEntry",   fieldbackground="white",
                    bordercolor=T_BORDER, lightcolor=T_BORDER, darkcolor=T_BORDER,
                    font=("Helvetica", 11), padding=5)
        s.configure("TSpinbox", fieldbackground="white",
                    bordercolor=T_BORDER, lightcolor=T_BORDER, darkcolor=T_BORDER,
                    font=("Helvetica", 11), padding=5)
        s.map("TEntry",   bordercolor=[("focus", T_ACCENT)])
        s.map("TSpinbox", bordercolor=[("focus", T_ACCENT)])

        s.configure("TLabelframe",
                    background=T_BG, bordercolor=T_BORDER,
                    lightcolor=T_BORDER, darkcolor=T_BORDER,
                    relief="solid", borderwidth=1)
        s.configure("TLabelframe.Label",
                    background=T_BG, foreground=T_LF_LBL,
                    font=("Helvetica", 10, "bold"), padding=(4, 0))

        _btn = dict(borderwidth=0, focusthickness=2,
                    focuscolor=T_BORDER, padding=(16, 8))

        s.configure("Accent.TButton", background=T_ACCENT, foreground="white",
                    font=("Helvetica", 11, "bold"), **_btn)
        s.map("Accent.TButton",
              background=[("active", "#1648c4"), ("disabled", "#a8bfe8")],
              foreground=[("disabled", "#dde5f5")])

        s.configure("Danger.TButton", background=T_RED, foreground="white",
                    font=("Helvetica", 11), **_btn)
        s.map("Danger.TButton",
              background=[("active", "#b52a20"), ("disabled", "#d0d0d0")],
              foreground=[("disabled", "#a0a0a0")])

        s.configure("Teal.TButton", background=T_TEAL, foreground="white",
                    font=("Helvetica", 11), **_btn)
        s.map("Teal.TButton",
              background=[("active", "#0a7a70"), ("disabled", "#d0d0d0")],
              foreground=[("disabled", "#a0a0a0")])

        s.configure("Muted.TButton", background="#e2e8f0", foreground=T_FG,
                    font=("Helvetica", 11), **_btn)
        s.map("Muted.TButton",
              background=[("active", "#cbd5e1"), ("disabled", "#e9ecef")])

    # ── Layout ───────────────────────────────────────────────────────────────

    def _center_window(self) -> None:
        self.update_idletasks()
        w = max(self.winfo_reqwidth(), 860)
        h = max(self.winfo_reqheight(), 760)
        x = (self.winfo_screenwidth()  - w) // 2
        y = (self.winfo_screenheight() - h) // 2
        self.geometry(f"{w}x{h}+{x}+{y}")

    def _build_ui(self) -> None:
        rp = dict(padx=16, pady=7)

        # ── Header ──────────────────────────────────────────────────────────
        hdr = tk.Frame(self, bg=T_HDR)
        hdr.pack(fill="x")

        tk.Label(
            hdr, text="  COMET  Spatial Pipeline Launcher",
            font=("Helvetica", 16, "bold"), fg="white", bg=T_HDR,
        ).pack(side="left", padx=16, pady=12)

        self._docker_lbl = tk.Label(
            hdr, text="⬤  Docker: checking…",
            font=("Helvetica", 10), fg="#8aafd4", bg=T_HDR,
        )
        self._docker_lbl.pack(side="right", padx=16)

        tk.Frame(self, bg=T_ACCENT, height=3).pack(fill="x")

        # ── Main area ────────────────────────────────────────────────────────
        main = ttk.Frame(self, padding=(16, 12))
        main.pack(fill="both", expand=True)

        # ── Input CSV file ───────────────────────────────────────────────────
        file_lf = ttk.LabelFrame(main, text="Input Data File", padding=(12, 8))
        file_lf.pack(fill="x", **rp)

        ttk.Entry(file_lf, textvariable=self._csv_var).grid(
            row=0, column=0, sticky="ew", padx=(0, 10), ipady=2)
        ttk.Button(file_lf, text="Browse…",
                   command=self._browse_csv, style="Muted.TButton").grid(
            row=0, column=1)
        file_lf.columnconfigure(0, weight=1)
        ttk.Label(
            file_lf,
            text="Select the HALO or Horizon per-cell CSV export for this run.",
            foreground=T_HINT, font=("Helvetica", 10),
        ).grid(row=1, column=0, columnspan=2, sticky="w", pady=(4, 0))

        # ── Project name ─────────────────────────────────────────────────────
        proj_lf = ttk.LabelFrame(main, text="Project", padding=(12, 8))
        proj_lf.pack(fill="x", **rp)

        ttk.Label(proj_lf, text="Project name:").grid(row=0, column=0, sticky="w")
        ttk.Entry(proj_lf, textvariable=self._name_var).grid(
            row=0, column=1, sticky="ew", padx=(12, 0), ipady=2)
        ttk.Label(
            proj_lf, text="Auto-filled from the filename — edit freely.",
            foreground=T_HINT, font=("Helvetica", 10),
        ).grid(row=1, column=1, sticky="w", pady=(4, 0))
        proj_lf.columnconfigure(1, weight=1)

        # ── Output folder ────────────────────────────────────────────────────
        out_lf = ttk.LabelFrame(main, text="Output Folder", padding=(12, 8))
        out_lf.pack(fill="x", **rp)

        ttk.Entry(out_lf, textvariable=self._out_var).grid(
            row=0, column=0, sticky="ew", padx=(0, 10), ipady=2)
        ttk.Button(out_lf, text="Browse…",
                   command=self._browse_out, style="Muted.TButton").grid(
            row=0, column=1)
        out_lf.columnconfigure(0, weight=1)
        ttk.Label(
            out_lf,
            text="Results folder will be created here. Defaults to the same folder as the CSV.",
            foreground=T_HINT, font=("Helvetica", 10),
        ).grid(row=1, column=0, columnspan=2, sticky="w", pady=(4, 0))

        # ── Pipeline parameters ──────────────────────────────────────────────
        param_lf = ttk.LabelFrame(main, text="Pipeline Parameters", padding=(12, 8))
        param_lf.pack(fill="x", **rp)

        ttk.Label(param_lf, text="PCA dimensions:").grid(
            row=0, column=0, sticky="w")
        self._pca_var = tk.IntVar(value=8)
        ttk.Spinbox(param_lf, from_=1, to=50,
                    textvariable=self._pca_var, width=8).grid(
            row=0, column=1, sticky="w", padx=(12, 18))
        ttk.Label(param_lf, text="5–15 is typical for 20–40 markers",
                  foreground=T_HINT, font=("Helvetica", 10)).grid(
            row=0, column=2, sticky="w")

        ttk.Label(param_lf, text="Clustering resolution:").grid(
            row=1, column=0, sticky="w", pady=(10, 0))
        self._res_var = tk.StringVar(value="0.5")
        ttk.Entry(param_lf, textvariable=self._res_var, width=8).grid(
            row=1, column=1, sticky="w", padx=(12, 18), pady=(10, 0), ipady=2)
        ttk.Label(param_lf, text="Higher → more clusters.  Try 0.2 – 1.0",
                  foreground=T_HINT, font=("Helvetica", 10)).grid(
            row=1, column=2, sticky="w", pady=(10, 0))

        ttk.Label(param_lf, text="Annotation portal port:").grid(
            row=2, column=0, sticky="w", pady=(10, 0))
        self._port_var = tk.IntVar(value=8501)
        ttk.Spinbox(param_lf, from_=1024, to=65535,
                    textvariable=self._port_var, width=8).grid(
            row=2, column=1, sticky="w", padx=(12, 18), pady=(10, 0))
        ttk.Label(param_lf, text="Change if 8501 is already in use on this PC",
                  foreground=T_HINT, font=("Helvetica", 10)).grid(
            row=2, column=2, sticky="w", pady=(10, 0))

        param_lf.columnconfigure(2, weight=1)

        # ── Log ──────────────────────────────────────────────────────────────
        log_lf = ttk.LabelFrame(main, text="Pipeline Log", padding=(12, 8))
        log_lf.pack(fill="both", expand=True, **rp)

        self._log = scrolledtext.ScrolledText(
            log_lf,
            height=11,
            state="disabled",
            font=("Consolas", 12),
            wrap="word",
            bg="#1e2637",
            fg="#cdd6f4",
            insertbackground="#cdd6f4",
            cursor="arrow",
            selectbackground="#3d5a9e",
            selectforeground="#ffffff",
            relief="flat",
            borderwidth=0,
            padx=10,
            pady=8,
        )
        self._log.pack(fill="both", expand=True)

        # ── Buttons ──────────────────────────────────────────────────────────
        btn_frame = ttk.Frame(main)
        btn_frame.pack(fill="x", pady=(10, 4))

        self._launch_btn = ttk.Button(
            btn_frame, text="▶  Launch Pipeline",
            command=self._on_launch, style="Accent.TButton")
        self._launch_btn.pack(side="left", padx=(0, 10))

        self._stop_btn = ttk.Button(
            btn_frame, text="■  Stop",
            command=self._on_stop, style="Danger.TButton", state="disabled")
        self._stop_btn.pack(side="left")

        self._portal_btn = ttk.Button(
            btn_frame, text="🌐  Open Annotation Portal",
            command=self._open_portal, style="Teal.TButton", state="disabled")
        self._portal_btn.pack(side="right")

        # ── Status bar ──────────────────────────────────────────────────────
        tk.Frame(self, bg=T_BORDER, height=1).pack(fill="x")
        self._status_lbl = tk.Label(
            self, text="",
            font=("Helvetica", 10), fg=T_HINT, bg=T_BG,
            anchor="w",
        )
        self._status_lbl.pack(fill="x", padx=16, pady=(5, 8))

        self._csv_var.trace_add("write", self._on_csv_change)
        threading.Thread(target=self._check_docker_status, daemon=True).start()

    # ── Docker status ────────────────────────────────────────────────────────

    def _check_docker_status(self) -> None:
        if _docker_is_running():
            self.after(0, lambda: self._docker_lbl.configure(
                text="⬤  Docker: running", fg="#5fba7d"))
        else:
            self.after(0, lambda: self._docker_lbl.configure(
                text="⬤  Docker: not running", fg="#e06c75"))

    # ── Event handlers ───────────────────────────────────────────────────────

    def _on_csv_change(self, *_) -> None:
        path = self._csv_var.get()
        if path and os.path.isfile(path):
            self._name_var.set(_safe_project_name(path))
            if not self._out_var.get():
                self._out_var.set(os.path.dirname(os.path.abspath(path)))

    def _browse_csv(self) -> None:
        path = filedialog.askopenfilename(
            title="Select HALO / Horizon CSV export",
            filetypes=[("CSV files", "*.csv"), ("All files", "*.*")],
        )
        if path:
            self._csv_var.set(os.path.normpath(path))

    def _browse_out(self) -> None:
        folder = filedialog.askdirectory(
            title="Choose output folder for results",
            initialdir=self._out_var.get() or None,
        )
        if folder:
            self._out_var.set(os.path.normpath(folder))

    def _on_launch(self) -> None:
        csv_path = self._csv_var.get().strip()
        if not csv_path or not os.path.isfile(csv_path):
            messagebox.showerror("No input file",
                                 "Please select a valid CSV file before launching.")
            return

        project_name = self._name_var.get().strip()
        if not project_name:
            messagebox.showerror("No project name", "Please enter a project name.")
            return

        try:
            resolution = float(self._res_var.get())
            if resolution <= 0:
                raise ValueError
        except ValueError:
            messagebox.showerror(
                "Invalid resolution",
                "Clustering resolution must be a positive number (e.g. 0.5).",
            )
            return

        pca_dims = int(self._pca_var.get())
        port     = int(self._port_var.get())

        if not os.path.isfile(COMPOSE_FILE):
            messagebox.showerror(
                "docker-compose.yml not found",
                f"Expected:\n{COMPOSE_FILE}\n\n"
                "Make sure you are running from the correct pipeline folder.",
            )
            return

        # Resolve run / output folder
        csv_abs    = os.path.abspath(csv_path)
        out_folder = self._out_var.get().strip()
        if not out_folder:
            out_folder = os.path.dirname(csv_abs)

        run_dir = out_folder
        os.makedirs(run_dir, exist_ok=True)

        # ── Determine the folder Docker will actually mount as /data ─────────
        #
        # Docker Desktop on Windows (WSL2 backend) can only mount LOCAL drives.
        # Network shares — whether accessed via UNC (\\server\share) or a mapped
        # drive letter (Z:\) — are invisible inside the WSL2 VM that runs Docker,
        # so the volume mount silently appears as an empty directory.
        #
        # Strategy:
        #   1. Convert any UNC path to a drive letter (net use), as before.
        #   2. If that drive letter is a network drive, create a LOCAL staging
        #      folder under %LOCALAPPDATA%\COMET\<project> and run Docker there.
        #   3. After the pipeline finishes, copy all results back to run_dir.

        try:
            docker_run_dir = _resolve_docker_path(run_dir)
        except RuntimeError as exc:
            messagebox.showerror("Network path not usable with Docker", str(exc))
            return

        # Detect whether the resolved path lands on a network drive.
        drive = docker_run_dir[:1].upper()
        needs_staging = (
            docker_run_dir.startswith("\\\\")  # still a UNC (shouldn't happen)
            or (drive.isalpha() and _is_network_drive(drive))
        )

        staging_dir: "str | None" = None
        docker_data_dir = docker_run_dir   # the path Docker will mount

        if needs_staging:
            local_app = os.environ.get("LOCALAPPDATA") or os.path.expanduser("~")
            staging_dir = os.path.join(local_app, "COMET", f"{project_name}_{port}")
            docker_data_dir = staging_dir

        # Copy CSV to the data dir Docker will see
        data_dir_for_csv = staging_dir if staging_dir else run_dir
        csv_in_run = os.path.join(data_dir_for_csv, os.path.basename(csv_abs))
        needs_copy = os.path.abspath(csv_in_run) != csv_abs

        self._log_clear()

        if staging_dir:
            self._log_write(
                f"[COMET] Network drive detected — Docker will use a local staging folder.\n"
                f"        Results will be copied to your output folder after the pipeline finishes.\n\n"
                f"[COMET] Staging folder → {staging_dir}\n\n"
            )
            try:
                os.makedirs(staging_dir, exist_ok=True)
            except OSError as exc:
                messagebox.showerror("Could not create staging folder", str(exc))
                return
        elif docker_run_dir != run_dir:
            self._log_write(
                f"[COMET] Network path mapped for Docker:\n"
                f"        {run_dir}\n"
                f"     →  {docker_run_dir}\n\n"
            )

        if needs_copy:
            dest_label = staging_dir if staging_dir else run_dir
            self._log_write(
                f"[COMET] Copying data file…\n"
                f"        {csv_abs}\n"
                f"     →  {csv_in_run}\n"
            )
            try:
                shutil.copy2(csv_abs, csv_in_run)
                self._log_write("[COMET] Copy done.\n\n")
            except OSError as exc:
                messagebox.showerror("File copy failed", str(exc))
                return

        # Also stage this project's own metadata_markers.csv, not just the CSV.
        #
        # Bug found 2026-09-25: when staging was needed (network-drive run
        # folders, which is the common case in practice), this function only
        # ever copied the HALO/Horizon CSV into staging_dir -- never
        # metadata_markers.csv. Docker's entrypoint (docker-compose.yml) only
        # seeds /data/metadata_markers.csv from the image's own bundled copy
        # when /data doesn't already have one, so with metadata_markers.csv
        # missing from the staging folder, every Docker run silently fell
        # back to the pipeline image's generic placeholder marker list
        # instead of this project's real one. That placeholder list uses
        # different marker names (e.g. "Cytokeratin" instead of this
        # project's "PanCK") and is missing most of a real panel, which
        # silently reclassifies most marker columns as "extra metadata"
        # instead of expression data -- confirmed to be the actual cause of
        # Leiden cluster counts and UMAP layouts differing from a native
        # (non-Docker) run of the identical CSV, which was otherwise wrongly
        # suspected to be a package-version or CPU-architecture difference.
        # Mirror the CSV-staging logic above so the project's real
        # metadata_markers.csv always travels with it.
        if staging_dir:
            metadata_src = os.path.join(run_dir, "metadata_markers.csv")
            metadata_dst = os.path.join(staging_dir, "metadata_markers.csv")
            if os.path.exists(metadata_src) and os.path.abspath(metadata_src) != os.path.abspath(metadata_dst):
                self._log_write(
                    f"[COMET] Copying metadata_markers.csv…\n"
                    f"        {metadata_src}\n"
                    f"     →  {metadata_dst}\n"
                )
                try:
                    shutil.copy2(metadata_src, metadata_dst)
                    self._log_write("[COMET] Copy done.\n\n")
                except OSError as exc:
                    messagebox.showerror("File copy failed", str(exc))
                    return

        try:
            cfg = _write_config(
                data_dir_for_csv, project_name, os.path.basename(csv_abs),
                pca_dims, resolution, port,
            )
        except OSError as exc:
            messagebox.showerror("Config write failed", str(exc))
            return

        # Always write a copy of config.yaml to the user-facing run_dir too.
        if staging_dir:
            try:
                _write_config(run_dir, project_name, os.path.basename(csv_abs),
                              pca_dims, resolution, port)
            except OSError:
                pass  # non-fatal — the staging copy is what Docker uses

        self._log_write(f"[COMET] Config written → {cfg}\n")
        self._log_write(f"[COMET] Run folder     → {run_dir}\n\n")

        self._launch_btn.configure(state="disabled")
        self._stop_btn.configure(state="normal")
        self._portal_btn.configure(state="disabled")
        self._browser_opened = False
        self._running        = True
        self._active_port    = port

        threading.Thread(
            target=self._docker_thread,
            args=(run_dir, docker_data_dir, staging_dir, port),
            daemon=True,
        ).start()

    def _on_stop(self) -> None:
        if self._process and self._running:
            self._log_write("\n[COMET] Stopping…\n")
            # Ask docker compose to stop the container gracefully first. Just
            # calling self._process.terminate() kills the "docker compose up"
            # CLI wrapper on Windows via TerminateProcess, which does NOT stop
            # the actual container -- comet_pipeline (and Snakemake inside it)
            # keeps running in the background, invisible to the GUI. Running
            # "docker compose stop" actually stops the container; once it
            # exits, the "up" process's stdout closes naturally and the log
            # loop in _docker_thread finishes on its own.
            def _graceful_stop() -> None:
                try:
                    subprocess.run(
                        ["docker", "compose", "-f", COMPOSE_FILE, "stop"],
                        cwd=SCRIPT_DIR,
                        capture_output=True,
                        text=True,
                        timeout=30,
                        **_popen_kwargs(),
                    )
                except Exception as exc:
                    self._log_write(f"[COMET] Warning: 'docker compose stop' failed: {exc}\n")
                # Fallback in case the container didn't exit and the log loop
                # is still blocked reading stdout.
                try:
                    if self._process is not None:
                        self._process.terminate()
                except Exception:
                    pass

            threading.Thread(target=_graceful_stop, daemon=True).start()

    def _await_streamlit(self, port: int) -> None:
        """
        Background thread: poll http://localhost:{port} every 2 s until Streamlit
        responds with any HTTP reply, then enable the portal button and open the
        browser.  Uses only stdlib so no extra dependencies are needed.
        """
        import urllib.request
        import urllib.error

        url = f"http://localhost:{port}"
        while self._running:
            try:
                urllib.request.urlopen(url, timeout=3)
                # Got a response — Streamlit is up
                self._log_write(f"[COMET] Annotation portal is ready → {url}\n")
                self.after(0, lambda u=url: self._portal_btn.configure(state="normal"))
                webbrowser.open(url)
                return
            except Exception:
                time.sleep(2)

    def _open_portal(self) -> None:
        webbrowser.open(f"http://localhost:{self._active_port}")

    # ── Docker thread ────────────────────────────────────────────────────────

    def _docker_thread(
        self,
        run_dir: str,
        docker_data_dir: str,
        staging_dir: "str | None",
        port: int,
    ) -> None:
        # Step 1 — ensure Docker is running
        if not _docker_is_running():
            self._log_write(
                "[COMET] Docker is not running — starting Docker Desktop…\n")
            self.after(0, lambda: self._docker_lbl.configure(
                text="⬤  Docker: starting…", fg="#e5c07b"))

            exe = _find_docker_desktop()
            if not exe:
                self._log_write(
                    "[ERROR] Docker Desktop not found.\n"
                    "Please install it from https://www.docker.com/products/docker-desktop/\n"
                    "then try again.\n"
                )
                self._finish(None)
                return

            try:
                subprocess.Popen([exe], **_popen_kwargs())
            except Exception as exc:
                self._log_write(f"[ERROR] Could not launch Docker Desktop: {exc}\n")
                self._finish(None)
                return

            self._log_write("[COMET] Waiting for Docker to be ready")
            deadline = time.time() + DOCKER_READY_TIMEOUT
            ready = False
            while time.time() < deadline:
                time.sleep(3)
                self._log_write(".")
                if _docker_is_running():
                    ready = True
                    break

            if not ready:
                self._log_write(
                    f"\n[ERROR] Docker did not start within {DOCKER_READY_TIMEOUT}s.\n"
                    "Open Docker Desktop manually, wait for 'Docker Desktop is running',\n"
                    "then click Launch again.\n"
                )
                self._finish(None)
                return

            self._log_write("\n[COMET] Docker is ready.\n\n")
            self.after(0, lambda: self._docker_lbl.configure(
                text="⬤  Docker: running", fg="#5fba7d"))

        # Step 2 — run the pipeline
        self._log_write("[COMET] Starting pipeline…\n\n")
        self._set_status("⬤  Running…", "#e5c07b")
        self.title("⬤  COMET Pipeline — Running")

        env = os.environ.copy()
        env["RUN_DIR"]        = docker_data_dir   # local path Docker can actually mount
        env["STREAMLIT_PORT"] = str(port)

        # --exit-code-from makes docker compose return the container's exit code.
        # Without it, docker compose up always returns 0 even if the container fails.
        cmd = [
            "docker", "compose", "-f", COMPOSE_FILE,
            "up", "--build",
            "--exit-code-from", "comet_pipeline",
        ]

        try:
            self._process = subprocess.Popen(
                cmd,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                text=True,
                env=env,
                cwd=SCRIPT_DIR,
                bufsize=1,
                **_popen_kwargs(),
            )
        except FileNotFoundError:
            self._log_write(
                "[ERROR] 'docker' command not found.\n"
                "Please install Docker Desktop from https://www.docker.com/\n"
            )
            self._finish(None)
            return
        except Exception as exc:
            self._log_write(f"[ERROR] Could not start Docker: {exc}\n")
            self._finish(None)
            return

        # Noise patterns: Streamlit deprecation warnings that spam the log on every
        # page refresh.  Filter them out before writing to the GUI so the log stays
        # readable and the Tkinter event queue doesn't get flooded.
        _NOISE = (
            "use_container_width",
            "will be removed after",
            "use `width=",
            "For `use_container_width",
        )

        paused_re = re.compile(r"PIPELINE PAUSED", re.IGNORECASE)

        for line in self._process.stdout:
            if any(p in line for p in _NOISE):
                continue
            self._log_write(line)
            # When the pipeline prints its "PAUSED" banner, start a background
            # thread that polls the port.  This is more reliable than parsing
            # Streamlit's startup messages (which may carry ANSI codes or arrive
            # in a large buffered chunk after Docker flushes).
            if not self._browser_opened and paused_re.search(line):
                self._browser_opened = True   # prevent duplicate threads
                threading.Thread(
                    target=self._await_streamlit,
                    args=(port,),
                    daemon=True,
                ).start()

        retcode = self._process.wait()

        # Copy results from the local staging folder back to the network drive.
        # This must happen regardless of success/failure so partial results are
        # preserved and researchers can inspect them.
        if staging_dir:
            self._log_write("\n[COMET] Copying results back to your output folder…\n")
            self._set_status("⬤  Copying results to network drive…", "#e5c07b")
            try:
                copied, errors = 0, []
                for item in os.listdir(staging_dir):
                    src_path = os.path.join(staging_dir, item)
                    dst_path = os.path.join(run_dir, item)
                    try:
                        if os.path.isdir(src_path):
                            if os.path.exists(dst_path):
                                shutil.rmtree(dst_path)
                            shutil.copytree(src_path, dst_path)
                        else:
                            shutil.copy2(src_path, dst_path)
                        copied += 1
                    except Exception as item_exc:
                        errors.append(f"{item}: {item_exc}")
                if errors:
                    self._log_write(
                        f"[COMET] Warning: {len(errors)} item(s) could not be copied:\n"
                        + "\n".join(f"  {e}" for e in errors) + "\n"
                    )
                self._log_write(
                    f"[COMET] Results ({copied} item(s)) copied to:\n"
                    f"        {run_dir}\n"
                )
            except Exception as exc:
                self._log_write(f"[COMET] Warning: could not copy results: {exc}\n")

        self._finish(retcode)

    def _set_status(self, text: str, color: str) -> None:
        """Update the persistent status bar (thread-safe)."""
        self.after(0, lambda: self._status_lbl.configure(text=text, fg=color))

    def _finish(self, retcode: "int | None") -> None:
        self._running = False
        self._process = None
        threading.Thread(target=self._check_docker_status, daemon=True).start()

        def _ui() -> None:
            self._launch_btn.configure(state="normal")
            self._stop_btn.configure(state="disabled")
            if retcode == 0:
                self._log_write("\n[COMET] Pipeline finished successfully.\n")
                self._set_status("✓  Completed successfully", "#5fba7d")
                self.title("✓  COMET Pipeline — Done")
                try:
                    import winsound
                    winsound.MessageBeep(winsound.MB_OK)
                except Exception:
                    pass
                messagebox.showinfo("Done", "The COMET pipeline completed successfully!")
            elif retcode is not None:
                msg = (
                    f"\n[COMET] ✖  Pipeline failed (exit code {retcode}).\n"
                    "Check the log above for details.\n"
                )
                self._log_write(msg)
                self._set_status(
                    f"✖  Failed (exit code {retcode}) — see log above for the error",
                    T_RED,
                )
                self.title("✖  COMET Pipeline — Failed")
                try:
                    import winsound
                    winsound.MessageBeep(winsound.MB_ICONEXCLAMATION)
                except Exception:
                    pass
                messagebox.showerror(
                    "Pipeline failed",
                    f"The COMET pipeline exited with code {retcode}.\n\n"
                    "Check the Annotation Log for the full error message.",
                )
            else:
                # User stopped
                self._set_status("■  Stopped", T_HINT)
                self.title("COMET Pipeline Launcher")

        self.after(0, _ui)

    # ── Log helpers ──────────────────────────────────────────────────────────

    def _log_write(self, text: str) -> None:
        def _do():
            self._log.configure(state="normal")
            self._log.insert("end", text)
            self._log.see("end")
            self._log.configure(state="disabled")
        self.after(0, _do)

    def _log_clear(self) -> None:
        self._log.configure(state="normal")
        self._log.delete("1.0", "end")
        self._log.configure(state="disabled")


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def main() -> None:
    app = CometLauncherApp()
    app.mainloop()


if __name__ == "__main__":
    main()
