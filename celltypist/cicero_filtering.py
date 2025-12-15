#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import os
import sys
from glob import glob
from typing import Set

import numpy as np
import pandas as pd
from tqdm import tqdm
from scipy import io
from scipy.sparse import csr_matrix
import shutil


def get_home() -> str:
    home = os.environ.get("HOME")
    if not home:
        raise RuntimeError("Environment variable HOME is not set")
    return home

HOME = get_home()


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description=(
            "Filter Cicero peaks using CellTypist core regions and coaccessibility threshold, "
            "and prepare a reduced matrix for scOpen input."
        )
    )
    p.add_argument("--data_name", required=True, help="Dataset name, e.g. PBMC10k")
    p.add_argument("--window", type=int, required=True, help="Window size in bp (must match core regions)")
    p.add_argument(
        "--coacc_thresh",
        type=float,
        required=True,
        help="Cicero coaccessibility threshold (e.g. 0.1)",
    )
    p.add_argument(
        "--quiet",
        action="store_true",
        help="Reduce logging",
    )
    return p.parse_args()


def log(msg: str, quiet: bool = False) -> None:
    if not quiet:
        print(msg, flush=True)


def main() -> int:
    args = parse_args()
    data_name = args.data_name
    window = args.window
    coacc_thresh = args.coacc_thresh
    quiet = args.quiet

    base = os.path.join(HOME, "Datasets", data_name)

    cicero_filtered_dir = os.path.join(base, "output", "cicero", "filtered")
    cicero_dir = os.path.join(base, "output", "cicero")
    out_dir = os.path.join(base, "output", "celltypist_cicero")

    core_file = os.path.join(
        base,
        "output",
        "celltypist",
        f"core_genes_regions_{window}bp_intersected_peaks.txt",
    )

    out_file = os.path.join(
        out_dir,
        f"selected_peaks_{window}bp_coacc{coacc_thresh}.txt",
    )

    # ---------- checks ----------
    for d in [cicero_filtered_dir, cicero_dir]:
        if not os.path.isdir(d):
            print(f"[ERROR] Required directory not found: {d}", file=sys.stderr)
            return 2

    if not os.path.isfile(core_file):
        print(f"[ERROR] Core peaks file not found: {core_file}", file=sys.stderr)
        return 3

    os.makedirs(out_dir, exist_ok=True)

    # ---------- load core peaks ----------
    with open(core_file) as f:
        core_peak_ids: Set[str] = {
            line.strip().strip('"') for line in f if line.strip()
        }

    log(f"[INFO] Core peaks: {len(core_peak_ids)}", quiet)

    # ---------- iterate Cicero coaccessibility files ----------
    files = sorted(glob(os.path.join(cicero_filtered_dir, "peaks_chr*.csv")))
    if not files:
        print(
            f"[ERROR] No peaks_chr*.csv found in {cicero_filtered_dir}",
            file=sys.stderr,
        )
        return 4

    log(f"[INFO] Found {len(files)} Cicero files", quiet)

    selected_peak_ids = set(core_peak_ids)

    total_connections = 0
    core_connections = 0

    for file in tqdm(files, disable=quiet):
        df = pd.read_csv(file)
        total_connections += df.shape[0]

        # normalize peak IDs
        df["Peak1"] = df["Peak1"].astype(str).str.replace('"', '')
        df["Peak2"] = df["Peak2"].astype(str).str.replace('"', '')

        # threshold by coaccessibility
        df_thr = df[df["coaccess"] >= coacc_thresh]

        # keep connections touching the core
        mask = df_thr["Peak1"].isin(core_peak_ids) | df_thr["Peak2"].isin(core_peak_ids)
        df_core = df_thr[mask]

        core_connections += df_core.shape[0]

        selected_peak_ids.update(df_core["Peak1"].tolist())
        selected_peak_ids.update(df_core["Peak2"].tolist())

    log(f"[INFO] Total connections: {total_connections}", quiet)
    log(f"[INFO] Core-related connections: {core_connections}", quiet)
    log(
        f"[INFO] Selected peaks (core + coaccessible): {len(selected_peak_ids)}",
        quiet,
    )

    # ---------- save selected peak IDs ----------
    with open(out_file, "w") as f:
        for pid in sorted(selected_peak_ids):
            f.write(pid + "\n")

    log(f"[INFO] Selected peak list saved to: {out_file}", quiet)

    # ---------- load original matrix ----------
    peaks_ids_path = os.path.join(cicero_dir, "peaks.txt")
    if not os.path.isfile(peaks_ids_path):
        print(f"[ERROR] peaks.txt not found: {peaks_ids_path}", file=sys.stderr)
        return 5

    with open(peaks_ids_path) as f:
        all_peak_ids = [line.strip() for line in f if line.strip()]

    log(f"[INFO] Peaks in original matrix: {len(all_peak_ids)}", quiet)

    idx = [i for i, pid in enumerate(all_peak_ids) if pid in selected_peak_ids]
    if not idx:
        raise RuntimeError(
            "No selected peaks found in peaks.txt — check peak ID format."
        )

    idx = np.array(idx, dtype=int)
    log(f"[INFO] Peaks retained: {len(idx)}", quiet)

    mat_path = os.path.join(cicero_dir, "matrix.mtx")
    if not os.path.isfile(mat_path):
        print(f"[ERROR] matrix.mtx not found: {mat_path}", file=sys.stderr)
        return 6

    mat = io.mmread(mat_path).tocsr()
    log(f"[INFO] Original matrix shape: {mat.shape}", quiet)

    mat_core = mat[idx, :]
    log(f"[INFO] Reduced matrix shape: {mat_core.shape}", quiet)

    # ---------- write scOpen input ----------
    core_peak_ids_ordered = [all_peak_ids[i] for i in idx]

    with open(os.path.join(out_dir, "peaks.txt"), "w") as f:
        for pid in core_peak_ids_ordered:
            f.write(pid + "\n")

    with open(os.path.join(out_dir, "peaks.bed"), "w") as f:
        for pid in core_peak_ids_ordered:
            f.write(pid.replace("_", "\t") + "\n")

    shutil.copy(
        os.path.join(cicero_dir, "barcodes.tsv"),
        os.path.join(out_dir, "barcodes.tsv"),
    )

    io.mmwrite(os.path.join(out_dir, "matrix.mtx"), mat_core)

    log(f"[DONE] scOpen input written to: {out_dir}", quiet)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
