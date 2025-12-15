#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import os
import sys
from typing import Optional

import numpy as np
import pandas as pd
from tqdm import tqdm
from celltypist import models


def get_home() -> str:
    home = os.environ.get("HOME")
    if not home:
        raise RuntimeError("Environment variable HOME is not set")
    return home

HOME = get_home()


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description=(
            "Build CellTypist core gene set (top-N genes per class by |coef|) "
            "and generate BED windows around TSS±WINDOW"
        )
    )
    p.add_argument("--data_name", required=True, help="Dataset name, e.g. PBMC10k")
    p.add_argument("--model_name", required=True, help="CellTypist model name, e.g. Healthy_COVID19_PBMC")
    p.add_argument("--top_n_genes_per_class", type=int, required=True, help="Top-N genes per class")
    p.add_argument("--window", type=int, required=True, help="TSS window size in bp (±WINDOW)")

    # Optional overrides (defaults aligned with your shell pipeline style: explicit $HOME)
    p.add_argument(
        "--gtf",
        default=os.path.join(HOME, "Datasets", "gencode.v49.basic.annotation.gtf"),
        help="Path to GTF file (default: $HOME/Datasets/gencode.v49.basic.annotation.gtf)",
    )
    p.add_argument(
        "--model_dir",
        default=os.path.join(HOME, ".celltypist", "data", "models"),
        help="Directory with CellTypist .pkl models (default: $HOME/.celltypist/data/models)",
    )
    p.add_argument(
        "--out_dir",
        default=None,
        help=(
            "Output directory (default: $HOME/Datasets/{DATA_NAME}/output/celltypist). "
            "If provided, overrides the default."
        ),
    )
    p.add_argument(
        "--quiet",
        action="store_true",
        help="Reduce logging.",
    )
    return p.parse_args()


def log(msg: str, quiet: bool = False) -> None:
    if not quiet:
        print(msg, flush=True)


def parse_attr(attrs: str, key: str) -> Optional[str]:
    """
    Parse GTF attributes field and return value for a given key.
    Example fragment: gene_name "ABCD"; gene_id "ENSG...";
    """
    if not isinstance(attrs, str):
        return None
    for part in attrs.split(";"):
        part = part.strip()
        if not part:
            continue
        # GTF style: key "value"
        if part.startswith(key + " "):
            val = part.split(" ", 1)[1].strip().strip('"')
            return val
    return None


def main() -> int:
    args = parse_args()

    data_name = args.data_name
    model_name = args.model_name
    top_n = args.top_n_genes_per_class
    window = args.window
    quiet = args.quiet

    # Resolve output paths (replicate your structure; avoid "~" for robustness)
    if args.out_dir is None:
        out_dir = os.path.join(HOME, "Datasets", data_name, "output", "celltypist")
    else:
        out_dir = os.path.expanduser(args.out_dir)  
    os.makedirs(out_dir, exist_ok=True)

    out_genes_list = os.path.join(out_dir, "celltypist_core_genes.txt")
    out_table = os.path.join(out_dir, "celltypist_core_genes_per_class.csv")
    out_bed = os.path.join(out_dir, f"core_genes_regions_{window}bp.bed")

    # Load CellTypist model
    model_path = os.path.join(os.path.expanduser(args.model_dir), model_name + ".pkl")
    if not os.path.isfile(model_path):
        print(f"[ERROR] Model not found: {model_path}", file=sys.stderr)
        return 2

    log(f"[INFO] Loading model: {model_path}", quiet)
    model = models.Model.load(model_path)
    clf = model.classifier

    # Expect sklearn-like linear classifier with coef_
    if not hasattr(clf, "coef_"):
        print("[ERROR] Model classifier has no coef_. Not a linear model?", file=sys.stderr)
        return 3
    if not hasattr(clf, "features"):
        print("[ERROR] Model classifier has no 'features' attribute with gene names.", file=sys.stderr)
        return 4
    if not hasattr(clf, "classes_"):
        print("[ERROR] Model classifier has no 'classes_' attribute with cell types.", file=sys.stderr)
        return 5

    weights = clf.coef_  # (n_classes, n_genes)
    genes = np.array(clf.features)
    celltypes = np.array(clf.classes_)

    log(f"[INFO] Model classes: {len(celltypes)}, genes: {len(genes)}", quiet)

    if weights.shape[0] != len(celltypes) or weights.shape[1] != len(genes):
        print(
            f"[ERROR] Shape mismatch: coef_={weights.shape}, classes={len(celltypes)}, genes={len(genes)}",
            file=sys.stderr,
        )
        return 6

    if top_n <= 0:
        print("[ERROR] top_n_genes_per_class must be > 0", file=sys.stderr)
        return 7
    if top_n > len(genes):
        log(f"[WARN] top_n_genes_per_class ({top_n}) > #genes ({len(genes)}). Will cap to #genes.", quiet)
        top_n = len(genes)

    # Build core genes and per-class table
    top_genes_per_class = []
    rows = []

    log("[INFO] Selecting top genes per class by |coef|...", quiet)
    for i, ct in tqdm(list(enumerate(celltypes)), disable=quiet):
        class_weights = weights[i, :]
        idx_sorted = np.argsort(np.abs(class_weights))[::-1]
        idx_top = idx_sorted[:top_n]

        genes_top = genes[idx_top]
        weights_top = class_weights[idx_top]

        top_genes_per_class.append(set(genes_top.tolist()))

        for g, w in zip(genes_top, weights_top):
            rows.append(
                {
                    "cell_type": str(ct),
                    "gene": str(g),
                    "weight": float(w),
                    "abs_weight": float(abs(w)),
                }
            )

    core_genes = sorted(set().union(*top_genes_per_class))
    log(f"[INFO] Core genes (unique union): {len(core_genes)}", quiet)

    # Save core gene list
    with open(out_genes_list, "w") as f:
        for g in core_genes:
            f.write(g + "\n")
    log(f"[INFO] Saved core gene list: {out_genes_list}", quiet)

    # Save detailed table
    df = pd.DataFrame(rows)
    df.sort_values(["cell_type", "abs_weight"], ascending=[True, False], inplace=True)
    df.to_csv(out_table, index=False)
    log(f"[INFO] Saved per-class table: {out_table}", quiet)

    # Load GTF and build BED around TSS
    gtf_file = os.path.expanduser(args.gtf)  
    if not os.path.isfile(gtf_file):
        print(f"[ERROR] GTF not found: {gtf_file}", file=sys.stderr)
        return 8

    log(f"[INFO] Loading GTF: {gtf_file}", quiet)
    gtf_cols = ["chr", "source", "feature", "start", "end", "score", "strand", "frame", "attributes"]
    gtf = pd.read_csv(
        gtf_file,
        sep="\t",
        comment="#",
        header=None,
        names=gtf_cols,
        dtype={
            "chr": str,
            "source": str,
            "feature": str,
            "start": int,
            "end": int,
            "strand": str,
            "attributes": str,
        },
        low_memory=False,
    )

    gtf_genes = gtf[gtf["feature"] == "gene"].copy()
    log(f"[INFO] GTF gene features: {len(gtf_genes)}", quiet)

    # Extract gene_name
    log("[INFO] Parsing gene_name from attributes...", quiet)
    gtf_genes["gene_name"] = gtf_genes["attributes"].apply(lambda x: parse_attr(x, "gene_name"))

    # Keep only core genes
    gtf_core = gtf_genes[gtf_genes["gene_name"].isin(core_genes)].copy()
    found = gtf_core["gene_name"].nunique(dropna=True)
    log(f"[INFO] Core genes found in GTF: {found} / {len(core_genes)}", quiet)
    if found == 0:
        print("[ERROR] None of the core genes were found in the provided GTF (gene_name match).", file=sys.stderr)
        return 9

    # Compute TSS (GTF is 1-based inclusive)
    def tss_row(row) -> int:
        return int(row["start"]) if row["strand"] == "+" else int(row["end"])

    gtf_core["tss"] = gtf_core.apply(tss_row, axis=1)

    bed_rows = []
    log(f"[INFO] Building BED windows: TSS±{window} bp...", quiet)
    for _, r in tqdm(gtf_core.iterrows(), total=len(gtf_core), disable=quiet):
        tss = int(r["tss"])
        # Convert to BED 0-based half-open:
        # TSS in GTF is 1-based position; BED start is 0-based.
        start = max(0, tss - 1 - window)
        end = tss + window  # == tss + window in 1-based coordinates
        bed_rows.append(
            {
                "chr": r["chr"],
                "start": start,
                "end": end,
                "name": r["gene_name"],
                "strand": r["strand"],
            }
        )

    bed = pd.DataFrame(bed_rows)
    bed.dropna(subset=["chr", "start", "end", "name"], inplace=True)
    bed.sort_values(["chr", "start", "end", "name"], inplace=True)

    bed.to_csv(out_bed, sep="\t", header=False, index=False)
    log(f"[INFO] Saved BED: {out_bed}", quiet)

    log("[DONE] Core genes + regions generated successfully.", quiet)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
