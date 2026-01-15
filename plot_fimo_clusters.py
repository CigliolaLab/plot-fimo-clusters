#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Plots clustered FIMO motif hits relative to TSS.

Author: Nicolas Noel
Repository: https://github.com/CigliolaLab/plot-fimo-clusters

Outputs:
  - <outprefix>_windows.csv
  - <outprefix>_clusters.csv
  - <outprefix>_tracks.png/.pdf
  - <outprefix>_heatmap.png/.pdf
  - <outprefix>_cluster<N>_hits.tsv  (one per cluster, all motif hits inside)

Typical usage (Gata4, danRer11, minus strand, relative coords):

  python3 plot_fimo_clusters.py \
    --fimo fimo_out_1e-3.tsv \
    --outprefix new_gata4_clusters_1e-3 \
    --title "Gata4 promoter: cluster detection" \
    --relative \
    --window 250 --step 50 \
    --min_hits 4 --min_fams 3 --q_thresh 0.05 \
    --bins 60
"""

import argparse
import re
from pathlib import Path
from typing import List, Tuple

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# ----------------------- Family categorization -----------------------

FAMILY_RULES = [
    ("SOX",        r"\bSOX"),
    ("GATA",       r"\bGATA"),
    ("Homeobox/TALE", r"(HOX[A-D]?|MEIS|PBX|PITX|PAX)"),
    ("STAT",       r"\bSTAT"),
    ("AP-1",       r"\b(JUN|FOS|AP-1)"),
    ("KLF",        r"\bKLF"),
    ("FOX",        r"\bFOX(A|O)?\b"),
    ("SP",         r"\bSP[1-9]\b"),
    ("TEAD",       r"\bTEAD"),
    ("ETS",        r"\b(ETS|ETV)\b"),
    ("RUNX",       r"\bRUNX"),
    ("CEBP",       r"\bCEBP\b"),
    ("HIF",        r"\b(HIF1A|HIF2A)\b"),
]
FAMILY_ORDER = [
    "SOX","GATA","STAT","AP-1","KLF","SP","TEAD",
    "RUNX","CEBP","HIF","FOX","Homeobox/TALE","Other",
]

def assign_family(name: str) -> str:
    n = ("" if pd.isna(name) else str(name)).upper()
    for fam, pat in FAMILY_RULES:
        if re.search(pat, n, flags=re.IGNORECASE):
            return fam
    return "Other"

# ----------------------- Utilities -----------------------

def strand_sign(strand: str) -> int:
    return 1 if str(strand).strip() in ["+","plus","pos","positive"] else -1

def to_midpoints(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    df["midpoint"] = (df["start"].astype(float) + df["stop"].astype(float)) / 2.0
    return df

def window_iter(start: float, end: float, win: int, step: int) -> List[Tuple[float,float]]:
    if end < start:
        start, end = end, start
    edges = []
    x = start
    while x <= end:
        edges.append((x, min(x + win, end)))
        x += step
    if not edges or edges[-1][1] < end:
        edges.append((max(end - win, start), end))
    return edges

def merge_intervals(intervals: List[Tuple[float,float]]) -> List[Tuple[float,float]]:
    if not intervals:
        return []
    intervals = sorted(intervals, key=lambda t: (t[0], t[1]))
    merged = [intervals[0]]
    for s, e in intervals[1:]:
        ps, pe = merged[-1]
        if s <= pe:  # overlap
            merged[-1] = (ps, max(pe, e))
        else:
            merged.append((s, e))
    return merged

# ----------------------- Cluster calling -----------------------

def call_clusters(df: pd.DataFrame,
                  win: int,
                  step: int,
                  q_thresh: float,
                  min_hits: int,
                  min_fams: int,
                  xcol: str = "midpoint") -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    df columns needed: [xcol, 'q-value', 'family'].
    Returns:
      windows_df: per-window stats
      clusters_df: merged cluster intervals
    """
    x = df[xcol].values
    q = df["q-value"].astype(float).values
    fams = df["family"].astype(str).values

    xmin, xmax = np.nanmin(x), np.nanmax(x)
    wins = window_iter(xmin, xmax, win, step)

    rows = []
    tiny = np.finfo(float).tiny

    for (ws, we) in wins:
        mask = (x >= ws) & (x <= we)
        if not mask.any():
            rows.append((ws, we, 0, 0, np.nan, 0.0, False))
            continue

        hits = int(mask.sum())
        uniq_fams = int(len(set(fams[mask])))
        sig = q[mask]
        mean_neglog10q = float(np.nanmean(-np.log10(np.clip(sig, tiny, 1.0))))
        strength = float(hits * uniq_fams * (mean_neglog10q if not np.isnan(mean_neglog10q) else 0.0))
        passed = bool((hits >= min_hits) and (uniq_fams >= min_fams) and (np.nanmin(sig) <= q_thresh))
        rows.append((ws, we, hits, uniq_fams, mean_neglog10q, strength, passed))

    windows_df = pd.DataFrame(rows, columns=[
        "w_start","w_end","hits","unique_families","mean_neglog10_q","strength","passed"
    ])

    # intervals for windows that passed
    pass_intervals = windows_df.loc[windows_df["passed"], ["w_start","w_end"]].itertuples(index=False, name=None)
    merged = merge_intervals(list(pass_intervals))

    clust_rows = []
    for cs, ce in merged:
        m = (x >= cs) & (x <= ce)
        if not m.any():
            continue
        subq = q[m]
        subf = fams[m]
        clust_rows.append({
            "c_start": float(cs),
            "c_end": float(ce),
            "length": float(ce - cs),
            "hits": int(m.sum()),
            "unique_families": int(len(set(subf))),
            "min_q": float(np.nanmin(subq)),
            "mean_neglog10_q": float(np.nanmean(-np.log10(np.clip(subq, tiny, 1.0)))),
            "center": float((cs + ce) / 2.0),
        })

    clusters_df = pd.DataFrame(clust_rows, columns=[
        "c_start","c_end","length","hits","unique_families","min_q","mean_neglog10_q","center"
    ]).sort_values(["min_q","hits","unique_families","length"], ascending=[True, False, False, True])

    return windows_df, clusters_df

# ----------------------- Plotting -----------------------

def plot_tracks(df: pd.DataFrame,
                clusters_df: pd.DataFrame,
                outprefix: Path,
                xlabel: str,
                tss_x: float,
                title: str):
    plt.figure(figsize=(13, 6))

    for fam in [f for f in FAMILY_ORDER if f in df["family"].unique()]:
        sub = df[df["family"] == fam]
        plt.scatter(sub["x"], sub["family"], s=26, alpha=0.85)

    # TSS line
    plt.axvline(tss_x, linestyle="--", linewidth=1, color="black")
    ymax = len(FAMILY_ORDER) - 0.5
    plt.text(tss_x, ymax, "Gata4 TSS", rotation=90, va="bottom", ha="right", fontsize=9)

    # Shade clusters
    for r in clusters_df.itertuples():
        plt.axvspan(r.c_start, r.c_end, color="orange", alpha=0.18)

    plt.xlabel(xlabel)
    plt.ylabel("Transcription Factor Family")
    plt.title(title or "Motif families with detected clusters")
    plt.tight_layout()

    png = outprefix.with_name(outprefix.name + "_tracks.png")
    pdf = outprefix.with_name(outprefix.name + "_tracks.pdf")
    plt.savefig(png, dpi=300)
    plt.savefig(pdf)
    plt.close()

def plot_heatmap(df: pd.DataFrame,
                 outprefix: Path,
                 bins: int,
                 xlabel: str,
                 title: str):
    fams = [f for f in FAMILY_ORDER if f in df["family"].unique()]
    if not fams:
        return

    fam_to_idx = {f: i for i, f in enumerate(fams)}
    x = df["x"].values
    xmin, xmax = np.nanmin(x), np.nanmax(x)
    bin_edges = np.linspace(xmin, xmax, bins + 1)

    mat = np.zeros((len(fams), bins), dtype=float)
    bid = np.clip(np.digitize(x, bin_edges) - 1, 0, bins - 1)
    for fam, b in zip(df["family"], bid):
        if fam in fam_to_idx:
            mat[fam_to_idx[fam], b] += 1

    plt.figure(figsize=(13, 6))
    plt.imshow(
        mat,
        aspect="auto",
        interpolation="nearest",
        extent=[bin_edges[0], bin_edges[-1], -0.5, len(fams) - 0.5],
    )
    plt.yticks(ticks=np.arange(len(fams)), labels=fams)
    plt.xlabel(xlabel)
    plt.ylabel("Transcription Factor Family")
    plt.title(title or "Motif density heatmap (counts per bin)")
    cbar = plt.colorbar()
    cbar.set_label("Count per bin")
    plt.tight_layout()

    png = outprefix.with_name(outprefix.name + "_heatmap.png")
    pdf = outprefix.with_name(outprefix.name + "_heatmap.pdf")
    plt.savefig(png, dpi=300)
    plt.savefig(pdf)
    plt.close()

# ----------------------- Export motifs inside clusters -----------------------

def export_cluster_hits(df: pd.DataFrame,
                        clusters_df: pd.DataFrame,
                        outprefix: Path,
                        tss: float,
                        strand: str,
                        use_relative: bool):
    """
    For each cluster, export a TSV listing ALL motif hits that fall inside that cluster window.
    One file per cluster: <outprefix>_cluster<N>_hits.tsv

    Columns include:
      - cluster_id, cluster_rel_start, cluster_rel_end
      - cluster_abs_start, cluster_abs_end
      - sequence name, start, stop, strand
      - motif name, motif_id, score, p-value, q-value
      - midpoint (genomic), rel_to_tss_bp
    """
    if clusters_df.empty:
        print("\n[cluster hits] No clusters to export.")
        return

    sign = strand_sign(strand)
    name_col = "motif_alt_id" if "motif_alt_id" in df.columns else "motif_id"

    for idx, c in clusters_df.iterrows():
        cid = idx + 1

        if use_relative:
            # clusters and df["x"] are in same relative frame
            mask = (df["x"] >= c["c_start"]) & (df["x"] <= c["c_end"])
            if strand == "-":
                abs_start = tss - c["c_end"]
                abs_end   = tss - c["c_start"]
            else:
                abs_start = tss + c["c_start"]
                abs_end   = tss + c["c_end"]
            if abs_start > abs_end:
                abs_start, abs_end = abs_end, abs_start
            rel_start = c["c_start"]
            rel_end   = c["c_end"]
        else:
            # clusters are in genomic coords; df["midpoint"] is genomic
            mask = (df["midpoint"] >= c["c_start"]) & (df["midpoint"] <= c["c_end"])
            abs_start = float(min(c["c_start"], c["c_end"]))
            abs_end   = float(max(c["c_start"], c["c_end"]))
            # derive relative just for convenience
            rel_start = sign * (abs_start - tss)
            rel_end   = sign * (abs_end   - tss)

        sub = df.loc[mask].copy()
        if sub.empty:
            print(f"[cluster hits] cluster {cid}: no motif hits found inside.")
            continue

        # relative distance to TSS (positive upstream on minus strand)
        sub["rel_to_tss_bp"] = sign * (sub["midpoint"] - tss)

        # annotate cluster metadata
        sub.insert(0, "cluster_id", cid)
        sub["cluster_rel_start"] = rel_start
        sub["cluster_rel_end"]   = rel_end
        sub["cluster_abs_start"] = int(abs_start)
        sub["cluster_abs_end"]   = int(abs_end)

        # column ordering
        cols = [
            "cluster_id",
            "cluster_rel_start",
            "cluster_rel_end",
            "cluster_abs_start",
            "cluster_abs_end",
            "sequence name",
            "start",
            "stop",
            "strand",
            name_col,
            "motif_id",
            "score",
            "p-value",
            "q-value",
            "midpoint",
            "rel_to_tss_bp",
        ]
        cols = [c for c in cols if c in sub.columns]  # keep only existing
        sub = sub[cols].sort_values(by="q-value", ascending=True)

        out_path = outprefix.with_name(f"{outprefix.name}_cluster{cid}_hits.tsv")
        sub.to_csv(out_path, sep="\t", index=False)
        print(f"[cluster hits] wrote {out_path}  ({len(sub)} hits)")

# ----------------------- Main -----------------------

def main():
    ap = argparse.ArgumentParser(description="Detect motif clusters from FIMO TSV and plot/annotate them.")
    ap.add_argument("--fimo", required=True, help="Path to FIMO TSV")
    ap.add_argument("--outprefix", required=True, help="Output prefix (no extension)")
    ap.add_argument("--title", default=None, help="Plot title")

    # TSS / coordinate handling
    ap.add_argument("--tss", type=float, default=52982587,
                    help="TSS genomic coordinate (e.g., Gata4 NM_131236.2).")
    ap.add_argument("--strand", default="-", help="Gene strand '+' or '-' (default '-')")
    ap.add_argument("--relative", action="store_true",
                    help="If set, plot positions relative to TSS (0 = TSS). If not, use genomic coords.")

    # window/cluster thresholds
    ap.add_argument("--window", type=int, default=200, help="Sliding window size (bp)")
    ap.add_argument("--step", type=int, default=50, help="Sliding window step (bp)")
    ap.add_argument("--min_hits", type=int, default=4, help="Minimum motif hits per window to pass")
    ap.add_argument("--min_fams", type=int, default=3, help="Minimum unique TF families per window to pass")
    ap.add_argument("--q_thresh", type=float, default=0.01,
                    help="At least one motif in window must have q <= this to pass.")

    # heatmap bins
    ap.add_argument("--bins", type=int, default=60, help="Number of x bins for heatmap")

    args = ap.parse_args()

    # load FIMO
    df = pd.read_csv(args.fimo, sep="\t")

    # sanity checks
    for col in ("start", "stop", "q-value"):
        if col not in df.columns:
            raise SystemExit(f"[ERROR] FIMO TSV must contain '{col}' column.")
    name_col = "motif_alt_id" if "motif_alt_id" in df.columns else "motif_id"
    df["family"] = pd.Categorical(
        [assign_family(v) for v in df[name_col]],
        categories=FAMILY_ORDER,
        ordered=True
    )

    df = to_midpoints(df)

    # coordinate system for plotting / cluster-calling
    if args.relative:
        sgn = strand_sign(args.strand)
        df["x"] = (df["midpoint"] - float(args.tss)) * sgn
        xlabel = "Position relative to TSS (bp)"
        tss_x = 0.0
    else:
        df["x"] = df["midpoint"]
        xlabel = "Genomic coordinate (bp)"
        tss_x = float(args.tss)

    # cluster calling operates in whatever frame 'x' is
    df_for_clust = df[["x", "q-value", "family"]].rename(columns={"x": "midpoint"})
    windows_df, clusters_df = call_clusters(
        df_for_clust,
        win=args.window,
        step=args.step,
        q_thresh=args.q_thresh,
        min_hits=args.min_hits,
        min_fams=args.min_fams,
        xcol="midpoint",
    )

    outprefix = Path(args.outprefix)

    # export basic cluster/window CSVs
    clusters_csv = outprefix.with_name(outprefix.name + "_clusters.csv")
    windows_csv  = outprefix.with_name(outprefix.name + "_windows.csv")
    clusters_df.to_csv(clusters_csv, index=False)
    windows_df.to_csv(windows_csv, index=False)
    print(f"[windows]  {windows_csv}")
    print(f"[clusters] {clusters_csv}")

    # export motif-level hits per cluster (ALL motifs in each window, no q filter)
    export_cluster_hits(df, clusters_df, outprefix, args.tss, args.strand, args.relative)

    # plots
    plot_tracks(
        df,
        clusters_df,
        outprefix,
        xlabel,
        tss_x,
        title=args.title or "Motif families with detected clusters",
    )
    plot_heatmap(
        df,
        outprefix,
        bins=args.bins,
        xlabel=xlabel,
        title="Motif density heatmap (counts per bin)",
    )

    print(f"[DONE] Wrote plots to {outprefix}_tracks.(png,pdf) and {outprefix}_heatmap.(png,pdf)")

if __name__ == "__main__":
    main()
