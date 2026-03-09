from datetime import datetime
from pathlib import Path
from typing import List

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from .models import GeneDensityResult, GCSkewResult, MobileElementAnalysis


def _save(fig: plt.Figure, output_dir: Path, prefix: str) -> Path:
    output_dir.mkdir(parents=True, exist_ok=True)
    ts = datetime.now().strftime("%Y%m%d_%H%M%S")
    path = output_dir / f"{prefix}_{ts}.png"
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return path


def plot_density_comparison(results: List[GeneDensityResult], output_dir: Path) -> Path:
    names = [r.genome_name for r in results]
    densities = [r.density_genes_per_kb for r in results]
    colors = [
        "skyblue" if "bacteria" in r.group else
        "salmon" if "organelle" in r.group else
        "gray"
        for r in results
    ]

    fig, ax = plt.subplots(figsize=(10, 6))
    ax.bar(names, densities, color=colors)
    ax.set_xticklabels(names, rotation=45, ha="right")
    ax.set_ylabel("Genes per kb")
    ax.set_title("Gene Density Comparison: Endosymbiosis Evidence")
    fig.tight_layout()

    return _save(fig, output_dir, "density_plot")


def plot_gc_skew_curve(result: GCSkewResult, output_dir: Path) -> Path:
    fig, ax = plt.subplots(figsize=(10, 4))
    ax.plot(result.genomic_positions, result.cumulative_skew, label="Cumulative Skew", color="purple")
    ax.axvline(
        x=result.predicted_ori_pos,
        color="red",
        linestyle="--",
        label=f"Predicted Ori ({result.predicted_ori_pos})",
    )
    ax.set_title(f"GC Skew Analysis: {result.genome_name}")
    ax.set_xlabel("Genomic Position (bp)")
    ax.set_ylabel("Cumulative Skew")
    ax.legend()
    ax.grid(True, alpha=0.3)
    fig.tight_layout()

    prefix = f"gc_skew_{result.genome_name}"
    return _save(fig, output_dir, prefix)


def plot_chromosome_map(result: MobileElementAnalysis, output_dir: Path) -> Path:
    if not result or not result.hits:
        print("No hits to plot.")
        return None

    ranges = [(hit.start, hit.end - hit.start) for hit in result.hits]

    fig, ax = plt.subplots(figsize=(15, 3))
    ax.broken_barh(ranges, (10, 10), facecolors="orange")
    ax.set_ylim(5, 25)
    ax.set_xlim(0, result.hits[-1].end + 1000)
    ax.set_yticks([])
    ax.set_xlabel("Chromosome Position (bp)")
    ax.set_title(
        f"Mobile Elements Map ({result.element_name}) on "
        f"{result.genome_name} — Total: {result.total_matches}"
    )
    fig.tight_layout()

    return _save(fig, output_dir, "chromosome_map")