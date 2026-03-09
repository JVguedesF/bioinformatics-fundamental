import csv
import json
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List

from rich.console import Console

from .models import GeneDensityResult, GCSkewResult, MobileElementAnalysis

console = Console()


class BioReporter:

    @staticmethod
    def _timestamped(output_dir: Path, prefix: str, ext: str) -> Path:
        ts = datetime.now().strftime("%Y%m%d_%H%M%S")
        return output_dir / f"{prefix}_{ts}.{ext}"

    @staticmethod
    def export_density(results: List[GeneDensityResult], output_dir: Path) -> None:
        if not results:
            return

        rows = [
            {
                "genome_name": r.genome_name,
                "group": r.group,
                "total_length_bp": r.total_length_bp,
                "gene_count": r.gene_count,
                "density_genes_per_kb": round(r.density_genes_per_kb, 4),
            }
            for r in results
        ]

        csv_path = BioReporter._timestamped(output_dir, "density", "csv")
        BioReporter._write_csv(rows, csv_path)

        json_path = BioReporter._timestamped(output_dir, "density", "json")
        BioReporter._write_json(rows, json_path)

        console.print(f"[dim]  ↳ density → {csv_path.name}, {json_path.name}[/dim]")

    @staticmethod
    def export_skew(results: List[GCSkewResult], output_dir: Path) -> None:
        if not results:
            return

        rows = [
            {
                "genome_name": r.genome_name,
                "group": r.group,
                "window_size": r.window_size,
                "step_size": r.step_size,
                "predicted_ori_pos": r.predicted_ori_pos,
            }
            for r in results
        ]

        csv_path = BioReporter._timestamped(output_dir, "gc_skew", "csv")
        BioReporter._write_csv(rows, csv_path)

        full = [
            {
                **row,
                "genomic_positions": r.genomic_positions,
                "skew_values": r.skew_values,
                "cumulative_skew": r.cumulative_skew,
            }
            for row, r in zip(rows, results)
        ]
        json_path = BioReporter._timestamped(output_dir, "gc_skew", "json")
        BioReporter._write_json(full, json_path)

        console.print(f"[dim]  ↳ gc_skew  → {csv_path.name}, {json_path.name}[/dim]")

    @staticmethod
    def export_mobile_elements(result: MobileElementAnalysis, output_dir: Path) -> None:
        if not result:
            return

        rows = [
            {
                "genome_name": result.genome_name,
                "element_name": result.element_name,
                "hit_index": i,
                "start": hit.start,
                "end": hit.end,
                "sequence_snippet": hit.sequence_snippet,
            }
            for i, hit in enumerate(result.hits)
        ]

        csv_path = BioReporter._timestamped(output_dir, "mobile_elements", "csv")
        BioReporter._write_csv(rows, csv_path)

        summary = {
            "genome_name": result.genome_name,
            "element_name": result.element_name,
            "total_matches": result.total_matches,
            "hits": [{"start": h.start, "end": h.end, "snippet": h.sequence_snippet} for h in result.hits],
        }
        json_path = BioReporter._timestamped(output_dir, "mobile_elements", "json")
        BioReporter._write_json(summary, json_path)

        console.print(f"[dim]  ↳ mobile   → {csv_path.name}, {json_path.name}[/dim]")

    @staticmethod
    def _write_csv(rows: List[Dict[str, Any]], path: Path) -> None:
        if not rows:
            return
        with open(path, "w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
            writer.writeheader()
            writer.writerows(rows)

    @staticmethod
    def _write_json(data: Any, path: Path) -> None:
        with open(path, "w", encoding="utf-8") as f:
            json.dump(data, f, indent=2, ensure_ascii=False)