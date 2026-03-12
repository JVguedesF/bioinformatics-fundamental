import csv
import json
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List

def _ts():
    return datetime.now().strftime("%Y%m%d_%H%M%S")

class BioReporter:
    @staticmethod
    def export_density(results, output_dir: Path) -> None:
        rows = [{"genome_name": r.genome_name, "group": r.group, "density": round(r.density_genes_per_kb, 4)} for r in results]
        path = output_dir / f"density_{_ts()}.json"
        with open(path, "w") as f: json.dump(rows, f, indent=2)
        print(f"    - Exportado: {path.name}")

    @staticmethod
    def export_skew(results, output_dir: Path) -> None:
        rows = [{"genome": r.genome_name, "ori": r.predicted_ori_pos} for r in results]
        path = output_dir / f"gc_skew_{_ts()}.json"
        with open(path, "w") as f: json.dump(rows, f, indent=2)
        print(f"    - Exportado: {path.name}")

    @staticmethod
    def export_mobile_elements(result, output_dir: Path) -> None:
        summary = {"genome": result.genome_name, "total": result.total_matches}
        path = output_dir / f"mobile_{_ts()}.json"
        with open(path, "w") as f: json.dump(summary, f, indent=2)
        print(f"    - Exportado: {path.name}")