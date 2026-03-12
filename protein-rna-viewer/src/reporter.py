import csv
import json
from dataclasses import asdict
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional

from src.models import AnalysisResult

def _ts() -> str:
    return datetime.now().strftime("%Y%m%d_%H%M%S")

def export_json(data: List[Dict], output_dir: Path, prefix: str = "data") -> Path:
    path = output_dir / f"{prefix}_{_ts()}.json"
    with open(path, 'w') as f:
        json.dump(data, f, indent=2)
    print(f"  JSON:  {path.name}")
    return path

def export_csv(data: List[Dict], output_dir: Path, prefix: str = "data") -> Optional[Path]:
    if not data:
        return None
    path = output_dir / f"{prefix}_{_ts()}.csv"
    with open(path, 'w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=list(data[0].keys()))
        writer.writeheader()
        writer.writerows(data)
    print(f"  CSV:   {path.name}")
    return path

def export_results(results: List[AnalysisResult], output_dir: Path):
    output_dir.mkdir(exist_ok=True)
    print("\nResults exported to:")
    data = [asdict(r) for r in results]
    export_json(data, output_dir)
    export_csv(data, output_dir)