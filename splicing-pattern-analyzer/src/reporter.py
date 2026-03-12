import csv
import json
from dataclasses import asdict
from datetime import datetime
from pathlib import Path
from typing import List

from src.models import Intron


class SplicingReporter:
    def __init__(self, output_dir: Path):
        self.output_dir = output_dir
        self.output_dir.mkdir(parents=True, exist_ok=True)

    def generate_reports(self, label: str, introns: List[Intron]) -> Path:
        clean_label = (
            label.replace(" ", "_").replace("(", "").replace(")", "").lower()
        )
        base_path = self.output_dir / clean_label

        self._save_json(base_path.with_suffix(".json"), label, introns)
        self._save_csv(base_path.with_suffix(".csv"), introns)

        return self.output_dir

    @staticmethod
    def _save_json(path: Path, label: str, introns: List[Intron]) -> None:
        data = {
            "analysis_target": label,
            "timestamp": datetime.now().isoformat(),
            "total_introns": len(introns),
            "introns": [asdict(i) for i in introns],
        }
        with open(path, "w", encoding="utf-8") as f:
            json.dump(data, f, indent=4)
        print(f"    - JSON: {path.name}")

    @staticmethod
    def _save_csv(path: Path, introns: List[Intron]) -> None:
        with open(path, "w", newline="", encoding="utf-8") as f:
            writer = csv.writer(f)
            writer.writerow([
                "ID", "Start", "End", "Length",
                "GC_Content", "Is_Canonical", "Donor", "Acceptor",
            ])
            for i in introns:
                writer.writerow([
                    i.id, i.start, i.end, i.length,
                    f"{i.gc_content:.2f}", i.is_canonical,
                    i.donor_site, i.acceptor_site,
                ])
        print(f"    - CSV:  {path.name}")