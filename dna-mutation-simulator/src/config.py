import os
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List

from dotenv import load_dotenv

load_dotenv()


@dataclass
class AppConfig:
    PROJECT_ROOT: Path = Path(__file__).resolve().parent.parent

    FILE_EXTENSIONS: List[str] = field(default_factory=lambda: ["*.fasta", "*.fa"])
    ENTREZ_EMAIL: str = field(default_factory=lambda: os.environ.get("ENTREZ_EMAIL", ""))

    DATASETS: List[Dict[str, str]] = field(default_factory=lambda: [
        {"id": "NC_000913", "name": "Escherichia_coli_K12.fasta",       "type": "fasta", "group": "replication_model"},
        {"id": "NM_000251", "name": "Homo_sapiens_MSH2.fasta",          "type": "fasta", "group": "repair_genes"},
        {"id": "NM_007294", "name": "Homo_sapiens_BRCA1.fasta",         "type": "fasta", "group": "recombination_genes"},
        {"id": "NC_001133", "name": "Saccharomyces_cerevisiae_Chr1.fasta","type": "fasta", "group": "crossing_over_model"},
    ])

    def __post_init__(self):
        self.DATA_DIR = self.PROJECT_ROOT / "data" / "sequences"
        self.RESULTS_DIR = self.PROJECT_ROOT / "results"
        self.DATA_DIR.mkdir(parents=True, exist_ok=True)
        self.RESULTS_DIR.mkdir(parents=True, exist_ok=True)


settings = AppConfig()
