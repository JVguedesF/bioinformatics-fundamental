import os
from pathlib import Path
from dataclasses import dataclass, field
from typing import List, Dict
from dotenv import load_dotenv

load_dotenv()


@dataclass
class AppConfig:
    PROJECT_ROOT: Path = Path(__file__).resolve().parent.parent

    FILE_EXTENSIONS: List[str] = field(default_factory=lambda: ["*.fasta", "*.fa"])
    ENTREZ_EMAIL: str = field(default_factory=lambda: os.environ.get("ENTREZ_EMAIL", ""))

    MIN_ORF_LENGTH: int = 100

    DATASETS: List[Dict[str, str]] = field(default_factory=lambda: [
        {"id": "NC_001416.1", "name": "lambda_phage.fasta", "type": "fasta"},
        {"id": "NC_000913.3", "name": "e_coli_k12.fasta", "type": "fasta"},
        {"id": "NM_000207.3", "name": "human_insulin.fasta", "type": "fasta"},
    ])

    def __post_init__(self):
        self.DATA_DIR = self.PROJECT_ROOT / "data" / "sequences"
        self.RESULTS_DIR = self.PROJECT_ROOT / "results"
        self.DATA_DIR.mkdir(parents=True, exist_ok=True)
        self.RESULTS_DIR.mkdir(parents=True, exist_ok=True)


settings = AppConfig()