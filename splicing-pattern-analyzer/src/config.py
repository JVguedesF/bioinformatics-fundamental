import os
from dataclasses import dataclass, field
from pathlib import Path
from typing import List, Dict

from dotenv import load_dotenv

load_dotenv()

CURRENT_PROJECT_DIR = Path(__file__).resolve().parent.parent


@dataclass
class AppConfig:
    PROJECT_ROOT: Path

    ENTREZ_EMAIL: str = field(default_factory=lambda: os.environ.get("ENTREZ_EMAIL"))
    MIN_ORF_LENGTH: int = 100

    DATASETS: List[Dict[str, str]] = field(default_factory=lambda: [
        {"id": "NG_017013.2", "name": "tp53_human_genomic_refseqgene.fasta", "type": "fasta"},
        {"id": "NM_000546.6", "name": "tp53_human_transcript_variant1.fasta", "type": "fasta"},
        {"id": "NM_001126112.2", "name": "tp53_human_transcript_variant2.fasta", "type": "fasta"},
        {"id": "NG_007557.1", "name": "trp53_mouse_genomic_ortholog.fasta", "type": "fasta"},
    ])

    def __post_init__(self):
        self.DATA_DIR = self.PROJECT_ROOT / "data" / "sequences"
        self.RESULTS_DIR = self.PROJECT_ROOT / "results"

        self.DATA_DIR.mkdir(parents=True, exist_ok=True)
        self.RESULTS_DIR.mkdir(parents=True, exist_ok=True)


settings = AppConfig(PROJECT_ROOT=CURRENT_PROJECT_DIR)