import os
from pathlib import Path
from dataclasses import dataclass, field
from typing import Dict, List
from dotenv import load_dotenv

load_dotenv()


@dataclass
class AppConfig:
    PROJECT_ROOT: Path = Path(__file__).resolve().parent.parent

    FILE_EXTENSIONS: List[str] = field(default_factory=lambda: ["*.fasta", "*.fa"])
    ENTREZ_EMAIL: str = field(default_factory=lambda: os.environ.get("ENTREZ_EMAIL", ""))

    MIN_PROTEIN_LENGTH: int = 5
    PREVIEW_LENGTH: int = 60

    TYPE_COLORS: Dict[str, str] = field(default_factory=lambda: {
        'DNA': 'blue',
        'RNA': 'magenta',
        'PROTEIN': 'green',
        'UNKNOWN': 'red'
    })

    DATASETS: List[Dict[str, str]] = field(default_factory=list)

    def __post_init__(self):
        self.DATA_DIR = self.PROJECT_ROOT / "data" / "sequences"
        self.RESULTS_DIR = self.PROJECT_ROOT / "results"
        self.DATA_DIR.mkdir(parents=True, exist_ok=True)
        self.RESULTS_DIR.mkdir(parents=True, exist_ok=True)


settings = AppConfig()