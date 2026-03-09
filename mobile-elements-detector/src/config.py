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

    GC_SKEW_WINDOW: int = 1000
    GC_SKEW_STEP: int = 100

    ALU_CONSENSUS_SEQ: str = "TGTAATCCCAGCACTTT|AAAGTGCTGGGATTACA"

    DATASETS: List[Dict[str, str]] = field(default_factory=lambda: [
        {"id": "NC_012920", "name": "Homo_sapiens_Mitochondrion", "type": "genbank", "group": "organelle_animal"},
        {"id": "NC_001807", "name": "Homo_sapiens_Mitochondrion_Ancient_Reference", "type": "genbank", "group": "organelle_animal"},
        {"id": "U00096", "name": "Escherichia_coli_K_12", "type": "genbank", "group": "bacteria_model"},
        {"id": "NC_000913", "name": "Escherichia_coli_MG1655", "type": "genbank", "group": "bacteria_model"},
        {"id": "NC_000963", "name": "Rickettsia_prowazekii", "type": "genbank", "group": "bacteria_parasite"},
        {"id": "NC_002947", "name": "Wolbachia_pipientis", "type": "genbank", "group": "bacteria_symbiont"},
        {"id": "NC_000932", "name": "Arabidopsis_thaliana_Chloroplast", "type": "genbank", "group": "organelle_plant"},
        {"id": "NC_001224", "name": "Saccharomyces_cerevisiae_Mitochondrion", "type": "genbank", "group": "organelle_fungi"},
        {"id": "NC_002607", "name": "Halobacterium_salinarum", "type": "genbank", "group": "archaea"},
        {"id": "NC_000854", "name": "Aeropyrum_pernix", "type": "genbank", "group": "archaea"},
        {"id": "NC_001416", "name": "Enterobacteria_phage_lambda", "type": "genbank", "group": "virus"},
        {"id": "NT_113818", "name": "Human_Chr19_Partial_Alu_Rich", "type": "fasta", "group": "nuclear_target"},
    ])

    def __post_init__(self):
        self.DATA_DIR = self.PROJECT_ROOT / "data" / "sequences"
        self.RESULTS_DIR = self.PROJECT_ROOT / "results"

        self.DATA_DIR.mkdir(parents=True, exist_ok=True)
        self.RESULTS_DIR.mkdir(parents=True, exist_ok=True)


settings = AppConfig(PROJECT_ROOT=CURRENT_PROJECT_DIR)