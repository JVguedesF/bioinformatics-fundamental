from dataclasses import dataclass
from typing import Dict, List, Tuple

class BioPipelineError(Exception):
    pass

class SequenceAnalysisError(BioPipelineError):
    pass

class EmptySequenceError(SequenceAnalysisError):
    pass

@dataclass
class ORFResult:
    id: str
    frame: int
    strand: int
    length_aa: int
    length_bp: int
    protein_seq: str
    start_pos: int


@dataclass
class CodonMetrics:
    total_codons: int
    gc_content: float
    codon_counts: Dict[str, int]
    amino_acid_counts: Dict[str, int]
    top_codons: List[Tuple[str, int]]