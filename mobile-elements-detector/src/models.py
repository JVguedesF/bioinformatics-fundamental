from dataclasses import dataclass
from typing import List


@dataclass
class GeneDensityResult:
    genome_name: str
    group: str
    total_length_bp: int
    gene_count: int
    density_genes_per_kb: float


@dataclass
class GCSkewResult:
    genome_name: str
    group: str
    window_size: int
    step_size: int
    genomic_positions: List[int]
    skew_values: List[float]
    cumulative_skew: List[float]
    predicted_ori_pos: int


@dataclass
class MobileElementHit:
    start: int
    end: int
    sequence_snippet: str


@dataclass
class MobileElementAnalysis:
    genome_name: str
    element_name: str
    total_matches: int
    hits: List[MobileElementHit]