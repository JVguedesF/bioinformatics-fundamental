from dataclasses import dataclass
from datetime import datetime
from typing import Optional

class BioPipelineError(Exception):
    pass

class SequenceAnalysisError(BioPipelineError):
    pass

class EmptySequenceError(SequenceAnalysisError):
    pass

@dataclass
class AnalysisResult:
    filename: str
    sequence_id: str
    molecule_type: str
    length: int
    gc_content: Optional[float] = None
    melting_temp: Optional[float] = None
    mfe: Optional[float] = None
    isoelectric_point: Optional[float] = None
    molecular_weight: Optional[float] = None
    stability_index: Optional[float] = None
    secondary_structure: Optional[str] = None
    timestamp: Optional[str] = None

    def __post_init__(self):
        if self.timestamp is None:
            self.timestamp = datetime.now().isoformat()