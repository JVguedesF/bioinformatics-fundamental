from dataclasses import dataclass

class BioPipelineError(Exception):
    pass

class SequenceAnalysisError(BioPipelineError):
    pass

class EmptySequenceError(SequenceAnalysisError):
    pass


@dataclass
class Intron:
    id: int
    start: int
    end: int
    length: int
    sequence: str
    donor_site: str
    acceptor_site: str
    is_canonical: bool
    gc_content: float