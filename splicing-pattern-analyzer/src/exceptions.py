class BioPipelineError(Exception):
    pass


class SequenceAnalysisError(BioPipelineError):
    pass


class EmptySequenceError(SequenceAnalysisError):
    pass