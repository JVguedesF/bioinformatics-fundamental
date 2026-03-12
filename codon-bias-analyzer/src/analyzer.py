from collections import Counter

from Bio.Data import CodonTable
from Bio.SeqUtils import gc_fraction

from src.models import CodonMetrics, EmptySequenceError, ORFResult, SequenceAnalysisError

def _scan_frame(sequence_id, seq, frame, strand, min_len, table_id, results):
    trans = seq[frame:].translate(table=table_id)
    current_aa_pos = 0
    for protein in trans.split("*"):
        prot_len = len(protein)
        if prot_len >= min_len:
            results.append(ORFResult(
                id=f"{sequence_id}_fr{frame}_st{strand}",
                frame=frame + 1,
                strand=strand,
                length_aa=prot_len,
                length_bp=prot_len * 3,
                protein_seq=str(protein),
                start_pos=(current_aa_pos * 3) + frame,
            ))
        current_aa_pos += prot_len + 1

def find_orfs(sequence_id: str, sequence, table_id: int = 1, min_len_aa: int = 100):
    if not sequence:
        raise EmptySequenceError(f"Sequence {sequence_id} is empty.")
    try:
        results = []
        targets = [(1, sequence), (-1, sequence.reverse_complement())]
        for strand, nuc_seq in targets:
            for frame in range(3):
                _scan_frame(sequence_id, nuc_seq, frame, strand, min_len_aa, table_id, results)
        return sorted(results, key=lambda x: x.length_aa, reverse=True)
    except Exception as e:
        raise SequenceAnalysisError(f"Error finding ORFs: {e}")

def analyze_codon_usage(sequence, table_id: int = 1):
    try:
        trim = len(sequence) // 3 * 3
        coding_seq = sequence[:trim]
        str_seq = str(coding_seq)
        codons = [str_seq[i:i + 3] for i in range(0, len(str_seq), 3)]
        codon_counts = Counter(codons)
        aa_counts = Counter(coding_seq.translate(table=table_id))
        return CodonMetrics(
            total_codons=sum(codon_counts.values()),
            gc_content=gc_fraction(sequence) * 100,
            codon_counts=dict(codon_counts),
            amino_acid_counts=dict(aa_counts),
            top_codons=codon_counts.most_common(5),
        )
    except Exception as e:
        raise SequenceAnalysisError(f"Error analyzing codon usage: {e}")