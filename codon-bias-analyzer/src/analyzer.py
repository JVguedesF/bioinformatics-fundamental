from collections import Counter

from Bio.Data import CodonTable
from Bio.SeqUtils import gc_fraction

from src.models import CodonMetrics, EmptySequenceError, ORFResult, SequenceAnalysisError


class SequenceAnalyzer:
    def __init__(self, sequence_id, sequence, table_id=1):
        if not sequence:
            raise EmptySequenceError(f"Sequence {sequence_id} is empty.")

        self.sequence_id = sequence_id
        self.sequence = sequence
        self.table_id = table_id
        self.table = CodonTable.unambiguous_dna_by_id[table_id]

    def find_orfs(self, min_len_aa=100):
        try:
            results = []
            targets = [(1, self.sequence), (-1, self.sequence.reverse_complement())]
            for strand, nuc_seq in targets:
                for frame in range(3):
                    self._scan_frame(nuc_seq, frame, strand, min_len_aa, results)
            return sorted(results, key=lambda x: x.length_aa, reverse=True)
        except Exception as e:
            raise SequenceAnalysisError(f"Error finding ORFs: {e}")

    def _scan_frame(self, seq, frame, strand, min_len, results):
        trans = seq[frame:].translate(table=self.table_id)
        current_aa_pos = 0
        for protein in trans.split("*"):
            prot_len = len(protein)
            if prot_len >= min_len:
                results.append(ORFResult(
                    id=f"{self.sequence_id}_fr{frame}_st{strand}",
                    frame=frame + 1,
                    strand=strand,
                    length_aa=prot_len,
                    length_bp=prot_len * 3,
                    protein_seq=str(protein),
                    start_pos=(current_aa_pos * 3) + frame,
                ))
            current_aa_pos += prot_len + 1

    def analyze_codon_usage(self):
        try:
            trim = len(self.sequence) // 3 * 3
            coding_seq = self.sequence[:trim]
            str_seq = str(coding_seq)
            codons = [str_seq[i:i + 3] for i in range(0, len(str_seq), 3)]
            codon_counts = Counter(codons)
            aa_counts = Counter(coding_seq.translate(table=self.table_id))
            return CodonMetrics(
                total_codons=sum(codon_counts.values()),
                gc_content=gc_fraction(self.sequence) * 100,
                codon_counts=dict(codon_counts),
                amino_acid_counts=dict(aa_counts),
                top_codons=codon_counts.most_common(5),
            )
        except Exception as e:
            raise SequenceAnalysisError(f"Error analyzing codon usage: {e}")