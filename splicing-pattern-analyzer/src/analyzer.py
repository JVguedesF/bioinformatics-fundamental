from typing import TYPE_CHECKING

from Bio.Seq import Seq
from Bio.SeqUtils import gc_fraction

from .models import Intron
from .exceptions import EmptySequenceError, SequenceAnalysisError

if TYPE_CHECKING:
    from Bio.Seq import Seq as SeqType


def _calculate_gc_percent(sequence: str) -> float:
    if not sequence:
        return 0.0
    return gc_fraction(sequence) * 100


class SplicingAnalyzer:
    CHUNK_SIZE: int = 15
    MIN_CHUNK_LEN: int = 5
    DONOR_SITE_LEN: int = 2
    ACCEPTOR_SITE_LEN: int = 2

    def analyze(self, mrna_seq: str | Seq, genomic_seq: str | Seq) -> list[Intron]:
        """
        Analyzes mRNA and genomic sequences to identify introns.

        Args:
            mrna_seq: The mRNA sequence.
            genomic_seq: The genomic DNA sequence.

        Returns:
            A list of identified Intron objects.
        """
        mrna = str(mrna_seq).upper()
        genomic = str(genomic_seq).upper()

        if not mrna:
            raise EmptySequenceError("mRNA sequence is empty.")
        if not genomic:
            raise EmptySequenceError("Genomic sequence is empty.")

        gen_idx = genomic.find(mrna[: self.CHUNK_SIZE])
        if gen_idx == -1:
            raise SequenceAnalysisError("mRNA start not found in genomic sequence.")

        rna_idx = 0
        introns: list[Intron] = []
        genomic_len = len(genomic)

        while rna_idx < len(mrna):
            # Safety check to prevent IndexError if genomic ends before mRNA
            if gen_idx >= genomic_len:
                break

            if genomic[gen_idx] == mrna[rna_idx]:
                rna_idx += 1
                gen_idx += 1
            else:
                chunk = mrna[rna_idx : rna_idx + self.CHUNK_SIZE]
                if len(chunk) < self.MIN_CHUNK_LEN:
                    break

                next_exon_start = genomic.find(chunk, gen_idx)
                if next_exon_start == -1:
                    break

                introns.append(self._create_intron(
                    intron_id=len(introns) + 1,
                    start=gen_idx,
                    end=next_exon_start,
                    genomic_seq=genomic,
                ))

                gen_idx = next_exon_start

        return introns

    def _create_intron(self, intron_id: int, start: int, end: int, genomic_seq: str) -> Intron:
        sequence = genomic_seq[start:end]
        donor = sequence[: self.DONOR_SITE_LEN]
        acceptor = sequence[-self.ACCEPTOR_SITE_LEN :]

        return Intron(
            id=intron_id,
            start=start,
            end=end,
            length=len(sequence),
            sequence=sequence,
            donor_site=donor,
            acceptor_site=acceptor,
            is_canonical=(donor == "GT" and acceptor == "AG"),
            gc_content=_calculate_gc_percent(sequence),
        )