from Bio.SeqRecord import SeqRecord
from Bio.SeqUtils import gc_fraction
from Bio.SeqUtils import MeltingTemp
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from Bio.Data.CodonTable import TranslationError

from src.config import settings
from src.models import AnalysisResult, SequenceAnalysisError, EmptySequenceError

try:
    import RNA
    HAS_VIENNA_RNA = True
except ImportError:
    RNA = None
    HAS_VIENNA_RNA = False

DNA_BASES = set("ATGCN")
RNA_BASES = set("AUGCN")


def _calculate_dna_metrics(seq) -> dict:
    return {
        "length": len(seq),
        "gc_content": gc_fraction(seq) * 100,
        "tm": MeltingTemp.Tm_NN(seq, strict=False),
    }


def _get_folding_data(sequence_str: str):
    if HAS_VIENNA_RNA:
        structure, mfe = RNA.fold(sequence_str)
        return structure, mfe
    return None, None


def _analyze_protein(sequence_str: str) -> dict | None:
    clean_seq = sequence_str.rstrip('*')
    if not clean_seq:
        return None
    try:
        protein_analysis = ProteinAnalysis(clean_seq)
        return {
            "length": len(clean_seq),
            "molecular_weight": protein_analysis.molecular_weight(),
            "isoelectric_point": protein_analysis.isoelectric_point(),
            "gravy": protein_analysis.gravy(),
            "instability_index": protein_analysis.instability_index(),
        }
    except ValueError:
        return None


# ── Analisador principal ──────────────────────────────────────────────────────

class CentralDogmaAnalyzer:

    @staticmethod
    def detect_molecule_type(sequence: str) -> str:
        unique_chars = set(sequence.upper()) - set("\n\r\t ")
        if len(unique_chars - DNA_BASES - RNA_BASES) > 0:
            return 'PROTEIN'
        if 'T' in unique_chars and 'U' not in unique_chars:
            return 'DNA'
        if 'U' in unique_chars and 'T' not in unique_chars:
            return 'RNA'
        return 'DNA'

    def process_sequence(self, record: SeqRecord, filename: str) -> AnalysisResult:
        if not record.seq:
            raise EmptySequenceError(f"Sequence in {filename} is empty.")

        mol_type = self.detect_molecule_type(str(record.seq)[:1000])
        result = AnalysisResult(
            filename=filename,
            sequence_id=record.id,
            molecule_type=mol_type,
            length=len(record.seq),
        )

        rna_seq = None
        protein_seq = record.seq

        try:
            if mol_type == 'DNA':
                dna = _calculate_dna_metrics(record.seq)
                result.gc_content = dna.get('gc_content')
                result.melting_temp = dna.get('tm')
                rna_seq = record.seq.transcribe()
            elif mol_type == 'RNA':
                rna_seq = record.seq

            if rna_seq and HAS_VIENNA_RNA:
                structure, mfe = _get_folding_data(str(rna_seq))
                if structure:
                    result.mfe = mfe
                    result.secondary_structure = structure
                try:
                    protein_seq = rna_seq.translate(to_stop=True)
                except (ValueError, TranslationError):
                    protein_seq = ""

            if len(protein_seq) >= settings.MIN_PROTEIN_LENGTH:
                prot = _analyze_protein(str(protein_seq))
                if prot:
                    result.isoelectric_point = prot.get('isoelectric_point')
                    result.molecular_weight = prot.get('molecular_weight')
                    result.stability_index = prot.get('instability_index')

        except Exception as e:
            raise SequenceAnalysisError(f"Failed to analyze {mol_type} sequence: {e}")

        return result