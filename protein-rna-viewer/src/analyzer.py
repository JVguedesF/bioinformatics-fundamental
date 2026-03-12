from Bio.SeqRecord import SeqRecord
from Bio.SeqUtils import gc_fraction, MeltingTemp
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from Bio.Data.CodonTable import TranslationError

from src.config import settings
from src.models import AnalysisResult, SequenceAnalysisError, EmptySequenceError

try:
    import RNA
    _RNA_AVAILABLE = True
except (ImportError, AttributeError):
    _RNA_AVAILABLE = False
    RNA = None # type: ignore

DNA_BASES = set("ATGCN")
RNA_BASES = set("AUGCN")

def _calculate_dna_metrics(seq) -> dict:
    return {
        "length": len(seq),
        "gc_content": gc_fraction(seq) * 100,
        "tm": MeltingTemp.Tm_NN(seq, strict=False),
    }

def _get_folding_data(sequence_str: str) -> tuple[str | None, float | None]:
    """Executa o folding de RNA apenas se a biblioteca ViennaRNA estiver disponível."""
    if _RNA_AVAILABLE and RNA:
        return RNA.fold(sequence_str)
    return None, None

def _analyze_protein(sequence_str: str) -> dict | None:
    clean_seq = sequence_str.rstrip('*')
    if not clean_seq:
        return None
    try:
        pa = ProteinAnalysis(clean_seq)
        return {
            "length": len(clean_seq),
            "molecular_weight": pa.molecular_weight(),
            "isoelectric_point": pa.isoelectric_point(),
            "gravy": pa.gravy(),
            "instability_index": pa.instability_index(),
        }
    except ValueError:
        return None

def detect_molecule_type(sequence: str) -> str:
    unique_chars = set(sequence.upper()) - set("\n\r\t ")
    if len(unique_chars - DNA_BASES - RNA_BASES) > 0:
        return 'PROTEIN'
    if 'T' in unique_chars and 'U' not in unique_chars:
        return 'DNA'
    if 'U' in unique_chars and 'T' not in unique_chars:
        return 'RNA'
    return 'DNA'

def _process_nucleotide_flow(record: SeqRecord, mol_type: str, result: AnalysisResult):
    """Extrai métricas de DNA/RNA e retorna a sequência proteica candidata."""
    if mol_type == 'PROTEIN':
        return record.seq

    rna_seq = record.seq
    if mol_type == 'DNA':
        dna_metrics = _calculate_dna_metrics(record.seq)
        result.gc_content = dna_metrics.get('gc_content')
        result.melting_temp = dna_metrics.get('tm')
        rna_seq = record.seq.transcribe()

    structure, mfe = _get_folding_data(str(rna_seq))
    if structure:
        result.mfe, result.secondary_structure = mfe, structure

    try:
        return rna_seq.translate(to_stop=True)
    except (ValueError, TranslationError):
        return ""


def _apply_protein_metrics(protein_seq, result: AnalysisResult):
    if len(protein_seq) < settings.MIN_PROTEIN_LENGTH:
        return

    prot = _analyze_protein(str(protein_seq))
    if prot:
        result.isoelectric_point = prot.get('isoelectric_point')
        result.molecular_weight = prot.get('molecular_weight')
        result.stability_index = prot.get('instability_index')


def analyze_sequence(record: SeqRecord, filename: str) -> AnalysisResult:
    if not record.seq:
        raise EmptySequenceError(f"Sequence in {filename} is empty.")

    mol_type = detect_molecule_type(str(record.seq)[:1000])
    result = AnalysisResult(
        filename=filename,
        sequence_id=record.id,
        molecule_type=mol_type,
        length=len(record.seq),
    )

    try:
        protein_candidate = _process_nucleotide_flow(record, mol_type, result)
        _apply_protein_metrics(protein_candidate, result)
    except Exception as e:
        raise SequenceAnalysisError(f"Failed to analyze {mol_type} sequence: {e}")

    return result