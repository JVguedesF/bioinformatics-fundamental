import re
from pathlib import Path
from typing import List, Tuple, Dict

from Bio import SeqIO
from Bio.Seq import UndefinedSequenceError

from .models import GeneDensityResult, GCSkewResult, MobileElementAnalysis, MobileElementHit
from .exceptions import EmptySequenceError, SequenceAnalysisError


def _calculate_density(record) -> Tuple[int, int, float]:
    allowed_types = {"CDS", "tRNA", "rRNA"}
    gene_count = sum(1 for f in record.features if f.type in allowed_types)

    try:
        total_length_bp = len(record.seq)
    except UndefinedSequenceError:
        total_length_bp = (
            max(f.location.end for f in record.features)
            if record.features
            else 0
        )

    if total_length_bp == 0:
        return 0, 0, 0.0

    density_genes_per_kb = gene_count / (total_length_bp / 1000)
    return total_length_bp, gene_count, density_genes_per_kb


def _get_sequence_string(sequence) -> str:
    try:
        if hasattr(sequence, "tobytes"):
            return sequence.tobytes().decode("utf-8").upper()
        return str(sequence).upper()
    except UndefinedSequenceError:
        return ""


def _calculate_skew_arrays(sequence, window: int, step: int) -> Tuple[List[int], List[float]]:
    sequence_upper = _get_sequence_string(sequence)

    if not sequence_upper:
        raise EmptySequenceError("Sequence is empty or undefined — cannot compute GC skew.")

    positions: List[int] = []
    skew_values: List[float] = []

    for i in range(0, len(sequence_upper) - window, step):
        sub = sequence_upper[i: i + window]
        g, c = sub.count("G"), sub.count("C")
        skew = (g - c) / (g + c) if (g + c) else 0.0
        positions.append(i)
        skew_values.append(skew)

    return positions, skew_values


def _accumulate_skew(skew_list: List[float]) -> List[float]:
    cumulative: List[float] = []
    total = 0.0
    for value in skew_list:
        total += value
        cumulative.append(total)
    return cumulative


def _scan_regex(sequence, pattern: str) -> List[MobileElementHit]:
    sequence_upper = _get_sequence_string(sequence)

    if not sequence_upper:
        raise EmptySequenceError("Sequence is empty or undefined — cannot scan for mobile elements.")

    return [
        MobileElementHit(
            start=m.start(),
            end=m.end(),
            sequence_snippet=sequence_upper[m.start(): m.end()],
        )
        for m in re.compile(pattern).finditer(sequence_upper)
    ]


def pipeline_density_analysis(datasets: List[Dict], data_dir: Path) -> List[GeneDensityResult]:
    results: List[GeneDensityResult] = []

    for data in (d for d in datasets if d.get("group") != "nuclear_target"):
        file_path = data_dir / f"{data['id']}.gb"
        if not file_path.exists():
            continue

        try:
            record = SeqIO.read(file_path, "genbank")
            total_len, count, density = _calculate_density(record)
            results.append(GeneDensityResult(
                genome_name=data["name"],
                group=data["group"],
                total_length_bp=total_len,
                gene_count=count,
                density_genes_per_kb=density,
            ))
        except SequenceAnalysisError as e:
            print(f"[density] Erro de análise em {data['name']}: {e}")
        except Exception as e:
            print(f"[density] Erro inesperado em {data['name']}: {e}")

    return sorted(results, key=lambda x: x.density_genes_per_kb, reverse=True)


def pipeline_skew_analysis(datasets: List[Dict], data_dir: Path, window: int, step: int) -> List[GCSkewResult]:
    valid_groups = {"bacteria_model", "organelle_animal", "organelle_plant", "bacteria_parasite"}
    results: List[GCSkewResult] = []

    for data in (d for d in datasets if d.get("group") in valid_groups):
        file_path = data_dir / f"{data['id']}.gb"
        if not file_path.exists():
            continue

        try:
            record = SeqIO.read(file_path, "genbank")
            positions, raw_skew = _calculate_skew_arrays(record.seq, window, step)
            cumulative_skew = _accumulate_skew(raw_skew)

            min_idx = cumulative_skew.index(min(cumulative_skew))
            results.append(GCSkewResult(
                genome_name=data["name"],
                group=data["group"],
                window_size=window,
                step_size=step,
                genomic_positions=positions,
                skew_values=raw_skew,
                cumulative_skew=cumulative_skew,
                predicted_ori_pos=positions[min_idx],
            ))
        except EmptySequenceError as e:
            print(f"[skew] Sequência vazia em {data['name']}: {e}")
        except SequenceAnalysisError as e:
            print(f"[skew] Erro de análise em {data['name']}: {e}")
        except Exception as e:
            print(f"[skew] Erro inesperado em {data['name']}: {e}")

    return results


def pipeline_mobile_elements(datasets: List[Dict], data_dir: Path, pattern: str) -> MobileElementAnalysis | None:
    target_data = next((d for d in datasets if d.get("group") == "nuclear_target"), None)

    if not target_data:
        return None

    file_path = data_dir / f"{target_data['id']}.fasta"
    if not file_path.exists():
        return None

    try:
        record = SeqIO.read(file_path, "fasta")
        hits = _scan_regex(record.seq, pattern)
        return MobileElementAnalysis(
            genome_name=target_data["name"],
            element_name="Alu_Consensus",
            total_matches=len(hits),
            hits=hits,
        )
    except EmptySequenceError as e:
        print(f"[mobile] Sequência vazia: {e}")
        return None
    except Exception as e:
        print(f"[mobile] Erro inesperado: {e}")
        return None