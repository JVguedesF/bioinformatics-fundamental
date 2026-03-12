import sys
import argparse
from pathlib import Path

from Bio import SeqIO

ROOT_DIR = Path(__file__).resolve().parents[1]
if str(ROOT_DIR) not in sys.path:
    sys.path.insert(0, str(ROOT_DIR))

from src.config import settings
from src.analyzer import SplicingAnalyzer
from src.downloader import download_datasets
from src.reporter import SplicingReporter
from src.models import BioPipelineError


def _parse_arguments():
    parser = argparse.ArgumentParser(description="Eukaryotic Splicing Analyzer")
    parser.add_argument("--input", "-i", type=Path)
    parser.add_argument("--output", "-o", type=Path)
    return parser.parse_args()


def _get_dataset_path(name_pattern: str) -> Path:
    for d in settings.DATASETS:
        if name_pattern in d["name"] or name_pattern in d["id"]:
            return settings.DATA_DIR / d["name"]
    raise FileNotFoundError(f"Dataset containing '{name_pattern}' not found.")


def print_visualization(label, introns, genomic_name, mrna_name):
    print("\n" + "=" * 60)
    print(f"▶ ANALISANDO: {label}")
    print("=" * 60)
    print(f"  Referência Genômica: {genomic_name}")
    print(f"  Alvo mRNA: {mrna_name}")
    print(f"  Íntrons Detectados: {len(introns)}")

    if not introns:
        print("\n  [ Aviso: Nenhum íntron mapeado nesta sequência ]")
        return

    print("\n  [ Mapa de Íntrons (Top 5) ]")

    canonical_count = sum(1 for i in introns if i.is_canonical)

    print(f"    - Cânonicos (GT-AG): {canonical_count}/{len(introns)}")
    print(f"    - Média de GC: {sum(i.gc_content for i in introns) / len(introns):.2f}%\n")

    print(f"    {'ID':<5} | {'INÍCIO':<10} | {'FIM':<10} | {'TAMANHO':<8} | SÍTIOS (5'...3')")
    print("    " + "-" * 55)

    for i in introns[:5]:
        sites = f"{i.donor_site}...{i.acceptor_site}"
        print(f"    {i.id:<5} | {i.start:<10} | {i.end:<10} | {i.length:<8} | {sites}")

    if len(introns) > 5:
        print(f"    ... e mais {len(introns) - 5} íntrons.")


def main():
    args = _parse_arguments()

    if args.input:
        settings.DATA_DIR = args.input
    if args.output:
        settings.RESULTS_DIR = args.output

    print("=" * 60)
    print("Eukaryotic Splicing Analyzer v1.0")
    print("=" * 60)

    print("Verificando arquivos...")
    download_datasets(settings.DATA_DIR)

    try:
        genomic_path = _get_dataset_path("genomic_ref")
        genomic_record = SeqIO.read(genomic_path, "fasta")

        analysis_map = [
            ("NM_000546", "Human TP53 (Variant 1)"),
            ("NM_001126112", "Human TP53 (Variant 2)"),
        ]

        analyzer = SplicingAnalyzer()
        reporter = SplicingReporter(settings.RESULTS_DIR)

        for mrna_pattern, label in analysis_map:
            mrna_path = _get_dataset_path(mrna_pattern)
            mrna_record = SeqIO.read(mrna_path, "fasta")

            introns = analyzer.analyze(mrna_record.seq, genomic_record.seq)

            print_visualization(label, introns, genomic_path.name, mrna_path.name)

            print("\n  [ Exportação de Dados ]")
            reporter.generate_reports(label, introns)

    except BioPipelineError as e:
        print(f"\n[ERRO DE PIPELINE] {e}")
        sys.exit(1)
    except Exception as e:
        print(f"\n[ERRO INESPERADO] {e}")
        sys.exit(1)

    print("\n" + "=" * 60)
    print("Pipeline Finalizado!")
    print("=" * 60)


if __name__ == "__main__":
    main()