import sys
import warnings
from pathlib import Path

from Bio import BiopythonWarning

warnings.simplefilter('ignore', BiopythonWarning)

ROOT_DIR = Path(__file__).resolve().parents[1]
if str(ROOT_DIR) not in sys.path:
    sys.path.insert(0, str(ROOT_DIR))

from src.analyzer import analyze_codon_usage, find_orfs
from src.config import settings
from src.downloader import download_sequences
from src.pipeline import run_batch_analysis, setup_arguments, setup_directories
from src.reporter import export_results

def print_visualization(filename, record_id, length, table_id, orfs, metrics):
    org_type = "BACTERIA (Table 11)" if table_id == 11 else "STANDARD (Table 1)"

    print("\n" + "=" * 60)
    print(f"▶ ANALISANDO: {filename}")
    print("=" * 60)
    print(f"  ID:   {record_id}")
    print(f"  TIPO: {org_type} ({length} resíduos)")

    print("\n  [ Métricas Genômicas ]")
    print(f"    - Conteúdo GC: {metrics.gc_content:.2f}%")
    print(f"    - Total de Códons: {metrics.total_codons}")

    print("\n  [ Predição de ORFs ]")
    print(f"    - Total detectado (>100aa): {len(orfs)}")
    if orfs:
        best = orfs[0]
        print(f"    - Maior Candidato: {best.length_aa} aa (Frame {best.frame} | Strand {best.strand})")
        print(f"    - Preview: {best.protein_seq[:25]}...")

    print("\n  [ Perfil de Uso de Códons (Top 3) ]")
    for c, n in metrics.top_codons[:3]:
        print(f"    - {c}: {n} ocorrências")

    print("\n" + "-" * 60)

def print_summary(results):
    print("\n\n" + "#" * 60)
    print("  RESUMO DA ANÁLISE DE CÓDONS E ORFs")
    print("#" * 60)

    for res in results:
        top_c = ", ".join(x[0] for x in res['top_codons'][:3])
        org = "BACTERIA" if res['table'] == 11 else "STANDARD"
        print(f"\n> Arquivo: {res['file']} | Tipo: {org} | Tamanho: {res['length']}")
        print(f"    GC%: {res['gc_percent']:.1f}%  |  ORFs: {res['orf_count']}  |  Top Códons: {top_c}")
    print()

def process_file(record, file_path):
    dataset_info = next((d for d in settings.DATASETS if d["name"] == file_path.name), None)
    table_id = dataset_info.get("table", 1) if dataset_info else 1

    orfs = find_orfs(record.id, record.seq, table_id, settings.MIN_ORF_LENGTH)
    metrics = analyze_codon_usage(record.seq, table_id)

    print_visualization(
        filename=file_path.name,
        record_id=record.id,
        length=len(record.seq),
        table_id=table_id,
        orfs=orfs,
        metrics=metrics,
    )

    return {
        "id": record.id,
        "file": file_path.name,
        "table": table_id,
        "length": len(record.seq),
        "orf_count": len(orfs),
        "gc_percent": metrics.gc_content,
        "top_codons": metrics.top_codons,
    }

def main():
    name = "GeneCodePro"
    version = "v1.0"
    args = setup_arguments(name, version)
    setup_directories(args, settings)

    run_batch_analysis(
        name=name,
        version=version,
        config=settings,
        process_func=process_file,
        pre_hook=lambda out_dir: download_sequences(settings.ENTREZ_EMAIL, settings.DATASETS, out_dir),
        post_hook=lambda results, out_dir: (
            print_summary(results),
            export_results(results, out_dir),
        ),
    )

if __name__ == "__main__":
    main()