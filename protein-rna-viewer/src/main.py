import sys
from pathlib import Path

ROOT_DIR = Path(__file__).resolve().parents[1]
if str(ROOT_DIR) not in sys.path:
    sys.path.insert(0, str(ROOT_DIR))

from src.analyzer import analyze_sequence
from src.config import settings
from src.pipeline import setup_arguments, setup_directories, run_batch_analysis
from src.reporter import export_results

# Garante que os caminhos do config sejam absolutos em relação à raiz do projeto
settings.BASE_DIR = ROOT_DIR
settings.RESULTS_DIR = ROOT_DIR / "results"
settings.DATA_DIR = ROOT_DIR / "data" / "sequences"


def format_analysis_metrics(res) -> str:
    """Centraliza a formatação de métricas por tipo de molécula para evitar duplicação."""
    if res.molecule_type == 'DNA':
        gc = f"{res.gc_content:.1f}%" if res.gc_content is not None else "N/A"
        tm = f"{res.melting_temp:.1f}°C" if res.melting_temp is not None else "N/A"
        return f"    GC%: {gc}  |  Tm: {tm}"
    
    if res.molecule_type == 'RNA':
        mfe = f"{res.mfe:.1f} kcal/mol" if res.mfe is not None else "N/A"
        structure = f"\n    Estrutura: {res.secondary_structure[:40]}..." if res.secondary_structure else ""
        return f"    MFE: {mfe}{structure}"
    
    # Default para Proteína
    pi = f"{res.isoelectric_point:.2f}" if res.isoelectric_point is not None else "N/A"
    mw = f"{res.molecular_weight:.0f} Da" if res.molecular_weight is not None else "N/A"
    stability = ""
    if res.stability_index is not None:
        stability = f"  |  Instabilidade: {res.stability_index:.2f}"
    return f"    pI: {pi}  |  MW: {mw}{stability}"


def print_visualization(result, sequence_preview):
    length_seq = 60
    seq_short = sequence_preview[:length_seq]
    if len(sequence_preview) > length_seq:
        seq_short += "..."

    print("\n" + "=" * 60)
    print(f"▶ ANALISANDO: {result.filename}")
    print("=" * 60)
    print(f"  ID:   {result.sequence_id}")
    print(f"  TIPO: {result.molecule_type} ({result.length} resíduos)")
    print(f"  SEQ:  {seq_short}")
    print(f"\n  [ Métricas Principais ]\n{format_analysis_metrics(result)}")
    print("\n" + "-" * 60)


def print_summary(results):
    print("\n\n" + "#" * 60)
    print("  RESUMO DA ANÁLISE COMPARATIVA")
    print("#" * 60)

    for res in results:
        print(f"\n> Arquivo: {res.filename} | Tipo: {res.molecule_type} | Tamanho: {res.length}")
        print(format_analysis_metrics(res))
    print("\n" + "#" * 60)


def process_file(record, file_path):
    result = analyze_sequence(record, file_path.name)
    print_visualization(result, str(record.seq))
    return result


def handle_analysis_completion(results, out_dir):
    print_summary(results)
    print(f"[*] Exportando resultados para: {out_dir.absolute()}")
    export_results(results, out_dir)


def main():
    name = "BioPipeline"
    version = "v1.0"
    args = setup_arguments(name, version)
    setup_directories(args, settings)

    run_batch_analysis(
        name=name,
        version=version,
        config=settings,
        process_func=process_file,
        post_hook=handle_analysis_completion,
    )


if __name__ == "__main__":
    main()