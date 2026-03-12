import sys
import argparse
from pathlib import Path

ROOT_DIR = Path(__file__).resolve().parents[1]
if str(ROOT_DIR) not in sys.path:
    sys.path.insert(0, str(ROOT_DIR))

from src.config import settings
from src.reporter import BioReporter
from src.printer import TerminalPrinter
from src import downloader, analyzer, view

def _check_missing(datasets, data_dir):
    missing = False
    for data in datasets:
        ext = "fasta" if data.get("group") == "nuclear_target" else "gb"
        if not (data_dir / f"{data['id']}.{ext}").exists():
            print(f"  [Aviso] Arquivo ausente: {data['id']}.{ext}")
            missing = True
    return missing

def _run_density_analysis(graphs_dir):
    TerminalPrinter.section("1. Análise de Densidade Gênica")
    density_data = analyzer.pipeline_density_analysis(settings.DATASETS, settings.DATA_DIR)
    if density_data:
        TerminalPrinter.density_table(density_data)
        BioReporter.export_density(density_data, settings.RESULTS_DIR)
        if input("\nGerar gráfico de densidade? (s/n): ").lower() == "s":
            view.plot_density_comparison(density_data, graphs_dir)

def _run_skew_analysis(graphs_dir):
    TerminalPrinter.section("2. Análise de GC Skew")
    skew_data = analyzer.pipeline_skew_analysis(settings.DATASETS, settings.DATA_DIR, settings.GC_SKEW_WINDOW, settings.GC_SKEW_STEP)
    if skew_data:
        for res in skew_data:
            TerminalPrinter.skew_result(res.genome_name, res.predicted_ori_pos)
        BioReporter.export_skew(skew_data, settings.RESULTS_DIR)
        if input("\nGerar gráficos de GC Skew? (s/n): ").lower() == "s":
            for res in skew_data:
                view.plot_gc_skew_curve(res, graphs_dir)

def _run_mobile_analysis(graphs_dir):
    TerminalPrinter.section("3. Rastreamento de Elementos Móveis")
    mobile_data = analyzer.pipeline_mobile_elements(settings.DATASETS, settings.DATA_DIR, settings.ALU_CONSENSUS_SEQ)
    if mobile_data:
        TerminalPrinter.mobile_summary(mobile_data.genome_name, mobile_data.total_matches)
        BioReporter.export_mobile_elements(mobile_data, settings.RESULTS_DIR)
        if input("\nGerar mapa cromossômico? (s/n): ").lower() == "s":
            view.plot_chromosome_map(mobile_data, graphs_dir)

def main():
    parser = argparse.ArgumentParser(description="Genomic Architecture Analyzer")
    parser.add_argument("--download", action="store_true")
    args = parser.parse_args()

    graphs_dir = settings.RESULTS_DIR / "graphs"
    graphs_dir.mkdir(parents=True, exist_ok=True)

    TerminalPrinter.header()

    if args.download or _check_missing(settings.DATASETS, settings.DATA_DIR):
        downloader.download_datasets(settings.DATASETS, settings.DATA_DIR)
    else:
        print("  ✓ Todos os arquivos verificados localmente.")

    _run_density_analysis(graphs_dir)
    _run_skew_analysis(graphs_dir)
    _run_mobile_analysis(graphs_dir)

    print("\nPipeline Finalizado!")

if __name__ == "__main__":
    main()