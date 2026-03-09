import logging
import sys
import argparse
from pathlib import Path

from Bio import SeqIO
from rich.console import Console
from rich.logging import RichHandler

ROOT_DIR = Path(__file__).resolve().parents[2]
if str(ROOT_DIR) not in sys.path:
    sys.path.insert(0, str(ROOT_DIR))

from .config import settings
from .analyzer import SplicingAnalyzer
from .downloader import download_datasets
from .reporter import SplicingReporter
from .view import SplicingView
from .exceptions import BioPipelineError

logging.basicConfig(
    level="INFO",
    format="%(message)s",
    datefmt="[%X]",
    handlers=[
        RichHandler(show_path=False),
        logging.FileHandler("pipeline.log"),
    ],
)

console = Console()


def _parse_arguments():
    parser = argparse.ArgumentParser(description="GeneCodePro Splicing Analyzer")
    parser.add_argument("--input", "-i", type=Path)
    parser.add_argument("--output", "-o", type=Path)
    return parser.parse_args()


def _get_dataset_path(name_pattern: str) -> Path:
    for d in settings.DATASETS:
        if name_pattern in d["name"] or name_pattern in d["id"]:
            return settings.DATA_DIR / d["name"]
    raise FileNotFoundError(f"Dataset containing '{name_pattern}' not found.")


def main():
    args = _parse_arguments()

    if args.input:
        settings.DATA_DIR = args.input
    if args.output:
        settings.RESULTS_DIR = args.output

    console.rule("[bold magenta]GeneCode v1.0[/bold magenta]")

    with console.status("[bold green]Verificando arquivos...[/]"):
        download_datasets(settings.DATA_DIR)

    try:
        genomic_path = _get_dataset_path("genomic_ref")
        console.print(f"[dim]Carregando referência: {genomic_path.name}[/dim]")
        genomic_record = SeqIO.read(genomic_path, "fasta")

        analysis_map = [
            ("NM_000546", "Human TP53 (Variant 1)"),
            ("NM_001126112", "Human TP53 (Variant 2)"),
        ]

        analyzer = SplicingAnalyzer()
        reporter = SplicingReporter(settings.RESULTS_DIR)
        view = SplicingView()

        for mrna_pattern, label in analysis_map:
            mrna_path = _get_dataset_path(mrna_pattern)

            console.rule(f"[bold cyan]Analysing: {label}[/]")
            console.print(f"[dim]Genomic: {genomic_path.name} | mRNA: {mrna_path.name}[/dim]")

            mrna_record = SeqIO.read(mrna_path, "fasta")

            with console.status(f"[bold yellow]Mapping introns...[/]"):
                introns = analyzer.analyze(mrna_record.seq, genomic_record.seq)

            view.display_results(label, introns)

            output_path = reporter.generate_reports(label, introns)
            console.print(f"   [green]✔[/] Reports saved to: [underline]{output_path}[/]")

    except BioPipelineError as e:
        console.print(f"[red]Pipeline error:[/red] {e}")
        raise
    except Exception as e:
        console.print(f"[red]Unexpected error:[/red] {e}")
        raise

    console.rule("[bold green]Pipeline Finalizado![/bold green]")


if __name__ == "__main__":
    main()