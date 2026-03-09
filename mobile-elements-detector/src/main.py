import sys
import argparse
from pathlib import Path

from rich.console import Console
from rich.table import Table
from rich.panel import Panel

ROOT_DIR = Path(__file__).resolve().parents[2]
if str(ROOT_DIR) not in sys.path:
    sys.path.insert(0, str(ROOT_DIR))

from .config import settings
from .reporter import BioReporter
from . import downloader, analyzer, view

console = Console()


def _has_missing_files(datasets, data_dir: Path) -> bool:
    if not data_dir.exists():
        return True

    for data in datasets:
        ext = "fasta" if data.get("group") == "nuclear_target" else "gb"
        if not (data_dir / f"{data['id']}.{ext}").exists():
            console.print(f"[yellow]Arquivo ausente: {data['id']}.{ext}[/yellow]")
            return True

    return False


def main():
    parser = argparse.ArgumentParser(description="Genomic Architecture Analyzer")
    parser.add_argument("--download", action="store_true", help="Force download of datasets")
    args = parser.parse_args()

    console.print(Panel.fit(
        "[bold cyan]Projeto 4: Endossimbiose & Elementos Móveis[/bold cyan]",
        subtitle="Bioinformatics Pipeline",
    ))

    if args.download or _has_missing_files(settings.DATASETS, settings.DATA_DIR):
        console.print("[bold yellow]Iniciando download dos datasets...[/bold yellow]")
        downloader.download_datasets(settings.DATASETS, settings.DATA_DIR)
        console.print("[bold green]✓ Downloads concluídos![/bold green]")
    else:
        console.print("[bold green]✓ Todos os arquivos verificados localmente.[/bold green]")

    console.print("\n[bold yellow]--- 1. Análise de Densidade Gênica ---[/bold yellow]")

    with console.status("[bold white]Processando arquivos..."):
        density_data = analyzer.pipeline_density_analysis(settings.DATASETS, settings.DATA_DIR)

    if not density_data:
        console.print("[bold red]Nenhum dado processado. Verifique os arquivos .gb em data/.[/bold red]")
    else:
        table = Table(title="Comparação de Densidade")
        table.add_column("Genoma", style="cyan")
        table.add_column("Grupo", style="magenta")
        table.add_column("Genes/kb", justify="right", style="green")

        for res in density_data:
            table.add_row(res.genome_name, res.group, f"{res.density_genes_per_kb:.2f}")

        console.print(table)
        console.print("[bold white]Exportando resultados...[/bold white]")
        BioReporter.export_density(density_data, settings.RESULTS_DIR)

        if console.input("\n[dim]Gerar gráfico de densidade? (s/n): [/dim]").lower() == "s":
            path = view.plot_density_comparison(density_data, settings.RESULTS_DIR)
            console.print(f"[dim]  ↳ gráfico  → {path.name}[/dim]")

    console.print("\n[bold yellow]--- 2. Análise de GC Skew (Replicação) ---[/bold yellow]")

    with console.status("[bold white]Calculando Skew..."):
        skew_data = analyzer.pipeline_skew_analysis(
            settings.DATASETS,
            settings.DATA_DIR,
            settings.GC_SKEW_WINDOW,
            settings.GC_SKEW_STEP,
        )

    if not skew_data:
        console.print("[bold red]Nenhum dado circular encontrado para análise de Skew.[/bold red]")
    else:
        for res in skew_data:
            console.print(
                f" > [cyan]{res.genome_name}[/cyan]: "
                f"Origem prevista na base [bold red]{res.predicted_ori_pos}[/bold red]"
            )
        console.print("[bold white]Exportando resultados...[/bold white]")
        BioReporter.export_skew(skew_data, settings.RESULTS_DIR)

        if console.input("\n[dim]Gerar gráficos de GC Skew? (s/n): [/dim]").lower() == "s":
            for res in skew_data:
                path = view.plot_gc_skew_curve(res, settings.RESULTS_DIR)
                console.print(f"[dim]  ↳ gráfico  → {path.name}[/dim]")

    console.print("\n[bold yellow]--- 3. Rastreamento de Elementos Móveis (Alu) ---[/bold yellow]")

    with console.status("[bold white]Escaneando Genoma Nuclear..."):
        mobile_data = analyzer.pipeline_mobile_elements(
            settings.DATASETS,
            settings.DATA_DIR,
            settings.ALU_CONSENSUS_SEQ,
        )

    if mobile_data:
        console.print(Panel(
            f"Alvo: [bold]{mobile_data.genome_name}[/bold]\n"
            f"Elemento Buscado: Alu Consensus\n"
            f"Ocorrências Encontradas: [bold red]{mobile_data.total_matches}[/bold red]",
            title="Resultado do Scan",
            border_style="green",
        ))
        console.print("[bold white]Exportando resultados...[/bold white]")
        BioReporter.export_mobile_elements(mobile_data, settings.RESULTS_DIR)

        if console.input("\n[dim]Gerar mapa cromossômico? (s/n): [/dim]").lower() == "s":
            path = view.plot_chromosome_map(mobile_data, settings.RESULTS_DIR)
            console.print(f"[dim]  ↳ gráfico  → {path.name}[/dim]")
    else:
        console.print("[red]Aviso: Dataset nuclear não encontrado ou erro na leitura.[/red]")

    console.print("\n[bold green]Pipeline Finalizado![/bold green]")


if __name__ == "__main__":
    main()