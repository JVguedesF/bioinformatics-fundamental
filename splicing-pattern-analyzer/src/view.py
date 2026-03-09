from typing import List

from rich.console import Console
from rich.panel import Panel
from rich.table import Table

from .models import Intron


class SplicingView:
    def __init__(self):
        self.console = Console()

    def display_gene_map(self, total_len: int, introns: List[Intron]) -> None:
        if not introns:
            return

        width = 60
        scale = width / total_len
        visual_map = ""
        current_pos = 0

        for intron in introns:
            exon_chars = max(1, int((intron.start - current_pos) * scale))
            visual_map += "[green]" + ("█" * exon_chars) + "[/]"

            intron_chars = max(1, int(intron.length * scale))
            visual_map += "[dim white]" + ("░" * intron_chars) + "[/]"

            current_pos = intron.end

        self.console.print(Panel(visual_map, title="[b]Gene Structure Map[/]", border_style="cyan"))

    def display_results(self, gene_name: str, introns: List[Intron]) -> None:
        if not introns:
            self.console.print(f"[yellow]Nenhum íntron encontrado para {gene_name}.[/]")
            return

        total_len = introns[-1].end
        avg_gc = sum(i.gc_content for i in introns) / len(introns)
        canonical_count = sum(1 for i in introns if i.is_canonical)

        summary = (
            f"\n        [bold]Alvo:[/bold] {gene_name}\n"
            f"        [bold]Total Íntrons:[/bold] {len(introns)}\n"
            f"        [bold]Padrão Canônico (GT-AG):[/bold] {canonical_count}/{len(introns)}\n"
            f"        [bold]Conteúdo GC Médio:[/bold] {avg_gc:.1f}%\n"
        )
        self.console.print(Panel(summary, title="🧬 Relatório de Splicing", border_style="magenta"))

        self.display_gene_map(total_len, introns)
        self.console.print()

        table = Table(show_header=True, header_style="bold white on blue", expand=True)
        table.add_column("ID", justify="center", style="dim")
        table.add_column("Localização", justify="center")
        table.add_column("Tamanho", justify="right")
        table.add_column("Doador...Aceitador", justify="center")
        table.add_column("GC Content", justify="right")
        table.add_column("Status", justify="center")

        for intron in introns:
            gc_color = "green" if intron.gc_content > 50 else "yellow"
            status = "[green]✔ OK[/]" if intron.is_canonical else "[red bold]⚠ ATENÇÃO[/]"
            sites_str = f"[bold white]{intron.donor_site}[/]...[bold white]{intron.acceptor_site}[/]"

            table.add_row(
                str(intron.id),
                f"{intron.start:,} - {intron.end:,}",
                f"{intron.length:,} bp",
                sites_str,
                f"[{gc_color}]{intron.gc_content:.1f}%[/]",
                status,
            )

        self.console.print(table)
        self.console.print()