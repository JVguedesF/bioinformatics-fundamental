from typing import List, Optional, Tuple, Union

from rich.console import Console, RenderableType
from rich.panel import Panel
from rich.table import Table
from rich.text import Text
from rich.tree import Tree

from src.config import settings
from src.models import AnalysisResult

console = Console()

def _format_sequence_preview(sequence: str, length: int = 60) -> Text:
    colors = {'A': 'green', 'T': 'red', 'U': 'red', 'G': 'yellow', 'C': 'blue', 'N': 'magenta'}
    text = Text()
    for char in sequence[:length].upper():
        text.append(char, style=colors.get(char, "white"))
    if len(sequence) > length:
        text.append("...", style="dim")
    return text


def _build_header_panel(filename: str, seq_id: str, seq_type: str, length: int, color: str) -> Panel:
    content = (
        f"[bold]FILE:[/bold] {filename}\n"
        f"[bold]ID:[/bold]   {seq_id}\n"
        f"[bold]TYPE:[/bold] [{color}]{seq_type}[/{color}] ({length} residues)"
    )
    return Panel(content, border_style=color, title="Sequence Information", expand=False)


def _build_tree(
    title: str,
    branches: List[Tuple[Optional[str], List[Union[str, RenderableType]]]]
) -> Tree:
    tree = Tree(f"[bold]{title}[/]", guide_style="dim")
    for branch_title, items in branches:
        target = tree.add(branch_title) if branch_title else tree
        for item in items:
            target.add(item)
    return tree

class StructuralView:

    @staticmethod
    def create_visualization(result: AnalysisResult, sequence_preview: str) -> Tuple[Panel, Tree]:
        color = settings.TYPE_COLORS.get(result.molecule_type, 'white')
        header = _build_header_panel(result.filename, result.sequence_id, result.molecule_type, result.length, color)

        branches: List[Tuple[Optional[str], List[Union[str, RenderableType]]]] = [
            (None, [Text("Sequence: ").append(_format_sequence_preview(sequence_preview))])
        ]

        if result.gc_content is not None:
            branches.append(("[blue]Genomic Analysis[/blue]", [
                f"GC Content: {result.gc_content:.2f}%",
                f"Melting Temperature: {result.melting_temp:.2f}°C",
            ]))

        if result.mfe is not None:
            items = [f"Minimum Free Energy: {result.mfe:.2f} kcal/mol"]
            if result.secondary_structure:
                items.append(f"Secondary Structure: {result.secondary_structure[:40]}...")
            branches.append(("[magenta]Transcriptomic Analysis[/magenta]", items))

        if result.molecular_weight is not None:
            items = [
                f"Molecular Weight: {result.molecular_weight:.2f} Da",
                f"Isoelectric Point: {result.isoelectric_point:.2f}",
            ]
            if result.stability_index:
                items.append(f"Instability Index: {result.stability_index}")
            branches.append(("[green]Proteomic Analysis[/green]", items))

        return header, _build_tree("Central Dogma Analysis", branches)

    @staticmethod
    def create_summary_table(results: List[AnalysisResult]) -> Table:
        table = Table(
            title="Comparative Analysis Summary",
            border_style="bright_black", header_style="bold cyan", show_lines=True,
        )
        for col, opts in [
            ("File",     {"style": "white", "no_wrap": True}),
            ("Type",     {"justify": "center"}),
            ("Length",   {"justify": "right"}),
            ("GC% / pI", {"justify": "right"}),
            ("Tm / MFE", {"justify": "right"}),
        ]:
            table.add_column(col, **opts)

        for res in results:
            def _fmt(v, spec, fallback="N/A"):
                return format(v, spec) if v is not None else fallback

            if res.molecule_type == 'DNA':
                m1 = _fmt(res.gc_content, ".1f", "N/A") + ("%" if res.gc_content is not None else "")
                m2 = _fmt(res.melting_temp, ".1f", "N/A") + ("°C" if res.melting_temp is not None else "")
                style = "blue"
            elif res.molecule_type == 'RNA':
                m1 = "N/A"
                m2 = (_fmt(res.mfe, ".1f") + " kcal/mol") if res.mfe is not None else "N/A"
                style = "magenta"
            else:
                m1 = _fmt(res.isoelectric_point, ".2f")
                m2 = (_fmt(res.molecular_weight, ".0f") + " Da") if res.molecular_weight is not None else "N/A"
                style = "green"
            table.add_row(res.filename, f"[{style}]{res.molecule_type}[/{style}]", str(res.length), m1, m2)
        return table