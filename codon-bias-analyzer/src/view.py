from typing import List, Tuple, Union

from rich.console import Console, RenderableType
from rich.panel import Panel
from rich.tree import Tree

from src.models import CodonMetrics, ORFResult
from src.visualizer import BioVisualizer

console = Console()


class CodonView:

    @staticmethod
    def create_analysis_view(
        filename: str,
        record_id: str,
        length: int,
        table_id: int,
        orfs: List[ORFResult],
        metrics: CodonMetrics,
    ) -> Tuple[Panel, Tree]:
        org_type = "BACTERIA (Table 11)" if table_id == 11 else "STANDARD (Table 1)"
        color = "cyan" if table_id == 11 else "green"

        header = BioVisualizer.create_header_panel(filename, record_id, org_type, length, color)
        branches = []

        gc_pct = (metrics.gc_content / (metrics.total_codons * 3)) * 100 if metrics.total_codons > 0 else 0
        branches.append((f"[{color}]Genomic Metrics[/{color}]", [
            f"GC Content: {gc_pct:.2f}%",
            f"Total Codons: {metrics.total_codons}",
        ]))

        orf_items: List[Union[str, RenderableType]] = [f"Total ORFs detected: {len(orfs)} (>100aa)"]
        if orfs:
            best = orfs[0]
            best_tree = Tree(f"[bold]Longest Candidate:[/bold] {best.length_aa} aa")
            best_tree.add(f"Frame: {best.frame} | Strand: {best.strand}")
            best_tree.add(f"Preview: [dim]{best.protein_seq[:25]}...[/dim]")
            orf_items.append(best_tree)
        branches.append(("[yellow]ORF Prediction[/yellow]", orf_items))

        top_codons = [f"[bold]{c}[/bold]: {n} occurrences" for c, n in metrics.top_codons[:3]]
        branches.append(("[magenta]Codon Bias Profile[/magenta]", top_codons))

        return header, BioVisualizer.create_result_tree("Genetic Code & ORF Analysis", branches)