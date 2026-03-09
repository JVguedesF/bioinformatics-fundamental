from typing import List, Optional, Tuple, Union

from rich.console import RenderableType
from rich.panel import Panel
from rich.tree import Tree


class BioVisualizer:

    @staticmethod
    def create_header_panel(
        filename: str,
        seq_id: str,
        seq_type: str,
        length: int,
        color: str = "white",
    ) -> Panel:
        content = (
            f"[bold]FILE:[/bold] {filename}\n"
            f"[bold]ID:[/bold]   {seq_id}\n"
            f"[bold]TYPE:[/bold] [{color}]{seq_type}[/{color}] ({length} residues)"
        )
        return Panel(content, border_style=color, title="Sequence Information", expand=False)

    @staticmethod
    def create_result_tree(
        title: str,
        branches: List[Tuple[Optional[str], List[Union[str, RenderableType]]]],
    ) -> Tree:
        tree = Tree(f"[bold]{title}[/]", guide_style="dim")
        for branch_title, items in branches:
            target = tree.add(branch_title) if branch_title else tree
            for item in items:
                target.add(item)
        return tree
