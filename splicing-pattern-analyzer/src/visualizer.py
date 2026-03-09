from typing import List, Tuple, Union, Optional
from rich.panel import Panel
from rich.tree import Tree
from rich.text import Text
from rich.console import RenderableType

class BioVisualizer:

    @staticmethod
    def format_sequence_preview(sequence: str, length: int = 60) -> Text:
        text = Text()
        colors = {
            'A': 'green', 'T': 'red', 'U': 'red',
            'G': 'yellow', 'C': 'blue',
            'N': 'magenta'
        }
        
        preview = sequence[:length].upper()
        for char in preview:
            color = colors.get(char, "white")
            text.append(char, style=color)
            
        if len(sequence) > length:
            text.append("...", style="dim")
            
        return text

    @staticmethod
    def create_header_panel(
        filename: str,
        seq_id: str,
        seq_type: str,
        length: int,
        color: str = "white"
    ) -> Panel:
        content = (
            f"[bold]FILE:[/bold] {filename}\n"
            f"[bold]ID:[/bold]   {seq_id}\n"
            f"[bold]TYPE:[/bold] [{color}]{seq_type}[/{color}] ({length} residues)"
        )
        return Panel(
            content,
            border_style=color,
            title="Sequence Information",
            expand=False
        )

    @staticmethod
    def create_result_tree(
        title: str,
        branches: List[Tuple[Optional[str], List[Union[str, RenderableType]]]]
    ) -> Tree:
        tree = Tree(f"[bold]{title}[/]", guide_style="dim")
        
        for branch_title, items in branches:
            target = tree.add(branch_title) if branch_title else tree
            for item in items:
                target.add(item)
                    
        return tree