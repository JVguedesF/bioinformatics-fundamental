import sys
from pathlib import Path
from typing import List, Dict

from Bio import Entrez
from rich.console import Console
from rich.progress import track

from .config import settings

console = Console()


def download_datasets(datasets: List[Dict], output_dir: Path) -> None:
    if not settings.ENTREZ_EMAIL:
        console.print("[red]Erro: Email do Entrez não fornecido![/]")
        console.print("[yellow]Verifique seu arquivo .env e as configurações.[/]")
        sys.exit(1)

    Entrez.email = settings.ENTREZ_EMAIL
    output_dir.mkdir(parents=True, exist_ok=True)

    for data in track(datasets, description="[green]Baixando arquivos...[/]"):
        ext = "fasta" if data.get("group") == "nuclear_target" else "gb"
        file_path = output_dir / f"{data['id']}.{ext}"

        if file_path.exists():
            console.print(f"[cyan]Já existe: {file_path.name} (pulando)[/]")
            continue

        r_type = "gbwithparts" if data.get("type") == "genbank" else "fasta"

        try:
            with Entrez.efetch(db="nucleotide", id=data["id"], rettype=r_type, retmode="text") as handle:
                file_path.write_text(handle.read())
            console.print(f"[green]✓ {file_path.name}[/]")
        except Exception as e:
            console.print(f"[bold red]Falha ao baixar {data['id']}: {e}[/]")