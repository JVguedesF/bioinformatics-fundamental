import sys
from pathlib import Path

from Bio import Entrez
from rich.console import Console
from rich.progress import track

from .config import settings

console = Console()


def download_datasets(output_dir: Path) -> None:
    if not settings.ENTREZ_EMAIL:
        console.print("[red]Erro: Email do Entrez não fornecido![/]")
        console.print("[yellow]Verifique seu arquivo .env e as configurações.[/]")
        sys.exit(1)

    Entrez.email = settings.ENTREZ_EMAIL
    output_dir.mkdir(parents=True, exist_ok=True)

    for data in track(settings.DATASETS, description="[green]Verificando arquivos...[/]"):
        file_path = output_dir / data["name"]

        if file_path.exists():
            console.print(f"[cyan]Já existe: {data['name']} (pulando)[/]")
            continue

        try:
            with Entrez.efetch(
                db="nucleotide",
                id=data["id"],
                rettype=data.get("type", "fasta"),
                retmode="text",
            ) as handle:
                file_path.write_text(handle.read())
            console.print(f"[green]✓ {data['name']}[/]")
        except Exception as e:
            console.print(f"[bold red]Falha ao baixar {data['id']}: {e}[/]")