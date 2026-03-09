import sys
from pathlib import Path
from typing import Dict, List
from Bio import Entrez
from rich.console import Console

console = Console()


def download_sequences(email: str, datasets: List[Dict], output_dir: Path):
    """Baixa sequências do NCBI Entrez e salva em output_dir."""
    if not email:
        console.print("[red]Erro: Email do Entrez não fornecido![/]")
        console.print("[yellow]Verifique seu arquivo .env e as configurações.[/]")
        sys.exit(1)

    Entrez.email = email
    output_dir.mkdir(parents=True, exist_ok=True)

    for data in datasets:
        file_path = output_dir / data["name"]
        if file_path.exists():
            console.print(f"[cyan]Já existe: {data['name']} (pulando)[/]")
            continue

        console.print(f"[yellow]Downloading {data['name']} ({data['id']})...[/]")
        try:
            with Entrez.efetch(
                db="nucleotide",
                id=data["id"],
                rettype=data.get("type", "fasta"),
                retmode="text",
            ) as handle:
                file_path.write_text(handle.read())
            console.print(f"[green]Salvo: {data['name']}[/]")
        except Exception as e:
            console.print(f"[red]Falha ao baixar {data['id']}: {e}[/]")