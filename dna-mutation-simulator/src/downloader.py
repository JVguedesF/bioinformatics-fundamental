import sys
from pathlib import Path
from typing import Dict, List

from Bio import Entrez

from src.config import settings

def download_sequences(email: str, datasets: List[Dict], output_dir: Path):
    if not email:
        print("\n[ERRO] Email do Entrez não fornecido!")
        print("Verifique seu arquivo .env e as configurações.\n")
        sys.exit(1)

    Entrez.email = email
    output_dir.mkdir(parents=True, exist_ok=True)

    for data in datasets:
        file_path = output_dir / data["name"]
        if file_path.exists():
            print(f"  - Já existe: {data['name']} (pulando)")
            continue

        print(f"  - Baixando {data['name']} ({data['id']})...")
        try:
            with Entrez.efetch(db="nucleotide", id=data["id"], rettype=data.get("type", "fasta"), retmode="text") as handle:
                file_path.write_text(handle.read())
            print(f"    Salvo: {data['name']}")
        except Exception as e:
            print(f"    [FALHA] Não foi possível baixar {data['id']}: {e}")

def download_datasets(output_dir: Path):
    download_sequences(
        email=settings.ENTREZ_EMAIL,
        datasets=settings.DATASETS,
        output_dir=output_dir,
    )