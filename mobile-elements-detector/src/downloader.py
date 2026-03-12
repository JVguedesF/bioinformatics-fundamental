import sys
from pathlib import Path
from typing import List, Dict
from Bio import Entrez
from src.config import settings

def download_datasets(datasets: List[Dict], output_dir: Path) -> None:
    if not settings.ENTREZ_EMAIL:
        print("  [Erro] Email do Entrez não fornecido!")
        sys.exit(1)
    Entrez.email = settings.ENTREZ_EMAIL
    output_dir.mkdir(parents=True, exist_ok=True)
    for data in datasets:
        ext = "fasta" if data.get("group") == "nuclear_target" else "gb"
        file_path = output_dir / f"{data['id']}.{ext}"
        if file_path.exists(): continue
        print(f"  - Baixando {data['id']}...")
        r_type = "gbwithparts" if data.get("type") == "genbank" else "fasta"
        try:
            with Entrez.efetch(db="nucleotide", id=data["id"], rettype=r_type, retmode="text") as handle:
                file_path.write_text(handle.read())
            print(f"    ✓ {file_path.name}")
        except Exception as e:
            print(f"    [Falha] {data['id']}: {e}")