import csv
import json
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional

from rich.console import Console

console = Console()


def _ts() -> str:
    return datetime.now().strftime("%Y%m%d_%H%M%S")


def export_json(data: List[Dict], output_dir: Path, prefix: str = "data") -> Path:
    path = output_dir / f"{prefix}_{_ts()}.json"
    with open(path, 'w') as f:
        json.dump(data, f, indent=2)
    console.print(f"  JSON:  {path.name}")
    return path


def export_csv(data: List[Dict], output_dir: Path, prefix: str = "data") -> Optional[Path]:
    if not data:
        return None
    path = output_dir / f"{prefix}_{_ts()}.csv"
    with open(path, 'w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=list(data[0].keys()))
        writer.writeheader()
        writer.writerows(data)
    console.print(f"  CSV:   {path.name}")
    return path


def generate_latex(
    title: str,
    summary_data: Dict[str, Any],
    headers: List[str],
    rows: List[List[str]],
    output_dir: Path,
    filename_prefix: str = "report",
) -> Path:
    path = output_dir / f"{filename_prefix}_{_ts()}.tex"
    summary_text = "\n".join(f"\\item \\textbf{{{k}}}: {v}" for k, v in summary_data.items())
    col_format = "l" * len(headers)
    header_row = " & ".join(f"\\textbf{{{h}}}" for h in headers)
    table_content = "".join(
        " & ".join(str(c).replace('_', r'\_').replace('%', r'\%') for c in row) + r" \\" + "\n"
        for row in rows
    )
    latex = rf"""\documentclass[11pt,a4paper]{{article}}
\usepackage[utf8]{{inputenc}}
\usepackage{{booktabs}}
\usepackage{{longtable}}
\usepackage{{geometry}}
\geometry{{margin=2.5cm}}

\title{{{title}}}
\author{{Bioinformatica Iniciante Pipeline}}
\date{{\today}}

\begin{{document}}
\maketitle

\section{{Summary}}
\begin{{itemize}}
{summary_text}
\end{{itemize}}

\section{{Detailed Results}}
\begin{{longtable}}{{{col_format}}}
\toprule
{header_row} \\
\midrule
{table_content}\bottomrule
\end{{longtable}}
\end{{document}}
"""
    path.write_text(latex, encoding='utf-8')
    return path

def export_results(results: List[Dict], output_dir: Path):
    output_dir.mkdir(parents=True, exist_ok=True)

    total = len(results)
    avg_gc = sum(r['gc_percent'] for r in results) / total if total else 0
    summary = {
        "Total Genomes Analyzed": total,
        "Average GC Content": f"{avg_gc:.2f}%",
        "Total ORFs Detected": sum(r['orf_count'] for r in results),
    }
    headers = ["Organism/ID", "Table", "Length", "GC%", "ORFs", "Top Codons"]
    rows = [
        [r['id'], str(r['table']), str(r['length']), f"{r['gc_percent']:.2f}%",
         str(r['orf_count']), ", ".join(x[0] for x in r['top_codons'][:3])]
        for r in results
    ]

    console.print("\n[green]Reports Generated:[/green]")
    export_json(results, output_dir, prefix="codon_usage")
    export_csv(results, output_dir, prefix="codon_usage")
    tex = generate_latex("GeneCodePro Analysis Report", summary, headers, rows, output_dir)
    console.print(f"  LaTeX: {tex.name}")