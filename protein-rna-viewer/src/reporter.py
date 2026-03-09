import csv
import json
from dataclasses import asdict
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional

from rich.console import Console

from src.models import AnalysisResult

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


def export_results(results: List[AnalysisResult], output_dir: Path):
    """Ponto de entrada principal: exporta todos os formatos de uma vez."""
    output_dir.mkdir(exist_ok=True)

    counts = {t: sum(1 for r in results if r.molecule_type == t) for t in ['DNA', 'RNA', 'PROTEIN']}
    summary = {"Total Sequences": len(results), **{f"{k} Count": v for k, v in counts.items()}}

    headers = ["File", "ID", "Type", "Length", "GC% / pI", "Tm / MFE"]
    def _fmt(value, spec, fallback="N/A") -> str:
        return format(value, spec) if value is not None else fallback

    rows = []
    for res in results:
        if res.molecule_type == 'DNA':
            m1 = _fmt(res.gc_content, ".2f")
            m2 = _fmt(res.melting_temp, ".2f")
        elif res.molecule_type == 'RNA':
            m1 = "-"
            m2 = _fmt(res.mfe, ".2f")
        else:
            m1 = _fmt(res.isoelectric_point, ".2f")
            m2 = _fmt(res.molecular_weight, ".0f")
        rows.append([res.filename, res.sequence_id, res.molecule_type, str(res.length), m1, m2])

    console.print("\n[green]Results exported to:[/green]")
    data = [asdict(r) for r in results]
    export_json(data, output_dir)
    export_csv(data, output_dir)
    tex = generate_latex("Structural Profiling Report", summary, headers, rows, output_dir)
    console.print(f"  LaTeX: {tex.name}")