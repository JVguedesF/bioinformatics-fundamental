import sys
from pathlib import Path

ROOT_DIR = Path(__file__).resolve().parents[1]
if str(ROOT_DIR) not in sys.path:
    sys.path.insert(0, str(ROOT_DIR))

from src.analyzer import SequenceAnalyzer
from src.config import settings
from src.downloader import download_sequences
from src.pipeline import BioPipeline
from src.reporter import export_results
from src.view import CodonView


def process_file(record, file_path, console):
    dataset_info = next((d for d in settings.DATASETS if d["name"] == file_path.name), None)
    table_id = dataset_info.get("table", 1) if dataset_info else 1

    analyzer = SequenceAnalyzer(record.id, record.seq, table_id)
    orfs = analyzer.find_orfs(min_len_aa=settings.MIN_ORF_LENGTH)
    metrics = analyzer.analyze_codon_usage()

    header, tree = CodonView.create_analysis_view(
        filename=file_path.name,
        record_id=record.id,
        length=len(record.seq),
        table_id=table_id,
        orfs=orfs,
        metrics=metrics,
    )
    console.print(header)
    console.print(tree)
    console.print()

    return {
        "id": record.id,
        "file": file_path.name,
        "table": table_id,
        "length": len(record.seq),
        "orf_count": len(orfs),
        "gc_percent": metrics.gc_content,
        "top_codons": metrics.top_codons,
    }


def main():
    pipeline = BioPipeline("BioPipeline", "v1.0", settings)
    pipeline.run_batch_analysis(
        process_func=process_file,
        pre_hook=lambda out_dir: download_sequences(settings.ENTREZ_EMAIL, settings.DATASETS, out_dir),
        post_hook=export_results,
    )


if __name__ == "__main__":
    main()