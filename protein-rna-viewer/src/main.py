import sys
from pathlib import Path

ROOT_DIR = Path(__file__).resolve().parents[1]
if str(ROOT_DIR) not in sys.path:
    sys.path.insert(0, str(ROOT_DIR))

from src.analyzer import CentralDogmaAnalyzer
from src.config import settings
from src.pipeline import BioPipeline
from src.reporter import export_results
from src.view import StructuralView


def process_file(record, file_path, console):
    result = CentralDogmaAnalyzer().process_sequence(record, file_path.name)
    header, tree = StructuralView.create_visualization(result, str(record.seq))
    console.print(header)
    console.print(tree)
    console.print()
    return result


def main():
    pipeline = BioPipeline("BioPipeline", "v1.0", settings)
    pipeline.run_batch_analysis(
        process_func=process_file,
        post_hook=lambda results, out_dir: (
            print(StructuralView.create_summary_table(results)),
            export_results(results, out_dir),
        ),
    )


if __name__ == "__main__":
    main()