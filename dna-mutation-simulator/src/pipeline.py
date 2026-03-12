import argparse
from pathlib import Path
from Bio import SeqIO

def setup_arguments(name: str, version: str):
    parser = argparse.ArgumentParser(description=f"{name} {version}")
    parser.add_argument("--input", "-i", type=Path)
    parser.add_argument("--output", "-o", type=Path)
    return parser.parse_args()

def setup_directories(args, config):
    if args.input:
        config.DATA_DIR = args.input
    if args.output:
        config.RESULTS_DIR = args.output

def run_batch_analysis(name, version, config, process_func, pre_hook=None, post_hook=None):
    print("=" * 60)
    print(f"{name} {version}")
    print("=" * 60)

    if pre_hook:
        print("Running pre-flight checks...")
        pre_hook(config.DATA_DIR)

    files = [f for ext in config.FILE_EXTENSIONS for f in config.DATA_DIR.glob(ext)]
    results = []

    print(f"Found {len(files)} files to analyze.\n")

    for file_path in files:
        try:
            record = SeqIO.read(file_path, "fasta")
            result = process_func(record, file_path)
            if result:
                results.append(result)
        except Exception as e:
            print(f"\n[ERRO] {file_path.name}: {e}")

    if post_hook and results:
        post_hook(results, config.RESULTS_DIR)

    print("\n" + "=" * 60)
    print("Pipeline Complete")
    print("=" * 60)