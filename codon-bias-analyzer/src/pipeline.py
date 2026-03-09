import argparse
import logging
from pathlib import Path

from Bio import SeqIO
from rich.console import Console
from rich.logging import RichHandler
from rich.progress import BarColumn, Progress, SpinnerColumn, TextColumn


def _init_logging() -> Console:
    logging.basicConfig(
        level="INFO",
        format="%(message)s",
        datefmt="[%X]",
        handlers=[RichHandler(show_path=False), logging.FileHandler("pipeline.log")],
    )
    return Console()


class BioPipeline:
    def __init__(self, name: str, version: str, config_obj):
        self.console = _init_logging()
        self.name = name
        self.version = version
        self.config = config_obj
        self.args = self._parse_arguments()
        self._setup_directories()

    def _parse_arguments(self):
        parser = argparse.ArgumentParser(description=f"{self.name} {self.version}")
        parser.add_argument("--input", "-i", type=Path)
        parser.add_argument("--output", "-o", type=Path)
        return parser.parse_args()

    def _setup_directories(self):
        if self.args.input:
            self.config.DATA_DIR = self.args.input
        if self.args.output:
            self.config.RESULTS_DIR = self.args.output

    def run_batch_analysis(self, process_func, pre_hook=None, post_hook=None):
        self.console.rule(f"[bold magenta]{self.name} {self.version}[/bold magenta]")

        if pre_hook:
            with self.console.status("[bold green]Running pre-flight checks...[/]"):
                pre_hook(self.config.DATA_DIR)

        files = [f for ext in self.config.FILE_EXTENSIONS for f in self.config.DATA_DIR.glob(ext)]
        results = []

        with Progress(
            SpinnerColumn(),
            TextColumn("[progress.description]{task.description}"),
            BarColumn(),
            TextColumn("[progress.percentage]{task.percentage:>3.0f}%"),
            console=self.console,
        ) as progress:
            task = progress.add_task("Analyzing...", total=len(files))
            for file_path in files:
                progress.update(task, description=f"Processing {file_path.name}")
                try:
                    record = SeqIO.read(file_path, "fasta")
                    result = process_func(record, file_path, self.console)
                    if result:
                        results.append(result)
                except Exception as e:
                    self.console.print(f"[red]Error {file_path.name}: {e}[/]")
                progress.advance(task)

        if post_hook and results:
            post_hook(results, self.config.RESULTS_DIR)

        self.console.rule("[bold green]Pipeline Complete[/bold green]")
