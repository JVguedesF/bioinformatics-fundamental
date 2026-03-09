import sys
import matplotlib.pyplot as plt
from pathlib import Path
from rich.console import Console
from rich.panel import Panel
from rich.table import Table

ROOT_DIR = Path(__file__).resolve().parents[1]
if str(ROOT_DIR) not in sys.path:
    sys.path.insert(0, str(ROOT_DIR))

from src.analyser import ReplicationSimulator
from src.config import settings
from src.downloader import download_datasets
from src.model import DamageSimulator
from src.reporter import export_json
from src.view import (
    plot_replication_fork,
    plot_repair_efficiency,
    plot_mutational_spectrum,
    plot_recombination_hotspots,
    plot_meselson_stahl,
    plot_msi_analysis,
)

console = Console()


def verify_datasets(datasets_config: list, data_dir: Path) -> bool:
    for ds in datasets_config:
        if not (data_dir / ds["name"]).exists():
            console.print(f"[yellow]Arquivo ausente: {ds['name']}[/]")
            return True
    return False


def section(title: str):
    console.rule(f"[bold cyan]{title}[/bold cyan]")


def display_repair_results(repair_results: dict):
    table = Table(title="Resultados de Reparo", show_lines=True)
    table.add_column("Sistema", style="bold cyan")
    table.add_column("Dano simulado")
    table.add_column("Fidelidade ANTES")
    table.add_column("Fidelidade APÓS")
    table.add_column("Melhora")
    for sys_name, data in repair_results.items():
        before = data['before']['identity_score']
        after  = data['after']['identity_score']
        delta  = after - before
        color  = "green" if delta >= 0 else "red"
        table.add_row(sys_name, data['damage'], f"{before:.2f}%", f"[{color}]{after:.2f}%[/{color}]", f"[{color}]{delta:+.2f}%[/{color}]")
    console.print(table)


def run_msi_validation(sim: ReplicationSimulator, data_path: Path, results_dir: Path) -> dict:
    section("5 · Instabilidade de Microssatélites (MSI) — Validação Cruzada")

    msh2_seq = sim.load_genome(data_path / "Homo_sapiens_MSH2.fasta")
    normal_segment = msh2_seq[:5000]

    ms_normal = sim.find_microsatellites(normal_segment, min_repeats=4)
    console.print(f"  Microssatélites encontrados no MSH2: [cyan]{len(ms_normal)}[/]")

    msi_control = sim.detect_msi(normal_segment, normal_segment)
    console.print(f"  [bold]Controle (normal vs normal):[/] Loci: {msi_control['total_loci']} | Instáveis: {msi_control['instable_loci']} | Status: [green]{msi_control['status']}[/green]")

    tumor_segment = list(normal_segment)
    modified = 0
    for ms in ms_normal[:8]:
        pos, unit, n_rep = ms['position'], ms['unit'], ms['n_repeats']
        extra = unit * 2
        tumor_segment = tumor_segment[:pos + len(unit) * n_rep] + list(extra) + tumor_segment[pos + len(unit) * n_rep:]
        modified += 1
        if modified >= 5:
            break

    tumor_segment = "".join(tumor_segment)
    msi_result = sim.detect_msi(normal_segment, tumor_segment)
    console.print(f"  [bold]Tumor (normal vs tumor):[/] Loci: {msi_result['total_loci']} | Instáveis: [red]{msi_result['instable_loci']}[/red] | Status: [bold red]{msi_result['status']}[/bold red]")

    if msi_control['status'] == 'MSS' and 'MSI' in msi_result['status']:
        console.print("  [bold green]✓ Validação aprovada: controle MSS, tumor MSI detectado corretamente.[/]")
    else:
        console.print("  [bold yellow]⚠ Resultado inesperado — revisar parâmetros de detecção.[/]")

    plot_msi_analysis(msi_result, save_path=str(results_dir / "fig5_msi_analysis.png"))
    plt.show()
    return msi_result


def main():
    results_dir = ROOT_DIR / "results"
    data_path   = ROOT_DIR / "data" / "sequences"
    results_dir.mkdir(parents=True, exist_ok=True)
    data_path.mkdir(parents=True, exist_ok=True)

    section("0 · Verificação de Datasets")
    if verify_datasets(settings.DATASETS, data_path):
        console.print("[bold green]Iniciando downloads...[/]")
        try:
            download_datasets(data_path)
        except Exception as e:
            console.print(f"[red]Erro no download: {e}[/]")
            return
    else:
        console.print("[green]✓ Todos os datasets existem localmente.[/]")

    sim = ReplicationSimulator()
    ecoli_path = data_path / "Escherichia_coli_K12.fasta"
    yeast_path = data_path / "Saccharomyces_cerevisiae_Chr1.fasta"
    test_pattern = "TTATCCACA"

    section("1 · Replicação Bidirecional (E. coli — oriC)")
    forks = sim.start_replication(ecoli_path, test_pattern)
    template = forks['right']['template']
    leading  = forks['right']['leading']
    lagging  = forks['right']['lagging']
    console.print(f"  oriC encontrado em: [cyan]{forks['origin_position']}[/]")
    console.print(f"  Leading strand:     [cyan]{len(leading)} nt[/]")
    console.print(f"  Lagging strand:     [cyan]{len(lagging)} nt[/]")
    console.print(f"  Primer RNA (início):[cyan] {leading[:15]}[/]")
    plot_replication_fork(template, title="Forquilha de Replicação — E. coli", save_path=str(results_dir / "fig1_replication_fork.png"))
    plt.show()

    section("2 · Sistemas de Reparo (MMR / BER / NER / NHEJ / HR)")
    repair_results = sim.simulate_all_repair_systems(template, num_lesions=10)
    display_repair_results(repair_results)
    plot_repair_efficiency(repair_results, title="Eficiência dos Sistemas de Reparo", save_path=str(results_dir / "fig2_repair_efficiency.png"))
    plt.show()

    section("3 · Experimento de Meselson-Stahl (3 gerações)")
    ms_stats = sim.simulate_meselson_stahl(template[:200], generations=3)
    for gen_data in ms_stats:
        console.print(f"  Gen {gen_data['generation']}  |  N15/N15: [red]{gen_data['N15_N15_pct']:.0f}%[/]  N15/N14: [yellow]{gen_data['N15_N14_pct']:.0f}%[/]  N14/N14: [green]{gen_data['N14_N14_pct']:.0f}%[/]")
    plot_meselson_stahl(ms_stats, save_path=str(results_dir / "fig3_meselson_stahl.png"))
    plt.show()

    section("4 · Assinaturas Mutacionais COSMIC")
    dmg = DamageSimulator(seed=42)
    _, uv_events = dmg.apply_uv_damage(template, num_lesions=30)
    mutations = [(ev.position, ev.original_base, ev.damaged_base) for ev in uv_events]
    spectrum = sim.compute_mutational_spectrum(template, mutations)
    classified = sim.classify_cosmic_signature(spectrum)
    console.print(Panel(
        f"[bold]Substituição dominante:[/] {classified.get('dominant_substitution', 'N/A')}\n"
        f"[bold]Assinatura provável:[/]    {classified.get('top_signature', 'N/A')}\n"
        f"[bold]Etiologia:[/]             {classified.get('aetiology', 'N/A')}",
        title="Resultado COSMIC", border_style="yellow",
    ))
    plot_mutational_spectrum(spectrum, classified_sig=classified, title="Espectro Mutacional — Dano UV (simulado)", save_path=str(results_dir / "fig4_mutational_spectrum.png"))
    plt.show()

    msi_result = run_msi_validation(sim, data_path, results_dir)

    section("6 · Hotspots de Recombinação (PRDM9 + GC)")
    brca1_seq = sim.load_genome(data_path / "Homo_sapiens_BRCA1.fasta")
    hotspots = sim.find_recombination_hotspots(brca1_seq, window_size=500)
    console.print(f"  Sítios PRDM9 encontrados: [cyan]{hotspots['n_prdm9']}[/]  |  Janelas GC ≥ 55 %: [cyan]{hotspots['n_high_gc']}[/]  |  Densidade: [cyan]{hotspots['hotspot_density_per_mb']:.1f} / Mb[/]")
    plot_recombination_hotspots(hotspots, gc_windows=hotspots['gc_windows'], save_path=str(results_dir / "fig6_recombination_hotspots.png"))
    plt.show()

    section("7 · Recombinação Genética — Crossing-Over (S. cerevisiae)")
    new_c1, new_c2 = sim.simulate_recombination(yeast_path)
    console.print("  Ponto de crossing-over: 50000")
    console.print(f"  Cromossomo 1 [49995:50005]: [cyan]{new_c1[49995:50005]}[/]")
    console.print(f"  Cromossomo 2 [49995:50005]: [cyan]{new_c2[49995:50005]}[/]")

    section("8 · Empacotamento em Cromatina (Nucleossomos)")
    packed = sim.simulate_packing(new_c1)
    nucleosomes = [b for b in packed if b[0] == 'nucleosome']
    linkers     = [b for b in packed if b[0] == 'linker']
    console.print(f"  Total de blocos: {len(packed)}  |  Nucleossomos: {len(nucleosomes)}  |  Linkers: {len(linkers)}")
    for block in packed[:3]:
        console.print(f"  [{block[0]}] {len(block[1])} nt → {block[1][:12]}...")

    section("9 · Correlação Defeito de Reparo ↔ Câncer")
    for sys_name in ['MMR', 'BER', 'NER', 'NHEJ', 'HR']:
        info = sim.correlate_repair_defect_with_cancer(sys_name)
        console.print(Panel(
            f"[bold]Genes:[/]      {', '.join(info['genes'])}\n"
            f"[bold]Cânceres:[/]   {', '.join(info['cancers'])}\n"
            f"[bold]Assinatura:[/] {info['signature']}\n"
            f"[bold]Terapia:[/]    {info['therapy']}",
            title=f"[red]{sys_name}[/red] — {info['biomarker']}",
            border_style="red" if sys_name in ('MMR', 'HR') else "yellow",
            expand=False,
        ))

    section("10 · Exportação de Resultados")
    export_data = [
        {"system": k, "damage": v['damage'], "identity_before": round(v['before']['identity_score'], 4),
         "identity_after": round(v['after']['identity_score'], 4), "errors_before": v['before']['total_errors'], "errors_after": v['after']['total_errors']}
        for k, v in repair_results.items()
    ]
    export_json(export_data, results_dir, prefix="repair_results")
    export_json(ms_stats, results_dir, prefix="meselson_stahl")
    export_json([{"status": msi_result['status'], "total_loci": msi_result['total_loci'], "instable_loci": msi_result['instable_loci'], "msi_ratio": round(msi_result['msi_ratio'], 4)}], results_dir, prefix="msi_analysis")

    console.print(f"\n[bold green]✓ Pipeline completo! Resultados em:[/] {results_dir}")
    console.rule("[bold green]Projeto 5 — Concluído[/bold green]")


if __name__ == "__main__":
    main()