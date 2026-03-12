import sys
import matplotlib.pyplot as plt
from pathlib import Path

ROOT_DIR = Path(__file__).resolve().parents[1]
if str(ROOT_DIR) not in sys.path:
    sys.path.insert(0, str(ROOT_DIR))

from src.analyser import ReplicationSimulator
from src.config import settings
from src.downloader import download_datasets
from src.model import DamageSimulator
from src.reporter import export_json
from src.printer import TerminalPrinter
from src.view import (
    plot_replication_fork,
    plot_repair_efficiency,
    plot_mutational_spectrum,
    plot_recombination_hotspots,
    plot_meselson_stahl,
    plot_msi_analysis,
)


def verify_datasets(datasets_config: list, data_dir: Path) -> bool:
    for ds in datasets_config:
        if not (data_dir / ds["name"]).exists():
            TerminalPrinter.dataset_status(ds['name'])
            return True
    return False


def run_msi_validation(sim: ReplicationSimulator, data_path: Path, graphs_dir: Path) -> dict:
    TerminalPrinter.section("5 · Instabilidade de Microssatélites (MSI) — Validação Cruzada")

    msh2_seq = sim.load_genome(data_path / "Homo_sapiens_MSH2.fasta")
    normal_segment = msh2_seq[:5000]

    ms_normal = sim.find_microsatellites(normal_segment, min_repeats=4)
    msi_control = sim.detect_msi(normal_segment, normal_segment)

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

    TerminalPrinter.msi_validation(len(ms_normal), msi_control, msi_result)

    plot_msi_analysis(msi_result, save_path=str(graphs_dir / "fig5_msi_analysis.png"))
    plt.show()
    return msi_result


def main():
    results_dir = ROOT_DIR / "results"
    graphs_dir = results_dir / "graphs"
    data_path = ROOT_DIR / "data" / "sequences"

    results_dir.mkdir(parents=True, exist_ok=True)
    graphs_dir.mkdir(parents=True, exist_ok=True)
    data_path.mkdir(parents=True, exist_ok=True)

    TerminalPrinter.header("DNA Mutation Simulator", "v1.0")

    TerminalPrinter.section("0 · Verificação de Datasets")
    if verify_datasets(settings.DATASETS, data_path):
        print("  A iniciar downloads...")
        try:
            download_datasets(data_path)
        except Exception as e:
            print(f"  [ERRO] Falha no download: {e}")
            return
    else:
        TerminalPrinter.dataset_status(None)

    sim = ReplicationSimulator()
    ecoli_path = data_path / "Escherichia_coli_K12.fasta"
    yeast_path = data_path / "Saccharomyces_cerevisiae_Chr1.fasta"
    test_pattern = "TTATCCACA"

    TerminalPrinter.section("1 · Replicação Bidirecional (E. coli — oriC)")
    forks = sim.start_replication(ecoli_path, test_pattern)
    template = forks['right']['template']
    leading = forks['right']['leading']
    lagging = forks['right']['lagging']
    TerminalPrinter.replication_fork(forks['origin_position'], len(leading), len(lagging), leading[:15])
    plot_replication_fork(template, title="Forquilha de Replicação — E. coli",
                          save_path=str(graphs_dir / "fig1_replication_fork.png"))
    plt.show()

    TerminalPrinter.section("2 · Sistemas de Reparo (MMR / BER / NER / NHEJ / HR)")
    repair_results = sim.simulate_all_repair_systems(template, num_lesions=10)
    TerminalPrinter.repair_results(repair_results)
    plot_repair_efficiency(repair_results, title="Eficiência dos Sistemas de Reparo",
                           save_path=str(graphs_dir / "fig2_repair_efficiency.png"))
    plt.show()

    TerminalPrinter.section("3 · Experimento de Meselson-Stahl (3 gerações)")
    ms_stats = sim.simulate_meselson_stahl(template[:200], generations=3)
    TerminalPrinter.meselson_stahl(ms_stats)
    plot_meselson_stahl(ms_stats, save_path=str(graphs_dir / "fig3_meselson_stahl.png"))
    plt.show()

    TerminalPrinter.section("4 · Assinaturas Mutacionais COSMIC")
    dmg = DamageSimulator(seed=42)
    _, uv_events = dmg.apply_uv_damage(template, num_lesions=30)
    mutations = [(ev.position, ev.original_base, ev.damaged_base) for ev in uv_events]
    spectrum = sim.compute_mutational_spectrum(template, mutations)
    classified = sim.classify_cosmic_signature(spectrum)
    TerminalPrinter.cosmic_signature(classified)
    plot_mutational_spectrum(spectrum, classified_sig=classified, title="Espectro Mutacional — Dano UV (simulado)",
                             save_path=str(graphs_dir / "fig4_mutational_spectrum.png"))
    plt.show()

    msi_result = run_msi_validation(sim, data_path, graphs_dir)

    TerminalPrinter.section("6 · Hotspots de Recombinação (PRDM9 + GC)")
    brca1_seq = sim.load_genome(data_path / "Homo_sapiens_BRCA1.fasta")
    hotspots = sim.find_recombination_hotspots(brca1_seq, window_size=500)
    TerminalPrinter.recombination_hotspots(hotspots)
    plot_recombination_hotspots(hotspots, gc_windows=hotspots['gc_windows'],
                                save_path=str(graphs_dir / "fig6_recombination_hotspots.png"))
    plt.show()

    TerminalPrinter.section("7 · Recombinação Genética — Crossing-Over (S. cerevisiae)")
    new_c1, new_c2 = sim.simulate_recombination(yeast_path)
    TerminalPrinter.crossing_over(new_c1[49995:50005], new_c2[49995:50005])

    TerminalPrinter.section("8 · Empacotamento em Cromatina (Nucleossomos)")
    packed = sim.simulate_packing(new_c1)
    nucleosomes = [b for b in packed if b[0] == 'nucleosome']
    linkers = [b for b in packed if b[0] == 'linker']
    TerminalPrinter.chromatin_packing(packed, nucleosomes, linkers)

    TerminalPrinter.section("9 · Correlação Defeito de Reparo ↔ Câncer")
    for sys_name in ['MMR', 'BER', 'NER', 'NHEJ', 'HR']:
        info = sim.correlate_repair_defect_with_cancer(sys_name)
        TerminalPrinter.repair_defects(sys_name, info)

    TerminalPrinter.section("10 · Exportação de Resultados")
    export_data = [
        {"system": k, "damage": v['damage'], "identity_before": round(v['before']['identity_score'], 4),
         "identity_after": round(v['after']['identity_score'], 4), "errors_before": v['before']['total_errors'],
         "errors_after": v['after']['total_errors']}
        for k, v in repair_results.items()
    ]
    export_json(export_data, results_dir, prefix="repair_results")
    export_json(ms_stats, results_dir, prefix="meselson_stahl")
    export_json([{"status": msi_result['status'], "total_loci": msi_result['total_loci'],
                  "instable_loci": msi_result['instable_loci'], "msi_ratio": round(msi_result['msi_ratio'], 4)}],
                results_dir, prefix="msi_analysis")

    TerminalPrinter.footer(results_dir)


if __name__ == "__main__":
    main()