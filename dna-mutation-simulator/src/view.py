from __future__ import annotations
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
from typing import Optional


COLORS = {
    'C>A': '#1e90ff',
    'C>G': '#222222',
    'C>T': '#e03030',
    'T>A': '#aaaaaa',
    'T>C': '#66cc44',
    'T>G': '#ff9933',
}
REPAIR_COLORS = {
    'MMR':  '#4e79a7',
    'BER':  '#f28e2b',
    'NER':  '#59a14f',
    'NHEJ': '#e15759',
    'HR':   '#76b7b2',
}


def plot_replication_fork(
    template: str,
    okazaki_size: int = 1000,
    display_length: int = 5000,
    title: str = "Forquilha de Replicação",
    save_path: Optional[str] = None,
) -> plt.Figure:
    fig, ax = plt.subplots(figsize=(14, 5))
    ax.set_xlim(0, display_length)
    ax.set_ylim(-3, 3)
    ax.axis('off')
    ax.set_title(title, fontsize=14, fontweight='bold', pad=15)

    length = min(display_length, len(template))

    ax.annotate(
        '', xy=(length, 1.5), xytext=(0, 1.5),
        arrowprops={'arrowstyle': '->', 'color': '#555555', 'lw': 2}
    )
    ax.text(-50, 1.5, "5'", va='center', ha='right', fontsize=9, color='#555555')
    ax.text(length + 30, 1.5, "3'", va='center', ha='left', fontsize=9, color='#555555')
    ax.text(length / 2, 1.85, f"Fita Molde (template)  [{length} nt]",
            ha='center', va='bottom', fontsize=9, color='#555555')

    ax.annotate(
        '', xy=(0, 0), xytext=(length, 0),
        arrowprops={'arrowstyle': '->', 'color': '#2196f3', 'lw': 2.5}
    )
    ax.text(length + 30, 0, "5'", va='center', ha='left', fontsize=9, color='#2196f3')
    ax.text(-50, 0, "3'", va='center', ha='right', fontsize=9, color='#2196f3')
    ax.text(length / 2, 0.35, "Fita Leading (contínua)", ha='center', fontsize=9, color='#2196f3')

    frag_colors = ['#e91e63', '#9c27b0']
    n_frags = length // okazaki_size
    for i in range(n_frags):
        start = i * okazaki_size
        end = min(start + okazaki_size, length)
        color = frag_colors[i % 2]
        ax.annotate(
            '', xy=(start + 5, -1.5), xytext=(end - 5, -1.5),
            arrowprops={'arrowstyle': '->', 'color': color, 'lw': 2}
        )
        ax.text((start + end) / 2, -1.15, f"Frag {i+1}", ha='center', fontsize=7, color=color)

    ax.text(-50, -1.5, "3'", va='center', ha='right', fontsize=9, color='#9c27b0')
    ax.text(length + 30, -1.5, "5'", va='center', ha='left', fontsize=9, color='#9c27b0')
    ax.text(length / 2, -2.0,
            f"Fita Lagging — {n_frags} fragmentos de Okazaki (~{okazaki_size} nt cada)",
            ha='center', fontsize=9, color='#9c27b0')

    for i in range(n_frags):
        primer_start = i * okazaki_size
        primer_end = primer_start + 10
        ax.axvspan(primer_start, primer_end, alpha=0.3, color='#ff9800', ymin=0.1, ymax=0.45)

    ax.axvline(0, color='red', lw=2, linestyle='--', alpha=0.6)
    ax.text(5, 2.4, "oriC", fontsize=9, color='red', fontweight='bold')

    legend_elements = [
        mpatches.Patch(color='#555555', label='Fita Molde'),
        mpatches.Patch(color='#2196f3', label='Leading (contínua)'),
        mpatches.Patch(color='#9c27b0', label='Lagging (Okazaki)'),
        mpatches.Patch(color='#ff9800', alpha=0.5, label='RNA Primer (Primase)'),
    ]
    ax.legend(handles=legend_elements, loc='lower right', fontsize=8, framealpha=0.9)

    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, dpi=150, bbox_inches='tight')
    return fig


def plot_repair_efficiency(
    repair_results: dict[str, dict],
    title: str = "Eficiência dos Sistemas de Reparo de DNA",
    save_path: Optional[str] = None,
) -> plt.Figure:
    systems = list(repair_results.keys())
    before = [repair_results[s]['before']['identity_score'] for s in systems]
    after  = [repair_results[s]['after']['identity_score']  for s in systems]

    x = np.arange(len(systems))
    width = 0.35

    fig, ax = plt.subplots(figsize=(10, 6))

    bars_before = ax.bar(x - width/2, before, width, label='Antes do reparo',
                         color='#ff9999', edgecolor='white')
    bars_after  = ax.bar(x + width/2, after,  width, label='Após reparo',
                         color=[REPAIR_COLORS.get(s, '#66b3ff') for s in systems],
                         edgecolor='white')

    for bar in bars_before:
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.1,
                f"{bar.get_height():.2f}%", ha='center', va='bottom', fontsize=8)
    for bar in bars_after:
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.1,
                f"{bar.get_height():.2f}%", ha='center', va='bottom', fontsize=8)

    for i, s in enumerate(systems):
        dmg = repair_results[s].get('damage', '')
        ax.text(i, 0.5, dmg, ha='center', va='bottom', fontsize=7,
                color='#555555', style='italic', rotation=0)

    ax.set_xticks(x)
    ax.set_xticklabels(systems, fontsize=11, fontweight='bold')
    ax.set_ylabel("Fidelidade Genômica (%)", fontsize=11)
    ax.set_title(title, fontsize=13, fontweight='bold')
    ax.set_ylim(0, 105)
    ax.legend(fontsize=10)
    ax.grid(axis='y', alpha=0.3)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, dpi=150, bbox_inches='tight')
    return fig


def plot_mutational_spectrum(
    spectrum: dict[str, int],
    title: str = "Espectro Mutacional — 96 Contextos Trinucleotídicos (COSMIC)",
    classified_sig: Optional[dict] = None,
    save_path: Optional[str] = None,
) -> plt.Figure:
    if not spectrum:
        fig, ax = plt.subplots()
        ax.text(0.5, 0.5, "Espectro vazio", ha='center', va='center')
        return fig

    sub_types = ['C>A', 'C>G', 'C>T', 'T>A', 'T>C', 'T>G']
    bases = ['A', 'C', 'G', 'T']

    labels, values, bar_colors = [], [], []
    for sub in sub_types:
        _, _ = sub[0], sub[2]
        for b5 in bases:
            for b3 in bases:
                key = f"{b5}[{sub}]{b3}"
                labels.append(key)
                values.append(spectrum.get(key, 0))
                bar_colors.append(COLORS[sub])

    fig, ax = plt.subplots(figsize=(18, 5))
    x = np.arange(len(labels))
    ax.bar(x, values, color=bar_colors, width=0.8, edgecolor='none')

    for i, sub in enumerate(sub_types):
        start = i * 16
        end = start + 16
        ax.axvspan(start - 0.5, end - 0.5, alpha=0.08,
                   color=COLORS[sub], zorder=0)
        ax.text((start + end) / 2, max(values) * 1.05 if values else 1,
                sub, ha='center', fontsize=11, fontweight='bold',
                color=COLORS[sub])

    ax.set_xticks(x[::4])
    ax.set_xticklabels(labels[::4], rotation=90, fontsize=6)
    ax.set_ylabel("Número de Mutações", fontsize=10)
    ax.set_xlim(-1, 96)
    ax.set_title(title, fontsize=12, fontweight='bold')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    if classified_sig:
        sig = classified_sig.get('top_signature', '')
        aet = classified_sig.get('aetiology', '')
        ax.text(0.98, 0.95, f"Assinatura: {sig}\n{aet}",
                transform=ax.transAxes, ha='right', va='top',
                fontsize=8, style='italic',
                bbox={'boxstyle': 'round,pad=0.3', 'facecolor': '#fffde7', 'alpha': 0.8})

    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, dpi=150, bbox_inches='tight')
    return fig


def plot_recombination_hotspots(
    hotspot_data: dict,
    gc_windows: list[dict],
    title: str = "Hotspots de Recombinação",
    save_path: Optional[str] = None,
) -> plt.Figure:
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(14, 7), sharex=True)

    if gc_windows:
        positions = [w['start'] for w in gc_windows]
        gc_values = [w['gc_pct'] for w in gc_windows]
        ax1.plot(positions, gc_values, color='#2196f3', lw=1.2, label='GC %')
        ax1.axhline(55, color='orange', lw=1.5, linestyle='--', label='Limiar 55 %')
        high_gc = [(p, g) for p, g in zip(positions, gc_values) if g >= 55]
        if high_gc:
            hx, hy = zip(*high_gc)
            ax1.scatter(hx, hy, color='orange', s=15, zorder=3, label='GC ≥ 55 %')
    ax1.set_ylabel("GC Content (%)", fontsize=10)
    ax1.set_ylim(0, 100)
    ax1.legend(fontsize=9, loc='upper right')
    ax1.set_title(title, fontsize=12, fontweight='bold')
    ax1.grid(alpha=0.2)
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)

    prdm9_positions = hotspot_data.get('hotspot_positions', [])
    for pos in prdm9_positions:
        ax2.axvline(pos, color='red', alpha=0.7, lw=1)
    ax2.set_ylabel("Motif PRDM9", fontsize=10)
    ax2.set_yticks([])
    ax2.set_xlabel("Posição no Genoma (pb)", fontsize=10)
    ax2.text(
        0.02, 0.85,
        f"n = {hotspot_data.get('n_prdm9', 0)} sítios PRDM9\n"
        f"Densidade: {hotspot_data.get('hotspot_density_per_mb', 0):.1f} / Mb",
        transform=ax2.transAxes, fontsize=9,
        bbox={'boxstyle': 'round', 'facecolor': '#fce4ec', 'alpha': 0.8}
    )
    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)

    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, dpi=150, bbox_inches='tight')
    return fig


def plot_meselson_stahl(
    stats: list[dict],
    title: str = "Experimento de Meselson-Stahl — Replicação Semiconservativa",
    save_path: Optional[str] = None,
) -> plt.Figure:
    generations = [s['generation'] for s in stats]
    heavy  = [s['N15_N15_pct'] for s in stats]
    hybrid = [s['N15_N14_pct'] for s in stats]
    light  = [s['N14_N14_pct'] for s in stats]

    fig, ax = plt.subplots(figsize=(8, 5))
    x = np.arange(len(generations))
    width = 0.5

    ax.bar(x, heavy,  width, label='N¹⁵/N¹⁵ (Pesado)',  color='#c0392b')
    ax.bar(x, hybrid, width, bottom=heavy,               label='N¹⁵/N¹⁴ (Híbrido)', color='#e67e22')
    ax.bar(x, light,  width, bottom=[h + hy for h, hy in zip(heavy, hybrid)],
           label='N¹⁴/N¹⁴ (Leve)',    color='#27ae60')

    ax.set_xticks(x)
    ax.set_xticklabels([f"Geração {g}" for g in generations], fontsize=10)
    ax.set_ylabel("Proporção de Moléculas (%)", fontsize=10)
    ax.set_title(title, fontsize=12, fontweight='bold')
    ax.set_ylim(0, 105)
    ax.legend(fontsize=9, loc='upper right')
    ax.grid(axis='y', alpha=0.3)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    for i, (h, hy, l) in enumerate(zip(heavy, hybrid, light)):
        if h > 3:
            ax.text(i, h/2,          f"{h:.0f}%", ha='center', va='center', fontsize=8, color='white', fontweight='bold')
        if hy > 3:
            ax.text(i, h + hy/2,     f"{hy:.0f}%", ha='center', va='center', fontsize=8, color='white', fontweight='bold')
        if l > 3:
            ax.text(i, h + hy + l/2, f"{l:.0f}%", ha='center', va='center', fontsize=8, color='white', fontweight='bold')

    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, dpi=150, bbox_inches='tight')
    return fig


def plot_msi_analysis(
    msi_result: dict,
    title: str = "Análise de MSI — Instabilidade de Microssatélites",
    save_path: Optional[str] = None,
) -> plt.Figure:
    loci = msi_result.get('loci_details', [])
    if not loci:
        fig, ax = plt.subplots()
        ax.text(0.5, 0.5, "Sem loci comparáveis", ha='center')
        return fig

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
    fig.suptitle(title, fontsize=12, fontweight='bold')

    instable = msi_result['instable_loci']
    stable   = msi_result['total_loci'] - instable
    status   = msi_result['status']
    status_color = {'MSI-H': '#e03030', 'MSI-L': '#ff9800', 'MSS': '#27ae60'}.get(status, 'gray')

    ax1.pie(
        [stable, instable],
        labels=[f'Estável ({stable})', f'Instável ({instable})'],
        colors=['#27ae60', '#e03030'],
        autopct='%1.1f%%',
        startangle=90,
        wedgeprops={'edgecolor': 'white', 'linewidth': 2},
    )
    ax1.set_title(f"Status: {status}", fontsize=12, color=status_color, fontweight='bold')
    ax1.text(0, -1.35,
             f"MSI ratio: {msi_result['msi_ratio']:.2%}\n"
             f"Total loci: {msi_result['total_loci']}",
             ha='center', fontsize=9)

    normal_rep = [l['normal_repeats'] for l in loci]
    tumor_rep  = [l['tumor_repeats']  for l in loci]
    colors     = ['#e03030' if l['instable'] else '#27ae60' for l in loci]

    ax2.scatter(normal_rep, tumor_rep, c=colors, alpha=0.7, s=60, edgecolors='white')
    max_val = max(max(normal_rep, default=1), max(tumor_rep, default=1)) + 2
    ax2.plot([0, max_val], [0, max_val], 'k--', alpha=0.4, lw=1, label='Normal = Tumor')
    ax2.set_xlabel("Repetições (Normal)", fontsize=10)
    ax2.set_ylabel("Repetições (Tumor)", fontsize=10)
    ax2.set_title("Repetições por Locus", fontsize=10)
    ax2.legend(fontsize=8)
    ax2.grid(alpha=0.2)
    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)

    red_patch = mpatches.Patch(color='#e03030', label='Instável')
    green_patch = mpatches.Patch(color='#27ae60', label='Estável')
    ax2.legend(handles=[red_patch, green_patch], fontsize=8)

    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, dpi=150, bbox_inches='tight')
    return fig