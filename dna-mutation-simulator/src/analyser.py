from __future__ import annotations
import re
from collections import Counter
from pathlib import Path
from typing import Any, Optional

from Bio import SeqIO

from src.model import (
    Helicase, DNAPolymerase, LaggingPolymerase, Ligase, Primase,
    Mutator, DamageSimulator, DamageEvent,
    MismatchRepair, BaseExcisionRepair, NucleotideExcisionRepair,
    NonHomologousEndJoining, HomologousRecombination,
    Recombinator, Chromatin, MeselsonStahl,
)

COSMIC_PROFILES = {
    'SBS1':  {'type': 'C>T', 'context': 'CpG',  'aetiology': 'Envelhecimento / desaminação de 5mC (clock-like)'},
    'SBS2':  {'type': 'C>T', 'context': 'TCA',  'aetiology': 'Atividade APOBEC (APOBEC3A/3B)'},
    'SBS4':  {'type': 'C>A', 'context': 'CCG',  'aetiology': 'Tabagismo (benzo[a]pireno)'},
    'SBS6':  {'type': 'indel','context': 'rep',  'aetiology': 'Defeito MMR / MSI'},
    'SBS7a': {'type': 'C>T', 'context': 'CC',   'aetiology': 'UV — dímeros de pirimidina (melanoma)'},
    'SBS11': {'type': 'G>A', 'context': 'GCC',  'aetiology': 'Agentes alquilantes (temozolomida)'},
    'SBS18': {'type': 'C>A', 'context': 'ACT',  'aetiology': 'Dano oxidativo (ROS / 8-oxo-G)'},
    'SBS36': {'type': 'C>A', 'context': 'ACG',  'aetiology': 'Defeito BER / MUTYH'},
    'SBS3':  {'type': 'flat','context': 'all',   'aetiology': 'Defeito HR (BRCA1/2) — padrão uniforme'},
}

CANONICAL_SUBS = {'C>A', 'C>G', 'C>T', 'T>A', 'T>C', 'T>G'}
COMPLEMENT = str.maketrans('ACGT', 'TGCA')


def _to_canonical(ref: str, alt: str, context: str) -> Optional[tuple[str, str]]:
    ref, alt, context = ref.upper(), alt.upper(), context.upper()
    sub = f"{ref}>{alt}"
    if sub in CANONICAL_SUBS:
        return sub, context
    ref_c = ref.translate(COMPLEMENT)
    alt_c = alt.translate(COMPLEMENT)
    ctx_c = context.translate(COMPLEMENT)[::-1]
    return f"{ref_c}>{alt_c}", ctx_c


class ReplicationSimulator:

    def __init__(self):
        self.helicase = Helicase()
        self.polymerase = DNAPolymerase()
        self.lagging_polymerase = LaggingPolymerase()
        self.ligase = Ligase()
        self.primase = Primase()
        self.mutator = Mutator()
        self.damage_sim = DamageSimulator(seed=42)
        self.mmr = MismatchRepair()
        self.ber = BaseExcisionRepair()
        self.ner = NucleotideExcisionRepair()
        self.nhej = NonHomologousEndJoining(seed=42)
        self.hr = HomologousRecombination()
        self.recombinator = Recombinator()
        self.chromatin = Chromatin()

    @staticmethod
    def load_genome(file_path: Path) -> str:
        record = SeqIO.read(file_path, "fasta")
        return str(record.seq)

    def _process_segment(self, template: str, primer_size: int, fragment_size: int) -> tuple[str, str]:
        rna_primer = self.primase.synthesize_primer(template, primer_size)
        dna_ext = self.polymerase.synthesize_complementary_strand(template[primer_size:])
        leading = rna_primer + dna_ext
        fragments = self.lagging_polymerase.generate_okazaki_fragments(leading, fragment_size)
        lagging = self.ligase.ligate_fragments(fragments)
        return leading, lagging

    def start_replication(
        self, fasta_path: Path, origin_pattern: str,
        segment_size: int = 5000, fragment_size: int = 1000, primer_size: int = 10,
    ) -> dict:
        genome = self.load_genome(fasta_path)
        origin = self.helicase.find_origin(genome, origin_pattern)
        seg_right = genome[origin: origin + segment_size]
        seg_left = genome[max(0, origin - segment_size): origin]
        leading_r, lagging_r = self._process_segment(seg_right, primer_size, fragment_size)
        leading_l, lagging_l = self._process_segment(seg_left, primer_size, fragment_size)
        return {
            "origin_position": origin,
            "right": {"leading": leading_r, "lagging": lagging_r, "template": seg_right},
            "left":  {"leading": leading_l, "lagging": lagging_l, "template": seg_left},
        }

    def start_replication_multi_origin(
        self, fasta_path: Path, origin_pattern: str,
        segment_size: int = 3000, fragment_size: int = 1000, primer_size: int = 10,
    ) -> list[dict]:
        genome = self.load_genome(fasta_path)
        origins = self.helicase.find_all_origins(genome, origin_pattern)
        if not origins:
            raise ValueError(f"Nenhuma origem encontrada com o padrão '{origin_pattern}'.")
        forks = []
        for origin in origins[:5]:
            seg = genome[origin: origin + segment_size]
            leading, lagging = self._process_segment(seg, primer_size, fragment_size)
            forks.append({"origin": origin, "leading": leading, "lagging": lagging, "template": seg})
        return forks

    @staticmethod
    def calculate_metrics(original: str, synthesized: str) -> dict:
        length = min(len(original), len(synthesized))
        matches = sum(1 for a, b in zip(original[:length], synthesized[:length]) if a == b)
        errors = length - matches
        return {"identity_score": matches / length * 100 if length else 0, "total_errors": errors, "length_compared": length}

    def simulate_all_repair_systems(self, template: str, num_lesions: int = 10, seed: int = 42) -> dict:
        dmg = DamageSimulator(seed=seed)
        results = {}

        strand, _ = dmg.apply_uv_damage(template, 0)
        mutated = strand
        for i in range(0, min(num_lesions * 500, len(template) - 1), 500):
            mutated = self.mutator.generate_point_mutation(mutated, i, 'X')
        m_before = self.calculate_metrics(template, mutated)
        repaired, _ = self.mmr.repair_all(mutated, template,
            [DamageEvent(i, template[i], 'X', 'mismatch') for i in range(0, min(num_lesions * 500, len(template) - 1), 500)])
        results['MMR'] = {'before': m_before, 'after': self.calculate_metrics(template, repaired), 'damage': 'mismatch'}

        damaged_ber, events_ber = dmg.apply_oxidative_damage(template, num_lesions)
        repaired_ber, _ = self.ber.excise_and_repair(damaged_ber, template, events_ber)
        results['BER'] = {'before': self.calculate_metrics(template, damaged_ber), 'after': self.calculate_metrics(template, repaired_ber), 'damage': 'oxidative (G→T)'}

        damaged_ner, events_ner = dmg.apply_uv_damage(template, num_lesions)
        repaired_ner, _ = self.ner.repair(damaged_ner, template, events_ner)
        results['NER'] = {'before': self.calculate_metrics(template, damaged_ner), 'after': self.calculate_metrics(template, repaired_ner), 'damage': 'UV (C→T em CC)'}

        damaged_nhej, events_nhej = dmg.apply_dsb_damage(template, min(num_lesions // 2, 3))
        repaired_nhej, _ = self.nhej.repair(damaged_nhej, events_nhej)
        results['NHEJ'] = {'before': self.calculate_metrics(template, damaged_nhej), 'after': self.calculate_metrics(template, repaired_nhej), 'damage': 'DSB (quebra dupla fita)'}

        damaged_hr, events_hr = dmg.apply_dsb_damage(template, min(num_lesions // 2, 3))
        repaired_hr, _ = self.hr.repair(damaged_hr, template, events_hr)
        results['HR'] = {'before': self.calculate_metrics(template, damaged_hr), 'after': self.calculate_metrics(template, repaired_hr), 'damage': 'DSB (quebra dupla fita)'}

        return results

    @staticmethod
    def compute_mutational_spectrum(sequence: str, mutations: list[tuple[int, str, str]]) -> dict[str, int]:
        spectrum: dict[str, int] = Counter()
        seq = sequence.upper()
        for pos, ref, alt in mutations:
            if pos == 0 or pos >= len(seq) - 1:
                continue
            if seq[pos] != ref.upper():
                continue
            ctx = seq[pos - 1: pos + 2]
            if len(ctx) != 3 or 'N' in ctx:
                continue
            result = _to_canonical(ref, alt, ctx)
            if result:
                sub, canon_ctx = result
                spectrum[f"{canon_ctx[0]}[{sub}]{canon_ctx[2]}"] += 1
        return dict(spectrum)

    @staticmethod
    def classify_cosmic_signature(spectrum: dict[str, int]) -> dict[str, Any]:
        if not spectrum:
            return {"signature": "Desconhecida", "aetiology": "Dados insuficientes"}
        type_counts: Counter = Counter()
        for key, count in spectrum.items():
            m = re.search(r'\[([A-Z]>[A-Z])]', key)
            if m:
                type_counts[m.group(1)] += count
        if not type_counts:
            return {"signature": "Desconhecida", "aetiology": "Formato inválido"}
        dominant_sub = type_counts.most_common(1)[0][0]
        total = sum(type_counts.values())
        ct_ratio = type_counts.get('C>T', 0) / total
        ca_ratio = type_counts.get('C>A', 0) / total
        ga_ratio = type_counts.get('G>A', 0) / total
        hits = []
        if ct_ratio > 0.4:
            cc_ctx = sum(v for k, v in spectrum.items() if '[C>T]' in k and k[0] == 'C')
            cpg_ctx = sum(v for k, v in spectrum.items() if '[C>T]' in k and k[2] == 'G')
            hits.append(('SBS1', COSMIC_PROFILES['SBS1']['aetiology']) if cpg_ctx > cc_ctx else ('SBS7a', COSMIC_PROFILES['SBS7a']['aetiology']))
        if ca_ratio > 0.3:
            hits.append(('SBS18', COSMIC_PROFILES['SBS18']['aetiology']))
        if ga_ratio > 0.3:
            hits.append(('SBS11', COSMIC_PROFILES['SBS11']['aetiology']))
        if not hits:
            hits.append(('SBS3', COSMIC_PROFILES['SBS3']['aetiology']))
        return {"dominant_substitution": dominant_sub, "signatures": hits, "type_counts": dict(type_counts), "top_signature": hits[0][0], "aetiology": hits[0][1]}

    @staticmethod
    def _find_repeats_for_unit(seq: str, unit_len: int, min_repeats: int) -> list[dict]:
        found = []
        len_seq = len(seq)
        for i in range(len_seq - unit_len * min_repeats):
            unit = seq[i: i + unit_len]
            if 'N' in unit:
                continue
            j = i
            while j + unit_len <= len_seq and seq[j: j + unit_len] == unit:
                j += unit_len
            n_repeats = (j - i) // unit_len
            if n_repeats >= min_repeats:
                found.append({'position': i, 'unit': unit, 'n_repeats': n_repeats, 'end': j, 'sequence': seq[i:j]})
        return found

    @staticmethod
    def _filter_overlapping_microsatellites(found: list[dict]) -> list[dict]:
        found.sort(key=lambda x: -x['n_repeats'])
        unique: list[dict] = []
        covered: set[int] = set()
        for ms in found:
            positions = set(range(ms['position'], ms['end']))
            if not positions & covered:
                unique.append(ms)
                covered |= positions
        unique.sort(key=lambda x: x['position'])
        return unique

    @staticmethod
    def find_microsatellites(sequence: str, min_unit: int = 1, max_unit: int = 6, min_repeats: int = 4) -> list[dict]:
        seq = sequence.upper()
        found = []
        for unit_len in range(min_unit, max_unit + 1):
            found.extend(ReplicationSimulator._find_repeats_for_unit(seq, unit_len, min_repeats))
        return ReplicationSimulator._filter_overlapping_microsatellites(found)

    def detect_msi(self, normal_seq: str, tumor_seq: str, min_unit: int = 1, max_unit: int = 4, min_repeats: int = 4) -> dict:
        normal_ms = {ms['position']: ms for ms in self.find_microsatellites(normal_seq, min_unit, max_unit, min_repeats)}
        tumor_ms  = {ms['position']: ms for ms in self.find_microsatellites(tumor_seq,  min_unit, max_unit, min_repeats)}
        common_positions = set(normal_ms.keys()) & set(tumor_ms.keys())
        instable = 0
        loci_details = []
        for pos in common_positions:
            n_rep_normal = normal_ms[pos]['n_repeats']
            n_rep_tumor  = tumor_ms[pos]['n_repeats']
            is_instable = n_rep_normal != n_rep_tumor
            if is_instable:
                instable += 1
            loci_details.append({'position': pos, 'unit': normal_ms[pos]['unit'], 'normal_repeats': n_rep_normal, 'tumor_repeats': n_rep_tumor, 'instable': is_instable})
        total = len(common_positions)
        msi_ratio = instable / total if total > 0 else 0
        status = 'MSI-H' if msi_ratio > 0.4 else ('MSI-L' if msi_ratio > 0.2 else 'MSS')
        return {'total_loci': total, 'instable_loci': instable, 'msi_ratio': msi_ratio, 'status': status, 'loci_details': sorted(loci_details, key=lambda x: x['position'])}

    PRDM9_MOTIF = r'CCNCCNTNNCCNC'

    @staticmethod
    def _prdm9_regex() -> re.Pattern:
        return re.compile(ReplicationSimulator.PRDM9_MOTIF.replace('N', '[ACGT]'), re.IGNORECASE)

    def find_recombination_hotspots(self, sequence: str, window_size: int = 1000, gc_threshold: float = 55.0) -> dict:
        seq = sequence.upper()
        prdm9_sites = [{'position': m.start(), 'sequence': m.group()} for m in self._prdm9_regex().finditer(seq)]
        gc_windows = []
        for i in range(0, len(seq) - window_size, window_size // 2):
            window = seq[i: i + window_size]
            gc = (window.count('G') + window.count('C')) / len(window) * 100
            gc_windows.append({'start': i, 'end': i + window_size, 'gc_pct': gc})
        high_gc = [w for w in gc_windows if w['gc_pct'] >= gc_threshold]
        density = len(prdm9_sites) / max(len(seq) / 1_000_000, 0.001)
        return {'prdm9_sites': prdm9_sites, 'n_prdm9': len(prdm9_sites), 'gc_windows': gc_windows, 'high_gc_windows': high_gc, 'n_high_gc': len(high_gc), 'hotspot_density_per_mb': density, 'hotspot_positions': [s['position'] for s in prdm9_sites]}

    def simulate_recombination(self, fasta_path: Path, crossover_point: int = 50000) -> tuple[str, str]:
        genome = self.load_genome(fasta_path)
        c1 = genome[:100000]
        c2 = c1.lower()
        return self.recombinator.crossing_over(c1, c2, crossover_point)

    def simulate_packing(self, dna: str, nuc_size: int = 146, linker_size: int = 50) -> list:
        return self.chromatin.pack(dna, nuc_size, linker_size)

    @staticmethod
    def simulate_meselson_stahl(template: str, generations: int = 3) -> list[dict]:
        return MeselsonStahl(template).run_generations(generations)

    @staticmethod
    def correlate_repair_defect_with_cancer(repair_system: str) -> dict:
        cancer_map = {
            'MMR': {'genes': ['MLH1', 'MSH2', 'MSH6', 'PMS2'], 'cancers': ['Câncer colorretal', 'Câncer endometrial', 'Síndrome Lynch'], 'signature': 'SBS6 / SBS15 / SBS20 (MSI)', 'biomarker': 'MSI-H / dMMR', 'therapy': 'Pembrolizumab (anti-PD1) — alta resposta em MSI-H'},
            'BER': {'genes': ['OGG1', 'MUTYH', 'NTHL1'], 'cancers': ['Câncer colorretal (MUTYH-Associated Polyposis)', 'Câncer de pulmão'], 'signature': 'SBS18 / SBS36 (C>A oxidativo)', 'biomarker': 'Mutações somáticas G>T', 'therapy': 'Sem terapia alvo específica; PARP inibidores em estudo'},
            'NER': {'genes': ['XPA', 'XPC', 'XPD', 'ERCC1'], 'cancers': ['Xeroderma Pigmentosum → carcinoma de pele', 'Câncer de pulmão (ERCC1)'], 'signature': 'SBS7a / SBS7b (UV — melanoma)', 'biomarker': 'Carga mutacional elevada em pele exposta ao sol', 'therapy': 'Imunoterapia; cisplatina reduzida (ERCC1 baixo)'},
            'NHEJ': {'genes': ['XRCC4', 'LIG4', 'DCLRE1C'], 'cancers': ['Linfomas', 'Leucemias (rearranjos cromossômicos)'], 'signature': 'Rearranjos estruturais / translocações', 'biomarker': 'Instabilidade cromossômica (CIN)', 'therapy': 'Inibidores DNA-PK (M3814 em estudo)'},
            'HR': {'genes': ['BRCA1', 'BRCA2', 'PALB2', 'RAD51'], 'cancers': ['Câncer de mama', 'Câncer de ovário', 'Câncer de próstata', 'Câncer pancreático'], 'signature': 'SBS3 (padrão uniforme — HRD)', 'biomarker': 'HRD score / LOH / TAI / LST', 'therapy': 'Olaparib / Niraparib / Rucaparib (PARP inibidores)'},
        }
        key = repair_system.upper()
        if key not in cancer_map:
            return {"error": f"Sistema '{repair_system}' não reconhecido. Use: {list(cancer_map.keys())}"}
        info = cancer_map[key].copy()
        info['repair_system'] = key
        return info
