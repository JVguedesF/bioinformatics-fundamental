from __future__ import annotations
import random
import re
from dataclasses import dataclass
from typing import Optional

class BioPipelineError(Exception):
    pass

class SequenceAnalysisError(BioPipelineError):
    pass

class EmptySequenceError(SequenceAnalysisError):
    pass

class Helicase:
    @staticmethod
    def find_origin(genome: str, origin_pattern: str) -> int:
        pos = genome.find(origin_pattern)
        if pos == -1:
            raise ValueError(f"Origem '{origin_pattern}' não encontrada.")
        return pos

    @staticmethod
    def find_all_origins(genome: str, origin_pattern: str) -> list[int]:
        return [m.start() for m in re.finditer(origin_pattern, genome)]


class DNAPolymerase:
    def __init__(self):
        self.complement_table = str.maketrans("ATCGatcg", "TAGCtagc")

    def synthesize_complementary_strand(self, template_strand: str) -> str:
        return template_strand.translate(self.complement_table)


class LaggingPolymerase:
    @staticmethod
    def generate_okazaki_fragments(strand: str, size: int) -> list[str]:
        return [strand[i: i + size] for i in range(0, len(strand), size)]


class Ligase:
    @staticmethod
    def ligate_fragments(fragments: list[str]) -> str:
        return "".join(fragments)


class Primase:
    def __init__(self):
        self.rna_table = str.maketrans("ATCGatcg", "UAGCuagc")

    def synthesize_primer(self, template_strand: str, size: int) -> str:
        return template_strand[:size].translate(self.rna_table)


@dataclass
class DamageEvent:
    position: int
    original_base: str
    damaged_base: str
    damage_type: str   # 'UV' | 'oxidative' | 'alkylation' | 'DSB'
    repaired: bool = False


class Mutator:
    @staticmethod
    def generate_point_mutation(sequence: str, position: int, new_base: str) -> str:
        if position >= len(sequence):
            return sequence
        return sequence[:position] + new_base + sequence[position + 1:]


class DamageSimulator:
    def __init__(self, seed: Optional[int] = None):
        self._rng = random.Random(seed)

    def apply_uv_damage(self, sequence: str, num_lesions: int = 5) -> tuple[str, list[DamageEvent]]:
        seq = list(sequence)
        events: list[DamageEvent] = []
        candidates = [i for i in range(len(seq) - 1) if seq[i].upper() == 'C' and seq[i + 1].upper() == 'C']
        chosen = self._rng.sample(candidates, min(num_lesions, len(candidates)))
        for pos in sorted(chosen):
            orig = seq[pos]
            seq[pos] = 'T'
            events.append(DamageEvent(pos, orig, 'T', 'UV'))
        return "".join(seq), events

    def apply_oxidative_damage(self, sequence: str, num_lesions: int = 5) -> tuple[str, list[DamageEvent]]:
        seq = list(sequence)
        events: list[DamageEvent] = []
        candidates = [i for i, b in enumerate(seq) if b.upper() == 'G']
        chosen = self._rng.sample(candidates, min(num_lesions, len(candidates)))
        for pos in sorted(chosen):
            orig = seq[pos]
            seq[pos] = 'T'
            events.append(DamageEvent(pos, orig, 'T', 'oxidative'))
        return "".join(seq), events

    def apply_alkylation_damage(self, sequence: str, num_lesions: int = 5) -> tuple[str, list[DamageEvent]]:
        seq = list(sequence)
        events: list[DamageEvent] = []
        candidates = [i for i, b in enumerate(seq) if b.upper() == 'G']
        chosen = self._rng.sample(candidates, min(num_lesions, len(candidates)))
        for pos in sorted(chosen):
            orig = seq[pos]
            seq[pos] = 'A'
            events.append(DamageEvent(pos, orig, 'A', 'alkylation'))
        return "".join(seq), events

    def apply_dsb_damage(self, sequence: str, num_breaks: int = 2) -> tuple[str, list[DamageEvent]]:
        seq = list(sequence)
        events: list[DamageEvent] = []
        candidates = list(range(1, len(seq) - 1))
        chosen = self._rng.sample(candidates, min(num_breaks, len(candidates)))
        for pos in sorted(chosen):
            orig = seq[pos]
            seq[pos] = 'X'
            events.append(DamageEvent(pos, orig, 'X', 'DSB'))
        return "".join(seq), events


class MismatchRepair:
    @staticmethod
    def repair_segment(error_strand: str, template: str, start: int, end: int) -> str:
        end = min(end, len(template), len(error_strand))
        return error_strand[:start] + template[start:end] + error_strand[end:]

    @staticmethod
    def repair_all(damaged: str, template: str, events: list[DamageEvent]) -> tuple[str, int]:
        seq = damaged
        count = 0
        for ev in events:
            p = ev.position
            if p < len(seq) and p < len(template):
                seq = seq[:p] + template[p] + seq[p + 1:]
                ev.repaired = True
                count += 1
        return seq, count


class BaseExcisionRepair:
    @staticmethod
    def excise_and_repair(sequence: str, template: str, events: list[DamageEvent]) -> tuple[str, list[DamageEvent]]:
        seq = list(sequence)
        remaining: list[DamageEvent] = []
        for ev in events:
            if ev.damage_type in ('oxidative', 'alkylation') and ev.position < len(seq):
                if ev.position < len(template):
                    seq[ev.position] = template[ev.position]
                    ev.repaired = True
                    continue
            remaining.append(ev)
        return "".join(seq), remaining


def _repair_window(seq: list[str], template: str, position: int, window_size: int) -> None:
    half = window_size // 2
    s = max(0, position - half)
    e = min(len(seq), position + half)
    for i in range(s, e):
        if i < len(template):
            seq[i] = template[i]


class NucleotideExcisionRepair:
    def __init__(self, window_size: int = 27):
        self.window_size = window_size

    def repair(self, sequence: str, template: str, events: list[DamageEvent]) -> tuple[str, list[DamageEvent]]:
        seq = list(sequence)
        remaining: list[DamageEvent] = []
        for ev in events:
            if ev.damage_type == 'UV':
                _repair_window(seq, template, ev.position, self.window_size)
                ev.repaired = True
            else:
                remaining.append(ev)
        return "".join(seq), remaining


class NonHomologousEndJoining:
    def __init__(self, error_rate: float = 0.3, seed: Optional[int] = None):
        self.error_rate = error_rate
        self._rng = random.Random(seed)

    def repair(self, sequence: str, events: list[DamageEvent]) -> tuple[str, list[DamageEvent]]:
        seq = sequence
        dsb = sorted([ev for ev in events if ev.damage_type == 'DSB'], key=lambda e: e.position, reverse=True)
        other = [ev for ev in events if ev.damage_type != 'DSB']
        for ev in dsb:
            pos = ev.position
            if pos < len(seq) and seq[pos] == 'X':
                if self._rng.random() < self.error_rate:
                    del_size = self._rng.randint(1, 3)
                    seq = seq[:pos] + seq[pos + del_size:]
                else:
                    seq = seq[:pos] + seq[pos + 1:]
                ev.repaired = True
        return seq, other


class HomologousRecombination:
    def __init__(self, window_size: int = 50):
        self.window_size = window_size

    def repair(self, damaged: str, homologous: str, events: list[DamageEvent]) -> tuple[str, list[DamageEvent]]:
        seq = list(damaged)
        remaining: list[DamageEvent] = []
        for ev in events:
            if ev.damage_type == 'DSB' and ev.position < len(seq):
                _repair_window(seq, homologous, ev.position, self.window_size)
                ev.repaired = True
            else:
                remaining.append(ev)
        return "".join(seq), remaining


class Recombinator:
    @staticmethod
    def crossing_over(c1: str, c2: str, point: int) -> tuple[str, str]:
        return c1[:point] + c2[point:], c2[:point] + c1[point:]

    @staticmethod
    def double_crossing_over(c1: str, c2: str, p1: int, p2: int) -> tuple[str, str]:
        new_c1 = c1[:p1] + c2[p1:p2] + c1[p2:]
        new_c2 = c2[:p1] + c1[p1:p2] + c2[p2:]
        return new_c1, new_c2


class Chromatin:
    @staticmethod
    def pack(dna: str, nuc_size: int = 146, linker_size: int = 50) -> list[tuple[str, str]]:
        result: list[tuple[str, str]] = []
        pos, total = 0, len(dna)
        while pos < total:
            nuc = dna[pos: pos + nuc_size]
            if nuc:
                result.append(('nucleosome', nuc))
            pos += nuc_size
            link = dna[pos: pos + linker_size]
            if link:
                result.append(('linker', link))
            pos += linker_size
        return result

    @staticmethod
    def compaction_ratio(dna: str, nuc_size: int = 146, linker_size: int = 50) -> float:
        unit = nuc_size + linker_size
        n = len(dna) / unit
        return len(dna) / max(n, 1)


@dataclass
class StrandLabel:
    sequence: str
    label: str  # 'N15' ou 'N14'


class MeselsonStahl:
    def __init__(self, template: str):
        self.population: list[tuple[StrandLabel, StrandLabel]] = [
            (StrandLabel(template, 'N15'), StrandLabel(template, 'N15'))
        ]

    def replicate_generation(self) -> list[tuple[StrandLabel, StrandLabel]]:
        new_pop: list[tuple[StrandLabel, StrandLabel]] = []
        for a, b in self.population:
            new_pop.append((a, StrandLabel(a.sequence, 'N14')))
            new_pop.append((b, StrandLabel(b.sequence, 'N14')))
        self.population = new_pop
        return self.population

    def run_generations(self, n: int) -> list[dict]:
        results = [self._stats(0)]
        for gen in range(1, n + 1):
            self.replicate_generation()
            results.append(self._stats(gen))
        return results

    def _stats(self, generation: int) -> dict:
        heavy, hybrid, light = 0, 0, 0
        for a, b in self.population:
            labels = {a.label, b.label}
            if labels == {'N15'}:
                heavy += 1
            elif labels == {'N14'}:
                light += 1
            else:
                hybrid += 1
        total = len(self.population)
        return {
            'generation': generation,
            'total': total,
            'N15_N15_pct': heavy / total * 100,
            'N15_N14_pct': hybrid / total * 100,
            'N14_N14_pct': light / total * 100,
        }
