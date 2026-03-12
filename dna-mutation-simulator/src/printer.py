class TerminalPrinter:
    @staticmethod
    def header(title, version):
        print("=" * 60)
        print(f"{title} {version}")
        print("=" * 60)

    @staticmethod
    def section(title):
        print("\n" + "=" * 60)
        print(f"▶ {title}")
        print("=" * 60)

    @staticmethod
    def dataset_status(missing):
        if missing:
            print(f"  [AVISO] Ficheiro ausente: {missing}")
        else:
            print("  ✓ Todos os datasets existem localmente.")

    @staticmethod
    def replication_fork(origin, leading_len, lagging_len, primer):
        print(f"  oriC encontrado em: {origin}")
        print(f"  Leading strand:     {leading_len} nt")
        print(f"  Lagging strand:     {lagging_len} nt")
        print(f"  Primer RNA (início): {primer}")

    @staticmethod
    def repair_results(repair_results):
        print("\n  [ Resultados de Reparo ]")
        print(f"  {'SISTEMA':<6} | {'DANO SIMULADO':<25} | {'ANTES':<8} | {'APÓS':<8} | MELHORA")
        print("  " + "-" * 65)
        for sys_name, data in repair_results.items():
            before = data['before']['identity_score']
            after  = data['after']['identity_score']
            delta  = after - before
            print(f"  {sys_name:<6} | {data['damage']:<25} | {before:>6.2f}% | {after:>6.2f}% | {delta:>+6.2f}%")

    @staticmethod
    def meselson_stahl(ms_stats):
        for gen_data in ms_stats:
            print(f"  Gen {gen_data['generation']}  |  N15/N15: {gen_data['N15_N15_pct']:.0f}%  |  N15/N14: {gen_data['N15_N14_pct']:.0f}%  |  N14/N14: {gen_data['N14_N14_pct']:.0f}%")

    @staticmethod
    def cosmic_signature(classified):
        print("\n  [ Resultado COSMIC ]")
        print(f"    - Substituição dominante: {classified.get('dominant_substitution', 'N/A')}")
        print(f"    - Assinatura provável:    {classified.get('top_signature', 'N/A')}")
        print(f"    - Etiologia:              {classified.get('aetiology', 'N/A')}")

    @staticmethod
    def msi_validation(ms_normal_count, msi_control, msi_result):
        print(f"  Microssatélites encontrados no MSH2: {ms_normal_count}")
        print(f"  [Controle (normal vs normal)] Loci: {msi_control['total_loci']} | Instáveis: {msi_control['instable_loci']} | Status: {msi_control['status']}")
        print(f"  [Tumor (normal vs tumor)] Loci: {msi_result['total_loci']} | Instáveis: {msi_result['instable_loci']} | Status: {msi_result['status']}")
        if msi_control['status'] == 'MSS' and 'MSI' in msi_result['status']:
            print("    ✓ Validação aprovada: controle MSS, tumor MSI detetado corretamente.")
        else:
            print("    ⚠ Resultado inesperado — rever parâmetros de deteção.")

    @staticmethod
    def recombination_hotspots(hotspots):
        print(f"  Sítios PRDM9 encontrados: {hotspots['n_prdm9']}  |  Janelas GC ≥ 55%: {hotspots['n_high_gc']}  |  Densidade: {hotspots['hotspot_density_per_mb']:.1f} / Mb")

    @staticmethod
    def crossing_over(c1, c2):
        print("  Ponto de crossing-over: 50000")
        print(f"  Cromossomo 1 [49995:50005]: {c1}")
        print(f"  Cromossomo 2 [49995:50005]: {c2}")

    @staticmethod
    def chromatin_packing(packed, nucleosomes, linkers):
        print(f"  Total de blocos: {len(packed)}  |  Nucleossomos: {len(nucleosomes)}  |  Linkers: {len(linkers)}")
        for block in packed[:3]:
            print(f"  [{block[0]}] {len(block[1])} nt → {block[1][:12]}...")

    @staticmethod
    def repair_defects(sys_name, info):
        print(f"\n  [ {sys_name} — {info.get('biomarker', '')} ]")
        print(f"    - Genes:      {', '.join(info.get('genes', []))}")
        print(f"    - Cânceres:   {', '.join(info.get('cancers', []))}")
        print(f"    - Assinatura: {info.get('signature', '')}")
        print(f"    - Terapia:    {info.get('therapy', '')}")

    @staticmethod
    def footer(results_dir):
        print(f"\n✓ Pipeline completo! Resultados em: {results_dir}")
        print("\n" + "=" * 60)
        print("Projeto 5 — Concluído")
        print("=" * 60)