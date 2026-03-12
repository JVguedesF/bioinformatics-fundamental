class TerminalPrinter:
    @staticmethod
    def header():
        print("=" * 60)
        print("ENDOSSIMBIOSE & ELEMENTOS MÓVEIS")
        print("=" * 60)

    @staticmethod
    def section(title):
        print("\n" + "=" * 60)
        print(f"▶ {title}")
        print("=" * 60)

    @staticmethod
    def density_table(results):
        print(f"\n  {'GENOMA':<35} | {'GRUPO':<20} | {'GENES/KB'}")
        print("  " + "-" * 65)
        for res in results:
            print(f"  {res.genome_name:<35} | {res.group:<20} | {res.density_genes_per_kb:>8.2f}")

    @staticmethod
    def skew_result(name, ori_pos):
        print(f"  - {name}: Origem prevista na base {ori_pos}")

    @staticmethod
    def mobile_summary(name, count):
        print(f"\n  [ Resultado do Scan ]")
        print(f"    - Alvo: {name}")
        print(f"    - Elemento: Alu Consensus")
        print(f"    - Ocorrências Encontradas: {count}")