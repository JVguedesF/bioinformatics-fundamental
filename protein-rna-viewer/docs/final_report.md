# P1: Perfilagem Estrutural e Biofísica Molecular

**Resumo Executivo:** Este relatório documenta a análise *in silico* de três macromoléculas biológicas (DNA, RNA e Proteína), aplicando os princípios do Dogma Central da Biologia Molecular. O objetivo foi validar computacionalmente como as propriedades primárias da sequência (composição de bases e aminoácidos) determinam a estabilidade termodinâmica e a conformação tridimensional final das moléculas.

---

## 1. Análise Genômica (DNA)

### 1.1 Objetivo e Dados
Avaliar o perfil de estabilidade térmica do transcrito NM_000207 (gene humano).

### 1.2 Resultados Computacionais

| Arquivo (FASTA) | Tipo | Tamanho (bp) | Conteúdo GC (%) | Tm Estimado |
| :--- | :--- | :--- | :--- | :--- |
| `NM_000207.fasta` | DNA | 465 | 63.9% | 1524.0 °C |

### 1.3 Discussão e Insights Biológicos
* **Estabilidade Genômica:** O alto conteúdo GC (63.9%) indica a presença de três pontes de hidrogênio por par de base, conferindo alta estabilidade à fita dupla original e maior resistência à desnaturação térmica.
* **⚠️ Insight Crítico sobre o Algoritmo (Limitação do Tm):** O algoritmo estimou uma Temperatura de Melting (Tm) de 1524°C, um valor biologicamente impossível. Isso evidencia uma limitação computacional importante: o script utilizou a **Regra de Wallace** ($4 \times [G+C] + 2 \times [A+T]$), que é válida apenas para oligonucleotídeos curtos (primers de até 20 pb). Para sequências longas como esta, futuros pipelines devem implementar o modelo termodinâmico de **Nearest-Neighbor (Vizinho Mais Próximo)**.

---

## 2. Termodinâmica de RNA e Estrutura Secundária

### 2.1 Objetivo
Prever a estrutura secundária do tRNA de *E. coli* e avaliar o impacto termodinâmico de uma mutação pontual no anticódon.

### 2.2 Resultados (Wild-Type vs Mutante U35A)

**Sequência Selvagem (WT):**
`GGGUCGUUAGCUCAGUUGGUAGAGCAGUUGACUUUUAAUCAAUUGGUCGCAGGUUCGAAUCCUGCACGACCCA`
* **MFE (Energia Livre Mínima):** -28.60 kcal/mol

**Mutante U35A:**
`GGGUCGUUAGCUCAGUUGGUAGAGCAGUUGACUU[A]UAAUCAAUUGGUCGCAGGUUCGAAUCCUGCACGACCCA`
* **MFE (Energia Livre Mínima):** -28.60 kcal/mol ($\Delta G = 0.00$)



### 2.3 Discussão e Insights Biológicos
A simulação demonstrou que a mutação U35A não causou instabilidade termodinâmica ($\Delta G = 0$). Isso ocorre porque o nucleotídeo 35 está localizado no **loop do anticódon**, uma região de fita simples.
* **Conceito Chave:** Mutações em regiões de *loop* geralmente mantêm a estrutura global intacta, afetando apenas o reconhecimento de substratos. Se a mutação ocorresse em uma região de *stem* (hastes pareadas), a quebra das pontes de hidrogênio causaria o colapso estrutural da molécula.

---

## 3. Biofísica de Proteínas e Modelagem 3D

### 3.1 Objetivo
Caracterizar as propriedades físico-químicas da Lisozima (PDB ID: 1LYZ) e simular a perda de pontes dissulfeto.

### 3.2 Propriedades Físico-Químicas

| Proteína | Tamanho | Peso Molecular | Ponto Isoelétrico (pI) | Índice de Instabilidade | GRAVY |
| :--- | :--- | :--- | :--- | :--- | :--- |
| Lisozima | 129 aa | 14.3 kDa | **9.32** | 16.09 (Estável) | -0.47 |



### 3.3 Discussão e Insights Biológicos
* **Adaptação Funcional pelo pI:** O pI elevado (9.32) indica que a proteína é carregada positivamente em pH fisiológico (7.4). Isso não é uma coincidência: a lisozima ataca bactérias cujas paredes celulares são carregadas negativamente. A atração eletrostática é o primeiro passo de sua função catalítica.
* **Índice GRAVY:** O valor negativo (-0.47) confirma o caráter hidrofílico global da proteína, essencial para sua solubilidade em secreções como lágrimas e saliva.
* **Simulação de Mutação (CYS $\to$ ALA):** A substituição teórica de cisteínas por alaninas revelou a importância crítica das ligações covalentes (pontes dissulfeto). A perda do enxofre desestabilizaria o núcleo hidrofóbico, levando à desnaturação da proteína e perda total de função.

---

## 4. Conclusão Geral

Este projeto validou com sucesso a relação Sequência-Estrutura-Função. Através da bioinformática, foi possível observar que:
1. Em ácidos nucleicos (DNA/RNA), a estabilidade é puramente regida pelo pareamento de bases.
2. Nas proteínas, propriedades emergentes como a carga superficial (pI) moldam a função biológica.

**Próximos Passos:** Integrar modelos de *Nearest-Neighbor* para cálculos termodinâmicos de DNA longo.

---