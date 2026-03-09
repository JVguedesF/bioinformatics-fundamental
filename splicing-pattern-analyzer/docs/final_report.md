# Relatório de Aprendizado e Análise Científica - P3: Dinâmica de Splicing Eucariótico

**Contexto Educacional:** Este relatório documenta a terceira etapa da trilha de estudos. O objetivo foi transitar da lógica contínua procariótica para a complexidade descontínua dos eucariotos, desenvolvendo um algoritmo capaz de mapear a estrutura de Exons e Íntrons e identificar isoformas.

---

## 1. Escopo e Metodologia

Para processar a informação genética segmentada, desenvolvi um **Alinhador de Splicing** (*Gap-Aware Aligner*). O estudo de caso foi o gene **TP53**, comparando o genoma de referência com variantes de transcritos.

| Alvo de Análise          | Tipo de Dado | Desafio de Aprendizado (Bioinformática)                   |
|:-------------------------|:-------------|:----------------------------------------------------------|
| **Genoma Humano (TP53)** | DNA Genômico | Processar sequências longas não-codificantes (Íntrons)    |
| **Variante 1**           | mRNA Maduro  | Implementar lógica de "salto" (*Chunking*) no alinhamento |
| **Variante 2**           | Isoforma     | Detectar micro-divergências estruturais entre transcritos |

---

## 2. Resultados da Análise e Interpretação

### 2.1 Arquitetura Gênica (O Abismo do Íntron 1)
* **Observação:** O algoritmo identificou que o primeiro íntron do TP53 possui ~10.700 pb, enquanto os exons adjacentes são curtos.
* **Aprendizado:** Validou computacionalmente a desproporção entre DNA não-codificante e codificante. O código teve que simular o **Spliceossomo**, "ancorando" trechos de exons e saltando grandes regiões genômicas para manter a continuidade do alinhamento.

### 2.2 Isoformas e Precisão (*TP53*)
* **Observação:** A comparação automatizada entre Variante 1 e 2 detectou uma diferença de 3 pb (um códon) na fronteira do primeiro íntron (10.754 pb vs 10.757 pb).
* **Aprendizado:**
    * A complexidade biológica reside na combinação de exons (*Splicing Alternativo*).
    * O algoritmo demonstrou precisão ao diferenciar isoformas reais de erros de alinhamento, capturando variações sutis na região 5' UTR ou no início da tradução.

### 2.3 Motivos de Splicing (Canônico vs. Não-Canônico)
* **Observação:** A maioria dos íntrons seguiu a regra **GT...AG**. O sistema de alertas flagrou sítios divergentes (ex: Íntrons 5 e 6 com GC...AG ou CT...GT).
* **Conexão Teórica:** Nem todo processamento de RNA é canônico. O software foi desenhado para registrar exceções biológicas sem interromper a execução, refletindo a flexibilidade necessária em ferramentas de bioinformática.

---

## 3. Competências Consolidadas

### 3.1 Bioinformática (Hard Skills)
1.  **Algoritmos de Alinhamento:** Implementação da lógica *Seed-and-Extend* (Ancoragem e Extensão) para definir fronteiras de íntrons com precisão de base única.
2.  **Engenharia de Software:** Aplicação robusta do padrão **MVC** (Model-View-Controller) e *Separation of Concerns*, isolando a lógica de análise (`analyzer`), visualização (`view`) e persistência de dados em múltiplos formatos (`reporter` gerando JSON, CSV e TeX).
3.  **Visualização de Dados:** Mapeamento visual da estrutura gênica (Exon/Íntron) no terminal utilizando a biblioteca *Rich*.

### 3.2 Biologia Molecular (Conceitos)
1.  **Sítios de Splicing:** Definição computacional de sítios Doadores (5') e Aceitadores (3').
2.  **Estabilidade Termodinâmica:** Análise comparativa de conteúdo GC, observando padrões distintos entre regiões codificantes e não-codificantes.

---