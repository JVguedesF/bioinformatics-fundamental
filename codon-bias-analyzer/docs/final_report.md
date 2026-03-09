# Relatório de Aprendizado e Análise Científica - P2: Decodificação Genômica

**Contexto Educacional:** Este relatório documenta a segunda etapa da minha trilha de estudos em Bioinformática. O objetivo deste módulo foi traduzir a teoria do **Dogma Central** para código Python, enfrentando desafios reais como a variação de tabelas genéticas e a análise estatística de viés de códons (CUB).

---

## 1. Escopo e Metodologia

Para compreender na prática como a informação flui do DNA para a Proteína em diferentes domínios da vida, desenvolvi um pipeline de **Tradução em 6 Quadros** (6-Frame Translation) e o apliquei aos seguintes organismos modelo:

| Organismo | Classificação | Desafio de Aprendizado (Tabela Genética) |
| :--- | :--- | :--- |
| **Fago Lambda** | Vírus | Lidar com alta densidade gênica (Table 11) |
| ***E. coli* K-12** | Bactéria | Identificar *Start Codons* alternativos (Table 11) |
| **Homo sapiens** | Eucarioto | Comparar perfil de códons com bactérias (Table 1) |

---

## 2. Resultados da Análise e Interpretação

### 2.1 Otimização Viral (Fago Lambda)
* **Observação:** O script detectou 122 ORFs longas em um genoma minúsculo (48.5 kb).
* **Aprendizado:** Isso validou computacionalmente o conceito de **Genes Sobrepostos**. Aprendi que vírus, limitados pelo tamanho do capsídeo, usam o mesmo DNA lido em frames diferentes para gerar proteínas distintas — uma "compressão de dados" biológica.

### 2.2 Viés de Códons e Wobble (*E. coli*)
* **Observação:** A bactéria prefere fortemente os códons CGC (Arg) e GCG (Ala), ignorando sinônimos.
* **Aprendizado (Wobble Pairing):**
    * Entendi que a redundância do código genético não é aleatória.
    * O viés reflete a **Hipótese do Pareamento Oscilante**: o genoma co-evolui com a disponibilidade de tRNAs da célula. O código confirmou que usar códons "raros" (para os quais há poucos tRNAs) é energeticamente custoso para a bactéria.

### 2.3 Biotecnologia e Insulina Humana
* **Observação:** O gene da insulina humana usa muito o códon **AGC** (Serina), que descobri ser raro na *E. coli*.
* **Conexão Teórica:** Isso explicou um problema real da indústria farmacêutica. Se eu tentasse expressar esse gene "bruto" na bactéria, a produção falharia.
* **Solução Estudada:** O conceito de **Otimização de Códons** (reescrever o gene trocando AGC por um sinônimo bacteriano) é necessário para contornar essa barreira evolutiva.

---

## 3. Competências Consolidadas

Ao finalizar este projeto, consolidei os seguintes conhecimentos:

### 3.1 Bioinformática (Hard Skills)
1.  **Manipulação de Tabelas:** Aprendi a programar a seleção dinâmica de tabelas de tradução (NCBI Table 1 vs 11) para evitar erros de predição em procariotos.
2.  **Algoritmos Genômicos:** Entendi a lógica matemática por trás da busca de ORFs (cálculo de coordenadas relativas vs absolutas).

### 3.2 Biologia Molecular (Conceitos)
1.  **Universalidade vs. Particularidade:** O código genético é universal na estrutura, mas "dialetal" na preferência de códons (CUB).
2.  **Eficiência Traducional:** A escolha do códon impacta diretamente a velocidade do ribossomo.

---