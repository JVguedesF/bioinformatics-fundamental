# P2 - Codon Usage Bias & ORF Finder

Ferramenta de genômica comparativa desenvolvida para analisar a eficiência traducional e identificar regiões codificantes em procariotos e eucariotos. O módulo automatiza a aquisição de dados do NCBI e aplica estatística descritiva para mapear o viés de uso de códons (CUB).

---

## 📂 Estrutura do Módulo

O projeto implementa uma arquitetura MVC com uma camada de serviço adicional para comunicação com APIs externas:

### 1. Modelagem de Dados (`src/models.py`)
Define os contratos de dados tipados que sustentam a análise.

* **`ORFResult`**: Dataclass que representa uma *Open Reading Frame* candidata.
    * *Campos:* Frame de leitura (+1 a +3, -1 a -3), coordenadas genômicas (bp), fita (*strand*) e sequência proteica traduzida.
* **`CodonMetrics`**: Estrutura para estatísticas de uso de códons.
    * *Campos:* Contagem total, conteúdo GC, dicionário de frequências absolutas e lista de códons preferenciais.

### 2. Lógica (`src/analyzer.py`)
O núcleo de processamento (`SequenceAnalyzer`), responsável pela manipulação de sequências.

* **`find_orfs()`**: Algoritmo de busca exaustiva em **6 frames de leitura** (3 *sense* + 3 *antisense*).
    * Utiliza tabelas genéticas específicas (ex: Tabela 11 para Bactérias, Tabela 1 Padrão) para tradução correta.
    * Filtra candidatos baseados em um limiar de comprimento mínimo (ex: >100aa).
* **`analyze_codon_usage()`**: Realiza a contagem de tripletos na região codificante.
    * Calcula o viés de frequência e o percentual de GC, essenciais para estudos de adaptação evolutiva.

### 3. Serviço de Dados (`src/downloader.py`)
Módulo de integração (Service Layer) com o **NCBI Entrez**.

* Gerencia a conexão com o banco de dados *Nucleotide*.
* Realiza o download automático e o cache local de genomas completos baseados nos Accession IDs configurados.

### 4. Visualização (`src/view.py`)
Camada de apresentação (`CodonView`) utilizando a biblioteca `Rich`.

* **`create_analysis_view()`**: Renderiza um painel interativo exibindo a métrica genômica e a "Melhor Candidata" a ORF detectada.
* **Relatórios**: Orquestra a exportação dos dados estatísticos para JSON, CSV e gera relatórios em LaTeX prontos para compilação.

### 5. Configuração (`src/config.py`)
Centraliza os parâmetros de execução.

* **`DATASETS`**: Lista de dicionários definindo os organismos alvo (Nome, ID do NCBI e Tipo).
* **`MIN_ORF_LENGTH`**: Define o corte de tamanho para considerar uma ORF válida.

---

## 📊 Inputs e Outputs

* **Input:**
    * IDs do GenBank configurados no `src/config.py` (ex: `NC_000913.3` para *E. coli*).
    * O sistema baixa automaticamente os arquivos `.fasta` para `P2_Codon_Usage_Analysis/data/sequences/`.
* **Output:**
    * **JSON:** Dados completos das ORFs e contagem de códons.
    * **CSV:** Tabela comparativa de uso de códons e conteúdo GC.
    * **LaTeX:** Relatório formatado destacando os códons mais frequentes e estatísticas gerais.