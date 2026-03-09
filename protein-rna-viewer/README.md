# P1 - Structural Profiling & Molecular Visualization

Módulo de bioinformática estrutural focado na caracterização físico-química de sequências biológicas. Implementa um pipeline de detecção automática de macromoléculas (DNA/RNA/Proteína) e aplica algoritmos específicos de análise termodinâmica e estrutural.

---

## 📂 Estrutura do Módulo

O projeto segue o padrão MVC, onde cada arquivo tem uma responsabilidade única na análise biológica:

### 1. Modelagem de Dados (`src/models.py`)
Define a estrutura de dados imutável que trafega pelo pipeline.

* **`AnalysisResult`**: Dataclass que padroniza os resultados.
    * *Campos:* Armazena métricas como `gc_content`, `melting_temp`, `mfe`, `molecular_weight` e `secondary_structure`.
    * *Automação:* Gera automaticamente um timestamp de análise na instanciação.

### 2. Lógica (`src/analyzer.py`)
O núcleo de processamento científico (`CentralDogmaAnalyzer`).

* **`detect_molecule_type()`**: Algoritmo que classifica a sequência baseada na presença de caracteres únicos (U vs T, aminoácidos específicos).
* **Pipeline de DNA**:
    * Calcula Tm (Temperatura de Melting) usando termodinâmica de vizinhos mais próximos (Nearest Neighbor).
* **Pipeline de RNA**:
    * Integração com a biblioteca C **ViennaRNA** para *folding* e cálculo de MFE (Minimum Free Energy).
* **Pipeline de Proteína**:
    * Calcula propriedades físico-químicas (pI, Peso Molecular, Instabilidade) via `Bio.ProtParam`.

### 3. Visualização (`src/view.py`)
Responsável pela apresentação dos dados (`StructuralView`), utilizando a biblioteca `Rich`.

* **`create_visualization()`**: Gera uma árvore de decisão visual, exibindo apenas os ramos relevantes para o tipo de molécula detectada.
* **`create_summary_table()`**: Compila múltiplos resultados em uma tabela comparativa para relatórios em lote.
* **Exportação**: Orquestra a geração de arquivos JSON, CSV e relatórios LaTeX através do módulo utilitário compartilhado.

### 4. Configuração (`src/config.py`)
Centraliza constantes e parâmetros do projeto.

* Define limiares biológicos (ex: `MIN_PROTEIN_LENGTH`).
* Mapeia cores para visualização no terminal.
* Gerencia caminhos de diretórios relativos à raiz do monorepo.

---

## 📊 Inputs e Outputs

* **Input:** Arquivos `.fasta` ou `.fa` localizados em `P1_Structural_Profiling_Viz/data/sequences/`. O sistema detecta o tipo de polímero automaticamente.
* **Output:**
    * **JSON:** Dados brutos de todas as sequências processadas.
    * **CSV:** Tabela sumarizada para análise externa.
    * **LaTeX:** Relatório técnico formatado com estatísticas do lote.