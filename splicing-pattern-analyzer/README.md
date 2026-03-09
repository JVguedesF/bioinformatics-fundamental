# Eukaryotic Splicing Dynamics

Ferramenta de análise transcriptômica projetada para mapear a estrutura exôn-intron de genes eucarióticos. O módulo realiza o alinhamento de sequências de mRNA maduro contra genomas de referência para identificar sítios de *splicing*, validar isoformas e calcular métricas de composição nucleotídica.

---

## 📂 Estrutura do Módulo

O projeto estende a arquitetura MVC base com componentes específicos para análise comparativa e relatórios detalhados:

### 1. Modelagem de Dados (`src/models.py`)
Estrutura de dados tipada para representar características genômicas.

* **`Intron`**: Dataclass que encapsula a biologia de um íntron.
    * *Campos:* Coordenadas genômicas (Start/End), sequência nucleotídica, sítios doadores/aceitadores (ex: GT...AG) e status de canonicidade.

### 2. Lógica (`src/analyzer.py`)
O motor de alinhamento (`SplicingAnalyzer`). Diferente dos módulos anteriores, este processa **pares de sequências** (Genoma vs. Transcrito).

* **Algoritmo "Chunk-Based Alignment"**:
    * Fragmenta o mRNA em pequenos blocos (*seeds*) para encontrar correspondências exatas no genoma.
    * Identifica "gaps" (lacunas) entre os blocos alinhados, classificando-os como íntrons.
    * **Validação de Splicing:** Verifica se as bordas do íntron possuem os dinucleotídeos canônicos (GT no doador, AG no aceitador).

### 3. Serviço de Dados (`src/downloader.py`)
Camada de integração com o **NCBI Entrez**.

* Automatiza o download de conjuntos de dados pareados: Gene Genômico (NG_*) e suas Variantes de Transcrito (NM_*).
* Garante que a análise seja feita sempre com sequências de referência atualizadas.

### 4. Visualização (`src/view.py`)
Camada de apresentação (`SplicingView`) focada na topologia do gene.

* **Mapa Genético Visual:** Renderiza uma barra gráfica no terminal representando a estrutura Exon-Intron (proporcional ao tamanho em pb).
* **Tabela de Status:** Exibe cada íntron encontrado com indicadores visuais de qualidade (✔ OK / ⚠ Atenção) para sítios não-canônicos.

### 5. Relatórios (`src/reporter.py`)
Módulo dedicado à persistência de dados complexos (`SplicingReporter`).

* Separa a lógica de salvar arquivos da lógica de exibição.
* Gera artefatos estruturados para cada isoforma analisada (JSON, CSV e tabelas LaTeX customizadas).

### 6. Configuração (`src/config.py`)
Parâmetros de execução e definição do dataset.

* **`DATASETS`**: Mapeia os IDs do GenBank para o gene TP53 humano (Genoma de Referência + Variantes 1 e 2).
* **`AppConfig`**: Define constantes de tolerância para o alinhamento.

---

## 📊 Inputs e Outputs

* **Input:**
    * Pares de arquivos `.fasta` baixados automaticamente (Genoma + mRNA).
    * O sistema detecta automaticamente qual variante de mRNA está sendo processada.
* **Output:**
    * **JSON:** Estrutura hierárquica contendo todos os íntrons mapeados.
    * **CSV:** Planilha detalhada com coordenadas e sequências dos sítios de *splice*.
    * **LaTeX:** Tabela formatada pronta para inserção em artigos científicos, sumarizando a qualidade do *splicing*.