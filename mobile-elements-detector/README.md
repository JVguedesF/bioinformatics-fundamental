# Mobile Elements Genomics 🧬

![Status](https://img.shields.io/badge/Status-Active_Development-green)
![Python](https://img.shields.io/badge/Python-3.11+-blue)
![Stack](https://img.shields.io/badge/Stack-Biopython_%7C_Rich_%7C_Matplotlib_%7C_Docker-2496ED)

Ferramenta de linha de comando para análise comparativa de **arquitetura genômica**, cobrindo três eixos: densidade gênica em genomas organelares e bacterianos, detecção de origens de replicação via **GC Skew**, e rastreamento de **elementos móveis Alu** em DNA nuclear humano.

---

## 🔬 O que o pipeline faz

1. Baixa automaticamente genomas do **NCBI Entrez** (formato GenBank e FASTA)
2. **Análise de Densidade Gênica** — compara genes/kb entre bactérias, organelas, arqueas e vírus como evidência de endossimbiose
3. **Análise de GC Skew** — calcula skew cumulativo ao longo do genoma e prediz a origem de replicação (oriC)
4. **Rastreamento de Elementos Alu** — varre o cromossomo 19 humano em busca de sequências consenso via regex
5. Exporta resultados em **CSV**, **JSON** e gráficos **PNG**

---

## 📂 Estrutura

```
.
├── src/
│   ├── analyzer.py      # Densidade gênica, GC Skew e scan de elementos móveis
│   ├── config.py        # Datasets e parâmetros de análise
│   ├── downloader.py    # Download via NCBI Entrez
│   ├── exceptions.py    # Exceções customizadas
│   ├── models.py        # Dataclasses de resultados
│   ├── reporter.py      # Exportação CSV / JSON
│   ├── view.py          # Geração de gráficos PNG (Matplotlib)
│   └── main.py          # Orquestrador do pipeline
├── data/
│   └── sequences/       # Sequências baixadas (gerado em runtime)
├── results/             # Saídas do pipeline (gerado em runtime)
├── Dockerfile
├── pyproject.toml
└── .env                 # Suas credenciais (não versionado)
```

---

## ⚙️ Configuração

Crie um arquivo `.env` na raiz do projeto:

```env
ENTREZ_EMAIL=seu_email@exemplo.com
```

> O NCBI exige um e-mail válido para uso da API Entrez. O pipeline encerra com erro caso não esteja configurado.

---

## 🚀 Execução Local

```bash
# 1. Criar e ativar ambiente virtual
python3 -m venv venv
source venv/bin/activate       # Linux/Mac
# venv\Scripts\activate        # Windows

# 2. Instalar o projeto
pip install -e .

# 3. Rodar
mobile-elements-detector
```

Na primeira execução, o pipeline baixa os genomas automaticamente. Nas execuções seguintes, reutiliza os arquivos em `data/sequences/`.

---

## 🐳 Execução via Docker

```bash
# Build
docker build -t mobile-elements-detector .

# Run
docker run \
  -e ENTREZ_EMAIL=seu@email.com \
  -v $(pwd)/data:/app/data \
  -v $(pwd)/results:/app/results \
  -it mobile-elements-detector
```

> Os volumes `-v` são necessários para persistir os dados baixados e os resultados após o container encerrar.

---

## 📊 Resultados

Todos os arquivos são exportados em `results/` com timestamp:

| Prefixo | Formatos | Conteúdo |
|---|---|---|
| `density_*` | `.csv` `.json` | Densidade gênica por genoma |
| `gc_skew_*` | `.csv` `.json` | Resumo e arrays completos de skew |
| `mobile_elements_*` | `.csv` `.json` | Hits do scan Alu com posição e sequência |
| `density_plot_*` | `.png` | Gráfico de barras comparativo |
| `gc_skew_<genome>_*` | `.png` | Curva de skew cumulativo com oriC previsto |
| `chromosome_map_*` | `.png` | Mapa cromossômico dos elementos Alu |

---

## 🧬 Datasets

| ID NCBI | Organismo | Grupo |
|---|---|---|
| `NC_012920` | *Homo sapiens* — Mitocôndria | organelle\_animal |
| `NC_001807` | *Homo sapiens* — Mitocôndria (referência antiga) | organelle\_animal |
| `U00096` | *Escherichia coli* K-12 | bacteria\_model |
| `NC_000913` | *Escherichia coli* MG1655 | bacteria\_model |
| `NC_000963` | *Rickettsia prowazekii* | bacteria\_parasite |
| `NC_002947` | *Wolbachia pipientis* | bacteria\_symbiont |
| `NC_000932` | *Arabidopsis thaliana* — Cloroplasto | organelle\_plant |
| `NC_001224` | *Saccharomyces cerevisiae* — Mitocôndria | organelle\_fungi |
| `NC_002607` | *Halobacterium salinarum* | archaea |
| `NC_000854` | *Aeropyrum pernix* | archaea |
| `NC_001416` | Bacteriófago lambda | virus |
| `NT_113818` | Cromossomo 19 humano (parcial, rico em Alu) | nuclear\_target |

---

*Desenvolvido como parte de uma trilha de especialização em Bioinformática.*