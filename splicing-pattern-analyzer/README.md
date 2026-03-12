# Eukaryotic Splicing Dynamics 🧬

![Status](https://img.shields.io/badge/Status-Active_Development-green)
![Python](https://img.shields.io/badge/Python-3.10+-blue)
![Stack](https://img.shields.io/badge/Stack-Biopython_%7C_Docker-2496ED)

Ferramenta de análise transcriptómica projetada para mapear a estrutura exão-intrão de genes eucarióticos. O pipeline realiza o alinhamento de sequências de mRNA maduro contra genomas de referência para identificar locais de *splicing*, validar isoformas e calcular métricas de composição nucleotídica.

---

## 🔬 O que o pipeline faz

1. Descarrega automaticamente conjuntos de dados emparelhados (Gene Genómico e Variantes de Transcrito) do **NCBI Entrez**.
2. **Alinhamento "Chunk-Based"** — fragmenta o mRNA em pequenos blocos para encontrar correspondências exatas no genoma.
3. **Mapeamento de Intrões** — identifica lacunas entre os blocos alinhados, classificando-os como intrões.
4. **Validação de Splicing** — verifica se as bordas do intrão possuem os dinucleótidos canónicos (GT no doador, AG no aceitador).
5. Renderiza um mapa genético visual no terminal e exporta relatórios estruturados (JSON e CSV).

---

## 📂 Estrutura

```text
.
├── src/
│   ├── __init__.py
│   ├── analyzer.py
│   ├── config.py
│   ├── downloader.py
│   ├── main.py
│   ├── models.py
│   └── reporter.py
├── data/
│   └── sequences/
├── results/
├── docs/
│   └── final_report.md
├── Dockerfile
├── pyproject.toml
└── .env


```

---

## ⚙️ Configuração

Crie um ficheiro `.env` na raiz do projeto:

```env
ENTREZ_EMAIL=seu_email@exemplo.com


```

---

## 🚀 Execução Local

```bash
python3 -m venv venv
source venv/bin/activate

pip install -e .

python src/main.py


```

---

## 🐳 Execução via Docker

O container isola o ambiente Debian e a instalação de dependências do sistema como o `gcc`.

```bash
docker build -t eukaryotic-splicing-analyzer .

docker run -it --rm \
  -e ENTREZ_EMAIL=seu@email.com \
  -v $(pwd)/data:/app/data \
  -v $(pwd)/results:/app/results \
  eukaryotic-splicing-analyzer


```

---

## 📊 Resultados

As saídas de dados e a documentação técnica são direcionadas para as seguintes pastas:

| Diretório / Arquivo | Formatos | Conteúdo |
| --- | --- | --- |
| `results/` | `.csv` `.json` | Estrutura hierárquica e planilhas detalhadas contendo todos os intrões mapeados e locais de *splice* |
| `docs/final_report.md` | `.md` | Relatório de aprendizado detalhando a dinâmica de splicing eucariótico e o abismo do primeiro intrão do gene TP53 |

---

## 🧬 Datasets

Configuração padrão de alvos para a análise do gene TP53:

| ID NCBI | Ficheiro | Tipo |
| --- | --- | --- |
| `NG_017013.2` | `tp53_human_genomic_refseqgene.fasta` | Genoma de Referência |
| `NM_000546.6` | `tp53_human_transcript_variant1.fasta` | Variante de Transcrito 1 |
| `NM_001126112.2` | `tp53_human_transcript_variant2.fasta` | Variante de Transcrito 2 |
| `NG_007557.1` | `trp53_mouse_genomic_ortholog.fasta` | Ortólogo Genómico |

```
