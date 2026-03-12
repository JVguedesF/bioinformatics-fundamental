Aqui está o `README.md` para o **protein-rna-viewer**, seguindo rigorosamente o mesmo padrão visual, técnico e estrutural que estabelecemos, e respeitando a sua regra de não incluir nenhum comentário dentro dos blocos de código.

Como este projeto envolve ferramentas pesadas de bioinformática estrutural (PyMOL e ViennaRNA), o Docker se torna ainda mais protagonista aqui.

```markdown
# Protein-RNA Viewer 🧬

![Status](https://img.shields.io/badge/Status-Active_Development-green)
![Python](https://img.shields.io/badge/Python-3.10+-blue)
![Stack](https://img.shields.io/badge/Stack-PyMOL_%7C_ViennaRNA_%7C_Biopython_%7C_Docker-2496ED)

Ferramenta de genômica estrutural e biologia computacional dedicada à análise de estruturas tridimensionais. O pipeline foca na visualização de complexos macromoleculares, análise de impacto de mutações em proteínas e predição termodinâmica de dobramento secundário de RNA.

---

## 🔬 O que o pipeline faz

1. Baixa automaticamente estruturas tridimensionais do **PDB** (Protein Data Bank) e sequências do **NCBI**
2. **Análise Estrutural Proteica** — renderiza modelos 3D, mapeia interações e avalia mutações utilizando o motor computacional do **PyMOL**
3. **Predição de Dobramento de RNA** — calcula a estrutura secundária de menor energia livre (MFE) utilizando os algoritmos da suíte **ViennaRNA**
4. **Exportação de Mídia** — gera figuras de alta resolução (PNG e SVG) das topologias 2D e 3D analisadas
5. Consolida as métricas termodinâmicas e estruturais em relatórios formatados

---

## 📂 Estrutura

```text
.
├── src/
│   ├── analyzer.py
│   ├── config.py
│   ├── downloader.py
│   ├── main.py
│   ├── models.py
│   ├── pipeline.py
│   └── reporter.py
├── data/
│   ├── sequences/
│   └── structures/
├── docs/
│   ├── figures/
│   └── final_report.md
├── Dockerfile
├── pyproject.toml
└── .env

```

---

## ⚙️ Configuração

Crie um arquivo `.env` na raiz do projeto:

```env
ENTREZ_EMAIL=seu_email@exemplo.com

```

---

## 🚀 Execução Local

Para rodar nativamente, é obrigatório ter as dependências de sistema do PyMOL e do ViennaRNA instaladas no seu SO (ex: via `apt-get`).

```bash
python3 -m venv venv
source venv/bin/activate

pip install -e ".[rna]"

python src/main.py

```

---

## 🐳 Execução via Docker (Recomendado)

O container já resolve de forma isolada todas as compilações em C/C++ exigidas pelo ViennaRNA e as bibliotecas do PyMOL.

```bash
docker build -t protein-rna-viewer .

docker run -it --rm \
  -e ENTREZ_EMAIL=seu@email.com \
  -v $(pwd)/results:/app/results \
  -v $(pwd)/data:/app/data \
  protein-rna-viewer

```

---

## 📊 Resultados

As saídas visuais e de modelagem são direcionadas para as pastas de dados e documentação:

| Diretório / Arquivo | Formatos | Conteúdo |
| --- | --- | --- |
| `results/` | `.csv` `.json` | Métricas tabulares de biofísica e termodinâmica exportadas pelo pipeline |
| `docs/figures/` | `.png` `.svg` | Renderizações 3D do PyMOL e gráficos de dobramento 2D de RNA |
| `data/structures/` | `.cif` `.pdb` | Arquivos de coordenadas espaciais brutos e manipulados |
| `docs/final_report.md` | `.md` | Relatório com análises de energia livre e impacto mutacional |

---

## 🧬 Datasets

Exemplos de alvos processados por padrão pela ferramenta:

| ID / Arquivo | Tipo | Grupo de Análise |
| --- | --- | --- |
| `1LYZ` | Estrutura PDB | Lisozima (modelo de visualização proteica) |
| `ecoli_tRNA_S88` | Sequência FASTA | Predição de dobramento (folding) de RNA |
| `NM_000207` | Sequência FASTA | Transcrito alvo para análise de mutações |

```

```