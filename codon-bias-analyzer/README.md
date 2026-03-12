# Codon Bias Analyzer 🧬

![Status](https://img.shields.io/badge/Status-Active_Development-green)
![Python](https://img.shields.io/badge/Python-3.10+-blue)
![Stack](https://img.shields.io/badge/Stack-Biopython_%7C_Docker-2496ED)

Ferramenta de genómica comparativa desenvolvida para analisar a eficiência traducional e identificar regiões codificantes (ORFs) em procariotas e eucariotas. O pipeline automatiza a aquisição de dados do NCBI e aplica estatística descritiva para mapear o viés de uso de codões (CUB).

---

## 🔬 O que o pipeline faz

1. Descarrega automaticamente genomas completos do **NCBI Entrez** (base de dados Nucleotide).
2. **Busca de ORFs** — executa algoritmos de busca exaustiva em 6 grelhas de leitura (3 *sense* e 3 *antisense*) utilizando tabelas genéticas específicas.
3. **Análise de Viés de Codões** — realiza a contagem de tripletos na região codificante, calcula o viés de frequência e a percentagem de conteúdo GC.
4. **Visualização no Terminal** — apresenta métricas estruturadas de forma limpa diretamente na saída padrão (stdout).
5. Exporta relatórios estatísticos em **JSON** e **CSV**.

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

O container isola as dependências do sistema operativo, como os compiladores `gcc` e `make` necessários.

```bash
docker build -t codon-bias-analyzer .

docker run -it --rm \
  -e ENTREZ_EMAIL=seu@email.com \
  -v $(pwd)/data:/app/data \
  -v $(pwd)/results:/app/results \
  codon-bias-analyzer


```

---

## 📊 Resultados

Os ficheiros são exportados para o diretório `results/`:

* **`.json`**: Dados completos das ORFs e contagem absoluta de codões.
* **`.csv`**: Sumário tabular de uso de codões e conteúdo GC.

---

## 🧬 Datasets

Os identificadores do GenBank a serem processados são configurados no ficheiro `src/config.py`. O sistema suporta a análise de múltiplos organismos em lote.