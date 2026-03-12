# Mobile Elements Genomics 🧬

![Status](https://img.shields.io/badge/Status-Active_Development-green)
![Python](https://img.shields.io/badge/Python-3.10+-blue)
![Stack](https://img.shields.io/badge/Stack-Biopython_%7C_Matplotlib_%7C_Docker-2496ED)

Ferramenta de linha de comando para análise comparativa de **arquitetura genômica**, cobrindo três eixos: densidade gênica em genomas organelares e bacterianos, detecção de origens de replicação via **GC Skew**, e rastreamento de **elementos móveis Alu** em DNA nuclear humano.

---

## 🔬 O que o pipeline faz

1. Baixa automaticamente genomas do **NCBI Entrez** (formato GenBank e FASTA).
2. **Análise de Densidade Gênica** — compara genes/kb entre bactérias, organelas, arqueas e vírus como evidência de endossimbiose.
3. **Análise de GC Skew** — calcula skew cumulativo ao longo do genoma e prediz a origem de replicação (oriC).
4. **Rastreamento de Elementos Alu** — varre o cromossomo 19 humano em busca de sequências consenso via regex.
5. Exporta resultados em **JSON** e gráficos **PNG**.

---

## 📂 Estrutura

.
├── src/
│   ├── analyzer.py      # Densidade gênica, GC Skew e scan de elementos móveis
│   ├── config.py        # Datasets e parâmetros de análise
│   ├── downloader.py    # Download via NCBI Entrez
│   ├── main.py          # Orquestrador do pipeline
│   ├── models.py        # Dataclasses e exceções customizadas
│   ├── printer.py       # Interface de terminal (prints padronizados)
│   ├── reporter.py      # Exportação JSON
│   └── view.py          # Geração de gráficos PNG (Matplotlib)
├── data/
│   └── sequences/       # Sequências baixadas (gerado em runtime)
├── results/             # Saídas do pipeline
│   └── graphs/          # Gráficos PNG gerados
├── Dockerfile
├── pyproject.toml
└── .env                 # Credenciais do Entrez (não versionado)

---

## ⚙️ Configuração

Crie um arquivo `.env` na raiz do projeto:
```
ENTREZ_EMAIL=seu_email@exemplo.com
```
---

## 🚀 Execução Local
```
python3 -m venv venv
source venv/bin/activate

pip install -e .

python src/main.py
```
---

## 🐳 Execução via Docker
```
docker build -t mobile-elements-detector .

docker run -it --rm \
  -e ENTREZ_EMAIL=seu@email.com \
  -v $(pwd)/data:/app/data \
  -v $(pwd)/results:/app/results \
  mobile-elements-detector
```
---

## 📊 Resultados

As saídas são organizadas em `results/`:

* **Subpasta `graphs/`**: Contém os arquivos `.png` (Gráfico de densidade, curvas de GC Skew por genoma e o mapa cromossômico de elementos Alu).
* **Arquivos JSON**: Relatórios estruturados com os dados de densidade gênica, coordenadas de origem de replicação e hits do scan Alu.

---

## 🧬 Datasets Principais

O pipeline processa diversos organismos para validar a teoria da endossimbiose e arquitetura genômica:

* **Organelas**: Mitocôndrias humanas, cloroplastos de Arabidopsis thaliana e mitocôndrias de levedura.
* **Bactérias**: Escherichia coli (Modelos), Rickettsia prowazekii (Parasita) e Wolbachia (Simbionte).
* **Arqueas e Vírus**: Halobacterium salinarum, Aeropyrum pernix e Bacteriófago lambda.
* **Nuclear**: Fragmentos do Cromossomo 19 humano para análise de elementos repetitivos.

---

*Desenvolvido como parte de uma trilha de especialização em Bioinformática.*