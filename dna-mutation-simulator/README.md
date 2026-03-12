# DNA Mutation Simulator 🧬

![Status](https://img.shields.io/badge/Status-Active_Development-green)
![Python](https://img.shields.io/badge/Python-3.10+-blue)
![Stack](https://img.shields.io/badge/Stack-Biopython_%7C_Matplotlib_%7C_Pandas_%7C_Docker-2496ED)

Ferramenta de biologia computacional concebida para simular a dinâmica genómica e a estabilidade do ADN. O pipeline modela processos de replicação, danos genotóxicos, eficiência de sistemas de reparação e rastreia assinaturas mutacionais.

---

## 🔬 O que o pipeline faz

1. **Dinâmica de Replicação** — simula a replicação bidirecional a partir de origens de replicação (ex: oriC) e emula a experiência histórica de Meselson-Stahl ao longo de várias gerações.
2. **Sistemas de Reparação** — modela lesões no ADN (ex: radiação UV) e avalia algoritmicamente a fidelidade de recuperação utilizando vias de reparação clássicas (MMR, BER, NER, NHEJ, HR).
3. **Assinaturas Mutacionais** — analisa o espetro de mutações geradas e classifica-as com base no catálogo COSMIC, inferindo a etiologia do dano.
4. **Instabilidade Genómica** — deteta e quantifica a Instabilidade de Microssatélites (MSI) cruzando dados de matrizes normais e simuladas (tumorais).
5. **Recombinação e Empacotamento** — mapeia *hotspots* de recombinação (sítios PRDM9), simula eventos de *crossing-over* e projeta o empacotamento em nucleossomas.

---

## 📂 Estrutura

```text
.
├── src/
│   ├── analyser.py
│   ├── config.py
│   ├── downloader.py
│   ├── main.py
│   ├── model.py
│   ├── pipeline.py
│   ├── printer.py
│   ├── reporter.py
│   └── view.py
├── data/
│   └── sequences/
├── results/
│   └── graphs/
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

O container isola todas as dependências do sistema necessárias para as bibliotecas de visualização e manipulação de matrizes (`numpy`, `matplotlib`).

```bash
docker build -t dna-mutation-simulator .

docker run -it --rm \
  -e ENTREZ_EMAIL=seu@email.com \
  -v $(pwd)/data:/app/data \
  -v $(pwd)/results:/app/results \
  dna-mutation-simulator


```

---

## 📊 Resultados

As saídas de dados e as representações visuais são direcionadas para a pasta de resultados:

| Diretório / Arquivo | Formatos | Conteúdo |
| --- | --- | --- |
| `results/graphs/` | `.png` | Visualizações geradas pelo Matplotlib (forquilha de replicação, espectro mutacional, análise de MSI, etc.) |
| `results/` | `.json` `.csv` | Métricas de fidelidade de reparo e análise quantitativa tabular da simulação |

---

## 🧬 Datasets

Para executar os cenários biológicos, o módulo orquestra o descarregamento automático dos seguintes organismos alvo:

| ID NCBI | Organismo | Foco da Análise |
| --- | --- | --- |
| `NC_000913` | *Escherichia coli K12* | Dinâmica da forquilha de replicação |
| `NM_000251` | *Homo sapiens* (MSH2) | Modelo genético para reparação |
| `NM_007294` | *Homo sapiens* (BRCA1) | *Hotspots* de recombinação e reparação HR |
| `NC_001133` | *Saccharomyces cerevisiae* | Simulação de eventos de *crossing-over* |
