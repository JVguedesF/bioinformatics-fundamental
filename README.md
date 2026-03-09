# Bioinformatics Fundamental 🧬

![Status](https://img.shields.io/badge/Status-Active_Development-green)
![Level](https://img.shields.io/badge/Level-Fundamentals-blue)
![Stack](https://img.shields.io/badge/Stack-Python_%7C_Biopython_%7C_Docker-2496ED)

Este repositório contém a **Fase 1** da minha trilha de especialização em Bioinformática. Aqui estão implementados os conceitos fundamentais da Biologia Molecular Computacional, estruturados com padrões de engenharia de software de mercado (MVC, Docker, Clean Code).

---

## 📂 Visão Geral dos Projetos

| ID                                          | Projeto                  | Foco Biológico                                        | Stack Principal                        |
|:--------------------------------------------|:-------------------------|:------------------------------------------------------|:---------------------------------------|
| **[P1](./P1_Structural_Profiling_Viz)**     | **Structural Profiling** | Biofísica, Termodinâmica (DNA/RNA) e Estrutura 3D     | `Biopython` `ViennaRNA` `PyMOL`        |
| **[P2](./P2_Codon_Usage_Analysis)**         | **Codon Usage Analysis** | Genômica Comparativa, ORFs e Viés de Códons           | `NCBI Entrez` `Rich` `Biopython`       |
| **[P3](./P3_Eukaryotic_Splicing_Dynamics)** | **Eukaryotic Splicing**  | Alinhamento, Mapeamento de Íntrons e Sítios de Splice | `Biopython` `Needleman-Wunsch` (Logic) |
| **[P4](./P4_Mobile_Elements_Genomics)**     | **Mobile Elements**      | Genomas Organelares e Transposons                     | *(Em breve)*                           |
| **[P5](./P5_Replication_Mutation_Sim)**     | **Replication Sim**      | Dinâmica de Replicação e Reparo                       | *(Em breve)*                           |

---

## 🛠️ Padrões de Engenharia

Para garantir reprodutibilidade e organização, todos os projetos seguem uma arquitetura padronizada:

1.  **Monorepo Modular:** Código compartilhado reside na pasta `utils/`, evitando duplicação.
2.  **Arquitetura MVC:** Separação estrita entre Lógica Biológica (`Model`), Interface CLI (`View`) e Orquestração (`Controller`).
3.  **Docker First:** Cada projeto possui seu próprio container isolado para resolver dependências de sistema.

---

## ⚙️ Configuração Obrigatória

Alguns projetos utilizam APIs externas (como o NCBI Entrez). Para identificação e controle de tráfego, você deve configurar suas credenciais.

1. Crie um arquivo `.env` na **raiz absoluta** do repositório.
2. Adicione seu e-mail:
```env
ENTREZ_EMAIL=seu_email@exemplo.com
```

---

## 🚀 Como Executar (Localmente)

Este repositório funciona como um pacote Python unificado. **Todos os comandos devem ser executados a partir da raiz.**

### 1. Instalação

```bash
# 1. Criar e ativar o ambiente virtual
python3 -m venv venv
source venv/bin/activate      # Linux/Mac
# ou: venv\Scripts\activate   # Windows

# 2. Instalar dependências e o módulo base 'utils'
# (O comando lê o arquivo pyproject.toml na raiz)
pip install .
```

### 2. Execução dos Projetos

Basta chamar o script principal de cada módulo:

```bash
# P1: Análise Estrutural (DNA/RNA/Proteína)
python P1_Structural_Profiling_Viz/src/main.py

# P2: Codon Usage e ORFs (Baixa genomas do NCBI)
python P2_Codon_Usage_Analysis/src/main.py

# P3: Splicing Dynamics (Alinhamento Genoma vs mRNA)
python P3_Eukaryotic_Splicing_Dynamics/src/main.py

```

---

## 🐳 Como Executar (Docker)

O uso do Docker é recomendado para garantir que ferramentas externas (como *ViennaRNA*) funcionem sem instalação manual no sistema.

**Importante:** Execute o build sempre a partir da **raiz** (note o `.` no final dos comandos) para que o Docker tenha acesso à pasta `utils`.

### P1 - Structural Profiling

```bash
docker build -t bio-p1 -f P1_Structural_Profiling_Viz/Dockerfile .
docker run --rm -v $(pwd)/P1_Structural_Profiling_Viz/results:/app/P1_Structural_Profiling_Viz/results bio-p1

```

### P2 - Codon Usage

```bash
docker build -t bio-p2 -f P2_Codon_Usage_Analysis/Dockerfile .
docker run --rm -v $(pwd)/P2_Codon_Usage_Analysis/results:/app/P2_Codon_Usage_Analysis/results bio-p2

```

### P3 - Splicing Dynamics

```bash
docker build -t bio-p3 -f P3_Eukaryotic_Splicing_Dynamics/Dockerfile .
docker run --rm -v $(pwd)/P3_Eukaryotic_Splicing_Dynamics/results:/app/P3_Eukaryotic_Splicing_Dynamics/results bio-p3

```

---

*Desenvolvido como parte do meu portfólio pessoal de estudos em Bioinformática.*
