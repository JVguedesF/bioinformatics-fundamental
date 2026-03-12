# Bioinformatics Fundamentals 🧬

![Status](https://img.shields.io/badge/Status-Active_Development-green)
![Level](https://img.shields.io/badge/Level-Fundamentals-blue)
![Stack](https://img.shields.io/badge/Stack-Python_%7C_Biopython_%7C_Docker-2496ED)

Este repositório contém a **Fase 1** da minha trilha de especialização em Bioinformática. O objetivo central destes projetos é o aprendizado prático e a consolidação de conceitos fundamentais da Biologia Molecular Computacional. Cada módulo foi construído como um exercício de estudo para aplicar teoria em implementações algorítmicas.

---

## 🔬 Visão Geral dos Projetos

Os projetos abaixo exploram diferentes temas da bioinformática, utilizando dados reais e simulações:

* **Protein & RNA Viewer (`protein-rna-viewer`):** Estudo de biofísica e termodinâmica de ácidos nucleicos, incluindo a predição de estruturas secundárias de RNA e visualização 3D de estruturas proteicas.
* **Codon Bias Analyzer (`codon-bias-analyzer`):** Focado na identificação de ORFs (Open Reading Frames), tradução utilizando tabelas genéticas alternativas e análise estatística de viés de uso de códons (CUB).
* **Splicing Pattern Analyzer (`splicing-pattern-analyzer`):** Análise da arquitetura gênica em eucariotos, realizando o mapeamento de íntrons e éxons e identificando isoformas geradas por splicing alternativo.
* **Mobile Elements Detector (`mobile-elements-detector`):** Investigação de genomas organelares e detecção de elementos móveis (como elementos Alu), além da identificação de origens de replicação via GC Skew.
* **DNA Mutation Simulator (`dna-mutation-simulator`):** Simulação algorítmica de processos de replicação enzimática, modelagem de sistemas de reparo de DNA e análise de assinaturas mutacionais.

---

## 🛠️ Padrões de Engenharia

Mesmo sendo projetos voltados ao aprendizado, todos seguem princípios de desenvolvimento robustos:

1.  **Arquitetura MVC:** Divisão clara entre a Lógica Biológica (`Model`), Interface (`View`) e Orquestração (`Controller`).
2.  **Dockerização:** Cada módulo contém seu próprio ambiente isolado via Docker para garantir reprodutibilidade.
3.  **Qualidade de Código:** Foco em código modular, tipado e legível (Clean Code).

---

## 🚀 Como Executar

Para manter a organização e a especificidade de cada ferramenta, **cada projeto possui suas próprias instruções de instalação e execução.**

Para rodar qualquer um dos projetos (seja localmente ou via Docker), entre no diretório correspondente e consulte o arquivo `README.md` específico daquela pasta.

Como diretriz geral:
* Os comandos devem ser executados a partir da raiz de cada projeto.
* Projetos que acessam o NCBI exigem um arquivo `.env` com a variável `ENTREZ_EMAIL`.

---

> [!NOTE]
> Este repositório é um portfólio de aprendizado pessoal para exploração acadêmica de conceitos biológicos.

*Desenvolvido com foco em Biologia Molecular Computacional e Engenharia de Software.*