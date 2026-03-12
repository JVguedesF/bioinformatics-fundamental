# Relatório de Aprendizado e Análise Científica - P4: Arquitetura Genômica e Evolução

**Objetivo:** Analisar a estrutura de diferentes genomas para validar a teoria da endossimbiose, localizar origens de replicação e mapear elementos repetitivos.

---

## 1. Densidade Gênica e Endossimbiose
Nesta etapa, comparamos a quantidade de genes por kilobase (kb) entre diversos organismos.

* **Resultados Principais**: A mitocôndria humana apresentou a maior densidade com 2,2331 genes/kb.
* **Comparação**: Organismos bacterianos como a *E. coli* apresentaram densidades próximas a 0,95 genes/kb.
* **Conclusão Biológica**: A alta densidade nas mitocôndrias reflete a perda massiva de DNA não essencial e a compactação genômica ocorrida durante a evolução como endossimbionte.

---

## 2. GC Skew e Origem de Replicação (OriC)
Utilizamos o cálculo de Skew cumulativo para identificar onde a replicação do DNA começa em genomas circulares.

* **E. coli**: O ponto mínimo da curva indicou a origem de replicação na base 3.925.400.
* **Mitocôndria**: A origem foi prevista na base 15.500.
* **Validação**: Os gráficos apresentaram a forma característica de "V", confirmando a transição entre as fitas *leading* e *lagging* no ponto de origem previsto.

---

## 3. Elementos Móveis (Sítios Alu)
Varremos um fragmento do Cromossomo 19 humano em busca de retrotransposons da família Alu.

* **Detecção**: O scan via regex identificou 3 ocorrências do consenso Alu no fragmento analisado.
* **Mapeamento**: O mapa cromossômico gerado demonstrou que os elementos estão concentrados em regiões específicas, evidenciando como essas sequências repetitivas colonizam o DNA nuclear.

---

## Conclusão
O pipeline integrou com sucesso ferramentas de busca por padrões e análises estatísticas. Os dados reforçam que a arquitetura genômica é o resultado de bilhões de anos de otimização energética e transferências horizontais de material genético.