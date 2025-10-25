# Análise de Sequências Genômicas de Espécies de Coronavírus Humanos
## Utilizando o Plano Complexidade–Entropia

Este repositório contém os scripts e dados utilizados para analisar as sequências genômicas completas e proteínas Spike de diferentes espécies de **coronavírus humanos**, com base na relação entre **complexidade** e **entropia**.

---

## Estrutura do Roteiro

O fluxo de análise está dividido em duas etapas principais:

1. **Full Sequence Data**
2. **Spike Protein Data**

Cada etapa possui scripts específicos para gerar tabelas de frequências, distribuições de probabilidades, gráficos de barras e planos de Complexidade–Entropia.

---

## 🔹 1. Full Sequence Data

### 1.1 Obter tabelas de frequências de *k*-mers
Usar o script abaixo para gerar as tabelas de frequências de *k*-mers para cada genoma completo:

```bash
Frequency_Count_of_k-mers.ipynb

Arquivos de entrada:

SARS-CoV-2

MERS-CoV

HCoV-229E

HCoV-HKU1

HCoV-NL63

HCoV-OC43

