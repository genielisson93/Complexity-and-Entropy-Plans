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
Usar o script ```Frequency_Count_of_k-mers.ipynb``` para gerar as tabelas de frequências de *k*-mers para cada genoma completo:

Entrada: arquivos da pasta ```Complete_Sequences```.

### 1.2 Obter tabelas de distribuição de probabilidades normalizadas 

Calcular as distribuições de probabilidade normalizadas com o script ```Normalization_k-mer_Counts.ipynb```:

Entrada: arquivos da pasta ```Frequencies_of_k-mers```.

### 1.3 Obter gráficos de barras das probabilidades versus códons 

Gerar gráficos de barras mostrando as probabilidades de ocorrência de cada códon usando o  script ```Bar_Chart_Distribution.ipynb```:

Entrada: arquivos da pasta ```Probability_Distributions```.

### 1.4 Obter planos de Complexidade–Entropia

Gerar os planos usando diferentes definições de entropia (Shannon, Tsallis e Rényi), usando os scripts:

```Shannon_Complexity-Entropy.ipynb```
```Tsallis_Complexity-Entropy.ipynb```
```Renyi_Complexity-Entropy.ipynb```

Entrada: arquivos da pasta ```Probability_Distributions```.

