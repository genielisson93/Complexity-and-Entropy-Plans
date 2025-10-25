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

## 1. Full Sequence Data

### 1.1 Obter tabelas de frequências de *k*-mers
Usar o script ```Frequency_Count_of_k-mers.ipynb``` para gerar as tabelas de frequências de *k*-mers para cada genoma completo.

Entrada: arquivos da pasta ```Complete_Sequences```

Saída: tabelas de frequências de k-mers para cada sequência genômica.


### 1.2 Obter tabelas de distribuição de probabilidades normalizadas 

Usar o script ```Normalization_k-mer_Counts.ipynb``` para Calcular as distribuições de probabilidade normalizadas.

Entrada: arquivos da pasta ```Frequencies of k-mers```.

Saída: tabelas contendo distribuições de probabilidades normalizadas para cada sequência genômica.



### 1.3 Obter gráficos de barras das probabilidades versus códons 

Usar o script ```Bar_Chart_Distribution.ipynb``` para gerar gráficos de barras mostrando as probabilidades de ocorrência de cada códon.

Entrada: arquivos da pasta ```Probability Distributions```.

Saída: painel contendo gráficos de barras.



### 1.4 Obter planos de Complexidade–Entropia

Usar os scripts ```Shannon_Complexity-Entropy.ipynb```, ```Tsallis_Complexity-Entropy.ipynb``` e ```Renyi_Complexity-Entropy.ipynb``` para gerar planos Complexidade-Entropia.


Entrada: arquivos da pasta ```Probability Distributions```.

Saída: plano Complexidade-Entropia para as entropias de Shannon, Tsallis e Rényi.



## 2. Spike Protein Data

### 2.1 Obter tabelas de frequências de aminoácidos e distribuições de probabilidade

Usar o script ```Generate AminoAcid Freq and Probability Dist.ipynb``` para gerar tabelas de contagem de aminoácidos das proteínas e probabilidades normalizadas.

Entrada: arquivos da pasta ```Sequences Spike Protein```. 

Saída: tabelas de frequências e distribuições de probabilidades para cada sequência.


### 2.2 Obter gráficos de barras das probabilidades versus aminoácidos

Usar o script ```Distribution_Bar_Chart.ipynb``` para criar o painel com os gráficos de barras das probabilidades versus aminoácidos:

Entrada: arquivos da pasta ```Probability Distributions```.

Saída: painel contendo gráficos de barras.



### 2.3 Obter planos de Complexidade–Entropia

Usar os scripts ```Shannon_Complexity-Entropy.ipynb```, ```Tsallis_Complexity-Entropy.ipynb``` e ```Renyi_Complexity-Entropy.ipynb``` para gerar planos de Complexidade–Entropia das proteínas Spike.

Entrada: arquivos da pasta ```Probability Distributions```.

Saída: plano Complexidade-Entropia para as entropias de Shannon, Tsallis e Rényi.
