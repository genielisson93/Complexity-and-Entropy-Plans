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
Usar o script ```Frequency_Count_of_k-mers.ipynb``` para gerar as tabelas de frequências de *k*-mers para cada genoma completo:

Entrada: arquivos da pasta ```Complete_Sequences```
Saída: tabelas de frequências de k-mers para cada sequência genômica.


### 1.2 Obter tabelas de distribuição de probabilidades normalizadas 

Calcular as distribuições de probabilidade normalizadas com o script ```Normalization_k-mer_Counts.ipynb```:

Entrada: arquivos da pasta ```Frequencies_of_k-mers```
Saída: tabelas contendo distribuições de probabilidades normalizadas para cada sequência genômica.


### 1.3 Obter gráficos de barras das probabilidades versus códons 

Gerar gráficos de barras mostrando as probabilidades de ocorrência de cada códon usando o  script ```Bar_Chart_Distribution.ipynb```:

Entrada: arquivos da pasta ```Probability_Distributions```
Saída: painel contendo gráficos de barras.


### 1.4 Obter planos de Complexidade–Entropia

Gerar os planos, para as definições de entropia de Shannon, Tsallis e Rényi, usando os scripts:

```Shannon_Complexity-Entropy.ipynb```
```Tsallis_Complexity-Entropy.ipynb```
```Renyi_Complexity-Entropy.ipynb```

Entrada: arquivos da pasta ```Probability_Distributions```
Saída: plano complexidade-entropia para as entropias de Shannon, Tsallis e Rényi.


## 2. Spike Protein Data

### 2.1 Obter tabelas de frequências de aminoácidos e distribuições de probabilidade

Gerar tabelas de contagem e probabilidades normalizadas dos aminoácidos das proteínas Spike usando o script ```Generate_AminoAcid_Freq_and_Probability_Dist.ipynb```:

Entrad: arquivos da pasta ```Sequences_Spike_Protein```
Saída: tabelas de frequências e distribuições de probabilidades para cada sequência.

### 2.2 Obter gráficos de barras das probabilidades versus aminoácidos

Usar o script ```Distribution_Bar_Chart.ipynb``` para criar o painel com os gráficos de barras das probabilidades versus aminoácidos:

Entrada: arquivos da pasta ```Probability_Distributions```
Saída: painel contendo gráficos de barras.

### 2.3 Obter planos de Complexidade–Entropia

Calcular e visualizar os planos de Complexidade–Entropia das proteínas Spike, usando os scripts:

```Shannon_Complexity-Entropy.ipynb```
```Tsallis_Complexity-Entropy.ipynb```
```Renyi_Complexity-Entropy.ipynb```

Entrada: arquivos da pasta ```Probability_Distributions```
Saída: plano complexidade-entropia para as entropias de Shannon, Tsallis e Rényi.
