# Genomic Sequence Analysis of Human Coronavirus Species
## Using the Complexity–Entropy Plane

This repository contains the scripts and data used to analyze the complete genomic sequences and Spike proteins of different human coronavirus species, based on the relationship between complexity and entropy.

---

## Structure of the Workflow

The analysis workflow is divided into two main steps:

1. **Full Sequence Data**
2. **Spike Protein Data**

Each step contains specific scripts to generate frequency tables, probability distributions, bar charts, and Complexity–Entropy plane plots.

---

## 1. Full Sequence Data

### 1.1 Get k-mer frequency tables

Use the script ```Frequency_Count_of_k-mers.ipynb``` to generate k-mer frequency tables for each complete genome.

Input: files from the folder ```Complete_Sequences```.

Output: k-mer frequency tables for each genomic sequence.


### 1.2 Obtain normalized probability distribution tables

Use the script ```Normalization_k-mer_Counts.ipynb``` to calculate the normalized probability distributions.

Input: files from the folder ```Frequencies of k-mers```.

Output: tables containing normalized probability distributions for each genomic sequence.



### 1.3 Obtain bar charts of probabilities versus codons

Use the script ```Bar_Chart_Distribution.ipynb``` to generate bar charts showing the occurrence probabilities of each codon.

Input: files from the folder ```Probability Distributions```.

Output: panel containing bar charts.


### 1.4 Obtain Complexity–Entropy planes

Use the script ```Shannon_Complexity-Entropy.ipynb```, ```Tsallis_Complexity-Entropy.ipynb``` and ```Renyi_Complexity-Entropy.ipynb``` to generate Complexity–Entropy planes.

Input: files from the folder ```Probability Distributions```.

Output: Complexity–Entropy plane for Shannon, Tsallis, and Rényi entropies.



## 2. Spike Protein Data

### 2.1 Obtain amino acid frequency tables and probability distributions

Use the script ```Generate AminoAcid Freq and Probability Dist.ipynb``` to generate amino acid count tables for the proteins and normalized probabilities.

Input: files from the folder ```Sequences Spike Protein```. 

Output: frequency tables and probability distributions for each sequence.


### 2.2 Obtain bar charts of probabilities versus amino acids

Use the script ```Distribution_Bar_Chart.ipynb``` para criar o painel com os gráficos de barras das probabilidades versus aminoácidos:

Input: files from the folder ```Probability Distributions```.

Output: panel containing bar charts.


### 2.3 Obtain Complexity–Entropy planes
Use the scripts ```Shannon_Complexity-Entropy.ipynb```, ```Tsallis_Complexity-Entropy.ipynb``` and ```Renyi_Complexity-Entropy.ipynb``` to generate Complexity–Entropy planes of the Spike proteins.

Input: files from the folder ```Probability Distributions```.

Output: Complexity–Entropy plane for Shannon, Tsallis, and Rényi entropies.
