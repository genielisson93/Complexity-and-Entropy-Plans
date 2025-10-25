# Complexity-and-Entropy-Plans


---

## 🎯 Objetivos do Projeto

O repositório busca compreender **padrões informacionais e estruturais** em genomas e proteínas virais através de métricas de **entropia e complexidade**.  
Os principais objetivos são:

1. **Quantificar a organização interna** das sequências genômicas completas e das proteínas Spike.  
2. **Comparar a estrutura estatística** entre diferentes espécies de coronavírus humanos.  
3. **Avaliar padrões de preferência de códons e aminoácidos** (bias de uso) a partir das distribuições de probabilidade.  
4. **Explorar relações não-lineares** entre entropia e complexidade, revelando o grau de ordem informacional em cada sequência.  

---

## 🧩 Seção 1 — Full Sequence Data

Análises baseadas no **genoma completo** das seis espécies de HCoV.

### 📘 Diretórios e Arquivos

#### 🔹 `Frequency of k-mers/`
Planilhas `.xlsx` com as **frequências absolutas** de cada k-mer genômico:

- `kmer_frequency_HCoV-229E.xlsx`  
- `kmer_frequency_HCoV-HKU1.xlsx`  
- `kmer_frequency_HCoV-NL63.xlsx`  
- `kmer_frequency_HCoV-OC43.xlsx`  
- `kmer_frequency_MERS-CoV.xlsx`  
- `kmer_frequency_SARS-CoV-2.xlsx`

#### 🔹 `Probability Distributions/`
Planilhas com as **probabilidades normalizadas** das combinações de bases (k-mers):

- `probabilities_HCoV-*.xlsx`  
- `Uniform Distribution.xlsx` (distribuição teórica uniforme).

#### 🔹 `IMG/`
Gráficos gerados a partir das distribuições de probabilidade e dos planos de complexidade:

- `Codon Probabilities Bar Chart Panel.png`  
- `Shannon Complexity-Entropy graph (Nucleotides).png`  
- `Rényi Complexity-Entropy graph (Nucleotides).png`  
- `Tsallis Complexity-Entropy Graph (Nucleotides).png`

#### 🔹 Notebooks Principais
- **`Frequency count of k-mers.ipynb`** → Lê FASTA e gera tabelas de frequências.  
- **`Normalizes k-mer counts.ipynb`** → Converte contagens em distribuições de probabilidade.  
- **`Bar Chart Distribution.ipynb`** → Cria painéis de barras comparando probabilidades de códons.  
- **`Shannon / Rényi / Tsallis Complexity Entropy.ipynb`** → Calcula entropia e complexidade segundo cada formalismo.

---

## 🧬 Seção 2 — Spike Protein Data

Análises **focadas nas proteínas Spike (S)**, responsáveis pela ligação viral à célula hospedeira.  
Permite comparação direta com os resultados genômicos.

### 📘 Diretórios e Arquivos

#### 🔹 `Amino Acid Frequency/`
Planilhas `.xlsx` com as **frequências absolutas de aminoácidos** para cada espécie:

- `kmer_frequency_Spike HCoV-229E.xlsx`  
- `kmer_frequency_Spike HCoV-HKU1.xlsx`  
- `kmer_frequency_Spike HCoV-NL63.xlsx`  
- `kmer_frequency_Spike HCoV-OC43.xlsx`  
- `kmer_frequency_Spike MERS-CoV.xlsx`  
- `kmer_frequency_Spike SARS-CoV-2.xlsx`

#### 🔹 `Probability Distributions/`
Planilhas de **probabilidades normalizadas** para as proteínas Spike:

- `probabilities_Spike HCoV-*.xlsx`  
- `probabilities_Uniform Distribution.xlsx`

#### 🔹 `IMG/`
Visualizações gráficas derivadas das análises:

- `Barplot Panel AminoAcid Probabilities.png`  
- `Shannon complexity-entropy graph (Amino Acids).png`  
- `Rényi complexity-entropy graph (Amino Acids).png`  
- `Tsallis complexity-entropy graph (Amino Acids).png`

#### 🔹 Notebooks
- **`Shannon Complexity Entropy.ipynb`**  
- **`Rényi Complexity Entropy.ipynb`**  
- **`Tsallis Complexity Entropy.ipynb`**

Cada notebook calcula entropia e complexidade sobre as probabilidades dos 20 aminoácidos canônicos, gerando gráficos comparativos e planos (H × C).

---

## 📊 Resultados Obtidos

| Nível de Análise | Tipo de Dado | Indicadores Principais |
|------------------|--------------|------------------------|
| **Genoma Completo** | k-mers | Preferência por códons terminados em A/U; estrutura estatística mais organizada em SARS-CoV-2 e HKU1 |
| **Proteína Spike** | Aminoácidos | Padrões de uso preferencial e alta complexidade em 229E e NL63; relação entre diversidade de aminoácidos e complexidade informacional |

Os resultados mostram que a **complexidade mede a presença de padrões estruturais**, enquanto a **entropia representa a dispersão probabilística** — juntas, caracterizam o nível de organização informacional das sequências.

---

## ⚙️ Dependências

```bash
pip install biopython pandas numpy matplotlib openpyxl
