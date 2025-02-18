# Import essential libraries
import math  # Imports the math module for mathematical operations
import warnings  # Imports warnings to display or suppress warnings
import numpy as np  # Imports NumPy for array manipulation and numerical calculations
import pandas as pd  # Imports Pandas for handling data tables
import matplotlib.pyplot as plt  # Imports Matplotlib for generating graphs


# Function to calculate Shannon entropy (or permutation entropy)
def shannon_entropy(data, k, normalized=True):
    """
    Computes Shannon entropy for a set of probabilities.

    Parameters:
    - data: array-like, list of probabilities.
    - k: int, order of k-mers (size of analyzed words).
    - normalized: bool, if True, normalizes entropy to the range [0,1].

    Returns:
    - Shannon entropy (normalized or not).
    """

    probabilities = np.asarray(data, dtype=np.float64)  # Converts input data to a NumPy array.
    probabilities = probabilities[probabilities > 0]  # Removes null probabilities to avoid log(0).

    # Compute Shannon entropy using the formula:
    # H = - Σ p_i * log2(p_i), where p_i are the probabilities of the distribution.
    s = -np.sum(probabilities * np.log2(probabilities))

    if normalized:  # Checks if entropy should be normalized to a range from 0 to 1.
        smax = k * np.log2(float(4))  # Computes the maximum possible entropy for normalization.
        return s / smax  # Returns normalized entropy.
    else:
        return s  # Returns entropy without normalization.


# Function to compute complexity-entropy based on Jensen-Shannon divergence
def complexity_entropy(data, k):
    """
    Computes Shannon entropy and the associated complexity using Jensen-Shannon divergence.

    Parameters:
    - data: array-like, list of probabilities.
    - k: int, order of k-mers (size of analyzed words).

    Returns:
    - h: normalized Shannon entropy.
    - complexity: complexity based on Jensen-Shannon divergence.
    """

    n = float(4 ** k)  # Defines the total number of possible states (4^k for DNA).

    probabilities = np.asarray(data, dtype=np.float64)  # Converts data to a NumPy double-precision array.
    probabilities = probabilities[probabilities > 0]  # Removes zero values to avoid logarithm errors.

    h = shannon_entropy(probabilities, k)  # Computes normalized Shannon entropy.

    n_states_not_occuring = n - len(probabilities)  # Counts how many states did not occur in the distribution.
    uniform_dist = 1 / n  # Defines the uniform distribution for comparison.

    # Computes the average distribution between the original and the uniform distribution.
    p_plus_u_over_2 = (uniform_dist + probabilities) / 2

    # Computes the entropy of the mixed distribution (real and uniform).
    s_of_p_plus_u_over_2 = (
        -np.sum(p_plus_u_over_2 * np.log2(p_plus_u_over_2))  # Main term
        - (0.5 * uniform_dist) * np.log2(0.5 * uniform_dist) * n_states_not_occuring  # Correction for unobserved states
    )

    # Computes entropy of the original distribution divided by 2.
    s_of_p_over_2 = -np.sum(probabilities * np.log2(probabilities)) / 2

    # Computes entropy of the uniform distribution divided by 2.
    s_of_u_over_2 = np.log2(n) / 2.0

    # Computes the theoretical maximum Jensen-Shannon divergence.
    js_div_max = -0.5 * (
        ((n + 1) / n) * np.log2(n + 1) + np.log2(n) - 2 * np.log2(2 * n)
    )

    # Computes Jensen-Shannon divergence of the original distribution relative to the uniform distribution.
    js_div = s_of_p_plus_u_over_2 - s_of_p_over_2 - s_of_u_over_2

    # Returns Shannon entropy and normalized complexity.
    return h, h * js_div / js_div_max


dist_prob_species = [
    "probabilities_SARS-CoV-2.xlsx",
    "probabilities_HCoV-229E.xlsx",
    "probabilities_HCoV-HKU1.xlsx",
    "probabilities_HCoV-NL63.xlsx",
    "probabilities_HCoV-OC43.xlsx",
    "probabilities_MERS-CoV.xlsx",
    "Uniform Distribution.xlsx"  # Uniform distribution for reference
]

# Define different colors for each file in the graph
colors = ['blue', 'green', 'black', 'red', 'orange', 'purple', 'lightblue', 'pink']

# Configure the figure size for the graph
plt.figure(figsize=(10, 7))

# Loop through each Excel file to process data and generate graphs
for i, file in enumerate(dist_prob_species):
    entropies = []  # List to store Shannon entropy values
    complexities = []  # List to store complexity values

    try:
        # Read data from the Excel file
        df = pd.read_excel(file)

        # Check if the 'Sequence_ID' column exists in the file
        if "Sequence_ID" not in df.columns:
            raise KeyError(f"The column 'Sequence_ID' was not found in the file {file}")

        # Select the columns representing k-mer frequencies (ignoring 'Sequence_ID')
        kmer_columns = df.columns[1:]

        # Iterate over each row in the DataFrame (each sequence/protein)
        for index, row in df.iterrows():
            # Get k-mer frequency values
            p = row[kmer_columns].tolist()

            # Check if probabilities sum approximately to 1.0 (normalization)
            if not np.isclose(sum(p), 1.0):
                print(f"Warning: The sum of probabilities in row {index} of file '{file}' is not equal to 1.0. Skipping this row.")
                continue  # Skip the row if probabilities are not normalized

            try:
                # Compute entropy and complexity using the complexity_entropy function
                entropy, complexity = complexity_entropy(p, k=3)

                # Store results
                entropies.append(entropy)
                complexities.append(complexity)
            except ValueError as e:
                print(f"Error unpacking the result: {e}")
            except Exception as e:
                print(f"Error computing entropy and complexity: {e}")

        # Plot Entropy vs Complexity points for each Excel file
        plt.scatter(entropies, complexities, color=colors[i % len(colors)], s=80, alpha=0.7, edgecolor='white',
                    label=file.replace('probabilities_', '').replace('.xlsx', ''))

    except FileNotFoundError:
        print(f"Error: The file '{file}' was not found. Check the file path and name.")
    except KeyError as e:
        print(e)

# Final graph configuration
plt.xticks(fontsize=18)  # Set font size for X-axis values
plt.yticks(fontsize=18)  # Set font size for Y-axis values
plt.xlabel("Shannon Entropy, H", fontsize=20)  # Label for the X-axis
plt.ylabel("Complexity, C", fontsize=20)  # Label for the Y-axis
plt.grid(True, linestyle='--', alpha=0.5)  # Add a dashed grid to the graph for better visualization
plt.legend(fontsize=10, loc='lower left')  # Add a legend to differentiate each file in the graph

# Save the generated graph as an image file
plt.savefig('Shannon Complexity-Entropy graph (Nucleotides).png')

# Display the generated graph
plt.show()

