# Import essential libraries
import math  # Imports the math module for mathematical operations
import warnings  # Imports warnings to display or suppress warnings
import numpy as np  # Imports NumPy for array manipulation and numerical calculations
import pandas as pd  # Imports Pandas for handling data tables
import matplotlib.pyplot as plt  # Imports Matplotlib for generating graphs


# Function to calculate the generalized logarithm log_q, used in Tsallis entropy
def logq(x, q=1):
    """
    Computes the generalized logarithm log_q of x for a given parameter q.

    Parameters:
    - x: array-like or scalar, input values.
    - q: float, Tsallis entropy parameter (when q → 1, log_q(x) → log2(x)).

    Returns:
    - log_q(x): value of the generalized logarithm.
    """
    x = np.asarray(x, dtype=float)  # Converts x to a NumPy floating-point array to avoid precision errors.

    if q == 1:  # Special case when q = 1, log_q(x) reduces to the standard log base 2.
        return np.log2(x)
    else:  # Tsallis generalized logarithm formula for q ≠ 1:
        return (x**(1 - q) - 1) / (1 - q)

# Function to compute Tsallis entropy
def tsallis_entropy(data, k, q=1):
    """
    Computes the Tsallis entropy for a set of probabilities.

    Parameters:
    - data: array-like, list of probabilities.
    - k: int, k-mer order.
    - q: float or array-like, Tsallis entropy parameter.

    Returns:
    - s: value (or array) of normalized Tsallis entropy.
    """
    probabilities = np.asarray(data)  # Converts data to a NumPy array.
    probabilities = probabilities[probabilities > 0]  # Removes zero values to avoid log(0).

    if isinstance(q, (tuple, list, np.ndarray)):  # If q is a list/array, compute for each value.
        s = []  # List to store entropy values.

        for q_ in q:
            smax = logq(float(4**k), q_)  # Computes log_q for the normalized maximum entropy.
            lnq_1_over_p = logq(1. / probabilities, q_)  # Applies log_q transformation to inverse probabilities.
            s.append(np.sum(probabilities * lnq_1_over_p) / smax)  # Computes normalized entropy.

        s = np.asarray(s)  # Converts the list to a NumPy array.

    else:  # If q is a single value, compute directly without a loop.
        smax = logq(float(4**k), q)  # Computes log_q for the normalized maximum entropy.
        lnq_1_over_p = logq(1. / probabilities, q)  # Applies log_q transformation to inverse probabilities.
        s = np.sum(probabilities * lnq_1_over_p) / smax  # Computes normalized entropy.

    return s  # Returns the computed entropy.

# Function to compute Tsallis complexity-entropy
def tsallis_complexity_entropy(data, k, q=1):
    """
    Computes Tsallis entropy and the associated complexity using Jensen-Tsallis divergence.

    Parameters:
    - data: array-like, list of probabilities.
    - k: int, k-mer order.
    - q: float or array-like, Tsallis entropy parameter.

    Returns:
    - Matrix with entropy and normalized complexity.
    """

    # Auxiliary function to compute the maximum Jensen-Tsallis divergence
    def jensen_tsallis_divergence_max(n_states, q):
        """
        Computes the maximum Jensen-Tsallis divergence value for a given number of states.

        Parameters:
        - n_states: number of possible states.
        - q: Tsallis entropy parameter.

        Returns:
        - Maximum Jensen-Tsallis divergence value.
        """
        if q == 1:  # Special case: limit for Jensen-Shannon divergence (q → 1)
            return -0.5 * (((n_states + 1) / n_states) * np.log2(n_states + 1) +
                           np.log2(n_states) - 2 * np.log2(2 * n_states))
        else:  # General formula for q ≠ 1
            return ((2**(2 - q)) * n_states - (1 + n_states)**(1 - q) -
                    n_states * (1 + 1 / n_states)**(1 - q) - n_states + 1) / ((1 - q) * (2**(2 - q)) * n_states)

    n = float(4**k)  # Computes the total number of possible states (4^k).

    probabilities = np.asarray(data)  # Converts to a NumPy array.
    probabilities = probabilities[probabilities != 0]  # Removes zero values.

    if isinstance(q, (tuple, list, np.ndarray)):  # If q is a list/array, compute for each value.
        h_q = tsallis_entropy(probabilities, k, q)  # Computes Tsallis entropy.
        jt_div = []  # List to store Jensen-Tsallis divergence values.
        jt_div_max = []  # List to store maximum Jensen-Tsallis divergence values.

        for q_, h_q_val in zip(q, h_q):
            n_states_not_occur = n - len(probabilities)  # Number of non-observed states.
            uniform_dist = 1 / n  # Defines the reference uniform distribution.
            p = probabilities  # Observed probabilities.

            # First part of the Jensen-Tsallis divergence calculation
            first_term = (uniform_dist + p) / (2 * p)
            first_term = -0.5 * np.sum(p * logq(first_term, q_))

            # Second part of the Jensen-Tsallis divergence calculation
            second_term = n * (uniform_dist + p) / 2
            second_term = -(0.5 / n) * (np.sum(logq(second_term, q_)) + logq(0.5, q_) * n_states_not_occur)

            jt_div_val = first_term + second_term  # Jensen-Tsallis divergence
            jt_div_max_val = jensen_tsallis_divergence_max(n, q_)  # Maximum Jensen-Tsallis divergence

            jt_div.append(jt_div_val)  # Stores divergence value.
            jt_div_max.append(jt_div_max_val)  # Stores maximum divergence value.

        jt_div = np.asarray(jt_div)  # Converts list to NumPy array.
        jt_div_max = np.asarray(jt_div_max)  # Converts list to NumPy array.

        # Returns the matrix containing entropy and normalized complexity
        return np.asarray([h_q, h_q * jt_div / jt_div_max]).T

    else:  # If q is a single value, compute directly without looping.
        h_q = tsallis_entropy(probabilities, k, q)  # Compute Tsallis entropy.
        n_states_not_occur = n - len(probabilities)  # Number of unobserved states.
        uniform_dist = 1 / n  # Define the uniform distribution for reference.
        p = probabilities  # Observed probabilities.

        # Calculation of Jensen-Tsallis divergence
        first_term = (uniform_dist + p) / (2 * p)
        first_term = -0.5 * np.sum(p * logq(first_term, q))

        second_term = n * (uniform_dist + p) / 2
        second_term = -(0.5 / n) * (np.sum(logq(second_term, q)) + logq(0.5, q) * n_states_not_occur)

        jt_div = first_term + second_term
        jt_div_max = jensen_tsallis_divergence_max(n, q)

        # Returns the matrix containing entropy and normalized complexity
        return np.asarray([h_q, h_q * jt_div / jt_div_max]).T


# List of Excel files containing probability distributions for each species
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
    entropies = []  # List to store Tsallis entropy values
    complexities = []  # List to store complexity values
    q_values = []  # List to store q values

    try:
        # Read data from the Excel file
        df = pd.read_excel(file)

        # Check if the 'Sequence_ID' column exists in the file
        if "Sequence_ID" not in df.columns:
            raise KeyError(f"The column 'Sequence_ID' was not found in the file {file}")

        # Select columns representing k-mer frequencies (ignoring 'Sequence_ID')
        kmer_columns = df.columns[1:]

        # Iterate over each row in the DataFrame (each sequence/protein)
        for index, row in df.iterrows():
            # Get k-mer frequency values, removing NaN values
            p = row[kmer_columns].dropna().tolist()

            # Check if probabilities sum approximately to 1.0 (normalization)
            if not np.isclose(sum(p), 1.0):
                print(f"Warning: The sum of probabilities in row {index} of file '{file}' is not equal to 1.0. Skipping this row.")
                continue  # Skip the row if probabilities are not normalized

            # Iterate over q values logarithmically ranging from 0.001 to 100
            for q in np.logspace(np.log10(0.001), np.log10(100), num=1000):
                try:
                    # Compute Tsallis entropy and complexity for each q value
                    entropy, complexity = tsallis_complexity_entropy(p, q=q, k=3)

                    # Store the results
                    entropies.append(entropy)
                    complexities.append(complexity)
                    q_values.append(q)
                except ValueError as e:
                    print(f"Error unpacking result for q={q}: {e}")
                except Exception as e:
                    print(f"Error computing entropy and complexity for q={q}: {e}")

        # Plot Entropy vs Complexity points for each Excel file
        plt.scatter(entropies, complexities, color=colors[i % len(colors)], s=20, alpha=0.7,
                    label=file.replace('probabilities_', '').replace('.xlsx', ''))

    except FileNotFoundError:
        print(f"Error: File '{file}' not found. Check the file path and name.")
    except KeyError as e:
        print(e)

# Final graph configuration
plt.xticks(fontsize=18)  # Set font size for X-axis values
plt.yticks(fontsize=18)  # Set font size for Y-axis values
plt.xlabel("Tsallis Entropy, $H_{q}$", fontsize=20)  # Label for the X-axis
plt.ylabel("Complexity, $C_{q}$", fontsize=20)  # Label for the Y-axis
plt.grid(True, linestyle='-', alpha=0.5)  # Add a grid to the graph for better visualization
plt.legend(fontsize=10, loc='lower left')  # Add a legend to differentiate each file in the graph

# Save the generated graph as an image file
plt.savefig('Tsallis Complexity-Entropy Graph (Nucleotides).png')

# Display the generated graph
plt.show()

