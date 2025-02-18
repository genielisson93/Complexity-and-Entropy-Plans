# Essential library imports
import math  # Imports the math module for mathematical operations
import warnings  # Imports warnings to display or suppress warnings
import numpy as np  # Imports NumPy for array manipulation and numerical calculations
import pandas as pd  # Imports Pandas for handling data tables
import matplotlib.pyplot as plt  # Imports Matplotlib for generating graphs


# Function to compute Rényi entropy
def renyi_entropy(data, k, alpha=1):
    """
    Computes the Rényi entropy for a set of probabilities.

    Parameters:
    - data: array-like, list of probabilities.
    - k: int, order of k-mers.
    - alpha: float or array-like, Rényi entropy parameter.

    Returns:
    - s: normalized Rényi entropy value (or array).
    """

    probabilities = np.asarray(data)  # Converts data to a NumPy array.
    probabilities = probabilities[probabilities > 0]  # Removes zero values to avoid log(0) issues.

    smax = np.log2(float(4**k))  # Computes the maximum possible entropy for normalization (4^k for DNA).

    # If alpha is a list/array, compute entropy for each alpha value
    if isinstance(alpha, (tuple, list, np.ndarray)):
        s = []  # List to store entropy values.

        for alpha_ in alpha:
            if alpha_ != 1:  # General case for alpha ≠ 1
                s += [(1 / (1 - alpha_)) * np.log2(np.sum(probabilities**alpha_)) / smax]
            else:  # Special case: Shannon entropy (limit when alpha → 1)
                s += [-np.sum(probabilities * np.log2(probabilities)) / smax]

        s = np.asarray(s)  # Convert the list to a NumPy array.

    else:  # If alpha is a single value, compute directly without looping.
        if alpha != 1:  # General case for alpha ≠ 1
            s = (1 / (1 - alpha)) * np.log2(np.sum(probabilities**alpha)) / smax
        else:  # Special case: Shannon entropy (alpha = 1)
            s = -np.sum(probabilities * np.log2(probabilities)) / smax

    return s  # Returns the computed entropy.



# Function to compute Rényi complexity-entropy
def renyi_complexity_entropy(data, k, alpha=1):
    """
    Computes Rényi entropy and the associated complexity using Jensen-Rényi divergence.

    Parameters:
    - data: array-like, list of probabilities.
    - k: int, order of k-mers.
    - alpha: float or array-like, Rényi entropy parameter.

    Returns:
    - Matrix containing entropy and normalized complexity.
    """

    # Auxiliary function to compute the maximum Jensen-Rényi divergence
    def jensen_renyi_divergence_max(n_states, q):
        """
        Computes the maximum value of Jensen-Rényi divergence for a given number of states.

        Parameters:
        - n_states: number of possible states.
        - q: Rényi entropy parameter.

        Returns:
        - Maximum Jensen-Rényi divergence value.
        """
        if q == 1:  # Special case: limit for Jensen-Shannon divergence (q → 1)
            return (-0.5 * (((n_states + 1) / n_states) * np.log2(n_states + 1) +
                            np.log2(n_states) - 2 * np.log2(2 * n_states)))
        else:  # General formula for q ≠ 1
            return (
                (np.log2(((n_states + 1.)**(1. - q) + n_states - 1.) / (2.**(1. - q) * n_states)) + (1. - q) *
                 np.log2((n_states + 1.) / (2. * n_states))) / (2. * (q - 1))
            )

    n = float(4**k)  # Computes the total number of possible states (4^k for DNA).

    probabilities = np.asarray(data)  # Converts to a NumPy array.
    probabilities = probabilities[probabilities != 0]  # Removes zero values.

    h_a = renyi_entropy(probabilities, k, alpha)  # Computes Rényi entropy.

    n_states_not_occur = n - len(probabilities)  # Number of unobserved states.
    uniform_dist = 1 / n  # Defines the uniform distribution for reference.
    p = probabilities  # Observed probabilities.

    if isinstance(alpha, (tuple, list, np.ndarray)):  # If alpha is a list/array, compute for each value.
        jr_div = []  # List to store Jensen-Rényi divergence values.
        jr_div_max = []  # List to store maximum Jensen-Rényi divergence values.

        for alpha_ in alpha:
            if alpha_ == 1:  # Special case: Shannon entropy
                p_plus_u_over_2 = (uniform_dist + p) / 2
                s_of_p_plus_u_over_2 = -np.sum(p_plus_u_over_2 * np.log2(p_plus_u_over_2)) - (0.5 * uniform_dist) * np.log2(0.5 * uniform_dist) * n_states_not_occur
                s_of_p_over_2 = -np.sum(p * np.log2(p)) / 2
                s_of_u_over_2 = np.log2(n) / 2.
                jr_div += [s_of_p_plus_u_over_2 - s_of_p_over_2 - s_of_u_over_2]
            else:
                first_term = ((p + uniform_dist) / 2)**(1 - alpha_)
                first_term = np.log2(np.sum(first_term * p**alpha_))
                second_term = np.log2(np.sum((1 / (n**alpha_)) * ((p + uniform_dist) / 2)**(1 - alpha_)) + 
                                      (n_states_not_occur * (1 / n**alpha_) * (1 / (2 * n))**(1 - alpha_)))

                jr_div += [(1 / (2 * (alpha_ - 1))) * (first_term + second_term)]

            jr_div_max += [jensen_renyi_divergence_max(n, alpha_)]

        jr_div = np.asarray(jr_div)  # Convert list to NumPy array.
        jr_div_max = np.asarray(jr_div_max)  # Convert list to NumPy array.

    else:  # If alpha is a single value, compute directly without looping.
        if alpha == 1:  # Special case: Shannon entropy
            p_plus_u_over_2 = (uniform_dist + p) / 2
            s_of_p_plus_u_over_2 = -np.sum(p_plus_u_over_2 * np.log2(p_plus_u_over_2)) - (0.5 * uniform_dist) * np.log2(0.5 * uniform_dist) * n_states_not_occur
            s_of_p_over_2 = -np.sum(p * np.log2(p)) / 2
            s_of_u_over_2 = np.log2(n) / 2.
            jr_div = s_of_p_plus_u_over_2 - s_of_p_over_2 - s_of_u_over_2
        else:
            first_term = ((p + uniform_dist) / 2)**(1 - alpha)
            first_term = np.log2(np.sum(first_term * p**alpha))
            second_term = np.log2(np.sum((1 / (n**alpha)) * ((p + uniform_dist) / 2)**(1 - alpha)) + 
                                  (n_states_not_occur * (1 / n**alpha) * (1 / (2 * n))**(1 - alpha)))

            jr_div = (1 / (2 * (alpha - 1))) * (first_term + second_term)

        jr_div_max = jensen_renyi_divergence_max(n, alpha)

    # Return the matrix containing entropy and normalized complexity
    return np.asarray([h_a, h_a * jr_div / jr_div_max]).T

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
    entropies = []  # List to store Rényi entropy values
    complexities = []  # List to store complexity values
    alpha_values = []  # List to store α (alpha) values

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

            # Iterate over alpha values logarithmically ranging from 0.001 to 100
            for alpha in np.logspace(np.log10(0.001), np.log10(100), num=1000):
                try:
                    # Compute Rényi entropy and complexity for each alpha value
                    result = renyi_complexity_entropy(p, alpha=alpha, k=3)

                    # Unpack entropy and complexity values
                    entropy, complexity = result
                    entropies.append(entropy)
                    complexities.append(complexity)
                    alpha_values.append(alpha)
                except ValueError as e:
                    print(f"Error unpacking the result for alpha={alpha}: {e}")
                except Exception as e:
                    print(f"Error computing entropy and complexity for alpha={alpha}: {e}")

        # Plot Entropy vs Complexity points for each Excel file
        plt.scatter(entropies, complexities, color=colors[i % len(colors)], s=8, alpha=0.7,
                    label=file.replace('probabilities_', '').replace('.xlsx', ''))

    except FileNotFoundError:
        print(f"Error: The file '{file}' was not found. Check the file path and name.")
    except KeyError as e:
        print(e)

# Final graph configuration
plt.xticks(fontsize=18)  # Set font size for X-axis values
plt.yticks(fontsize=18)  # Set font size for Y-axis values
plt.xlabel("Rényi Entropy, $H_{\\alpha}$", fontsize=20)  # Label for the X-axis
plt.ylabel("Complexity, $C_{\\alpha}$", fontsize=20)  # Label for the Y-axis
plt.grid(True, linestyle='-', alpha=0.5)  # Add a grid to the graph for better visualization
plt.legend(fontsize=10, loc='lower left')  # Add a legend to differentiate each file in the graph

# Save the generated graph as an image file
plt.savefig('Renyi Complexity-Entropy graph (Nucleotides).png')

# Display the generated graph
plt.show()

