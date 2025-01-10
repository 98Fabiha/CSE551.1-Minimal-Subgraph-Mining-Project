import pandas as pd
import os
from sudaR2 import find_unique_genes, find_minimal_unique_subsets, calculate_fk_msu, calculate_all_subsets


def main():
    # Input and output file paths
    input_file = "C:/Users/HP/PycharmProjects Dataset/sudaaaAlgo/MSU/msu/lib/datasets/cleaned_binary_gene_fpkm_matrix_dataset.csv"
    output_file = "C:/Users/HP/PycharmProjects Dataset/sudaaaAlgo/MSU/msu/lib/fk_msu_results.csv"

    # Check if the input file exists
    if not os.path.exists(input_file):
        raise FileNotFoundError(f"The input file '{input_file}' does not exist. Please verify the path.")

    # Load dataset
    try:
        df = pd.read_csv(input_file)
    except Exception as e:
        raise ValueError(f"Error loading input file: {e}")

    # Rename the first column to "Cancer" for clarity
    if 'Unnamed: 0' in df.columns:
        df.rename(columns={"Unnamed: 0": "Cancer"}, inplace=True)
    elif 'Cancer' not in df.columns:
        raise ValueError("The dataset does not contain the expected first column for 'Cancer'.")

    # Step 1: Identify unique genes
    print("Step 1: Identifying unique genes for each cancer type...")
    unique_genes = find_unique_genes(df)
    print("Unique Genes for each cancer type:\n", unique_genes)  # Debugging: Check unique genes

    # Step 2: Calculate minimal unique subsets
    print("Step 2: Calculating minimal unique subsets...")
    minimal_subsets = find_minimal_unique_subsets(unique_genes)
    print("Minimal Unique Subsets for each cancer type:\n", minimal_subsets)  # Debugging: Check minimal subsets

    # Step 3: Calculate all possible subsets (using a limited or optimized approach)
    print("Step 3: Calculating all possible subsets...")
    try:
        all_subsets = calculate_all_subsets(unique_genes, sample_size=1000)  # Example: Sampling subsets
        print("Sampled Subsets for each cancer type:\n",
              {k: v[:5] for k, v in all_subsets.items()})  # Display first 5 for brevity
    except MemoryError:
        print("MemoryError occurred while calculating all subsets. Skipping this step.")
        all_subsets = {}

    # Step 4: Calculate FK and MSU
    print("Step 4: Calculating FK and MSU values...")
    fk_values, msu_values = calculate_fk_msu(unique_genes, minimal_subsets)
    print("FK Values:\n", fk_values)  # Debugging: Check FK values
    print("MSU Values:\n", msu_values)  # Debugging: Check MSU values

    # Step 5: Compile results into a DataFrame
    print("\nCompiling results into a DataFrame...")
    results = []
    for cancer in unique_genes.keys():
        unique_genes_str = ', '.join(unique_genes[cancer]) if unique_genes[cancer] else 'None'
        all_subsets_str = ', '.join(['{' + ', '.join(subset) + '}' for subset in all_subsets[cancer]]) if all_subsets[cancer] else '{}'
        minimal_subsets_str = ', '.join(['{' + ', '.join(subset) + '}' for subset in minimal_subsets[cancer]]) if minimal_subsets[cancer] else '{}'
        fk = fk_values[cancer]
        msu = msu_values[cancer]

        results.append([cancer, unique_genes_str, all_subsets_str, minimal_subsets_str, fk, msu])

    results_df = pd.DataFrame(results, columns=['Cancer', 'Unique Genes', 'All Subsets', 'Minimal Unique Subsets', 'FK', 'MSU'])

    # Step 6: Save results to a CSV file
    try:
        results_df.to_csv(output_file, index=False)
        print(f"\nResults successfully saved to '{output_file}'.")
    except Exception as e:
        raise IOError(f"Error saving results to file: {e}")

    # Print the results for the first 10 cancer types only
    print("\nResults for the first 10 cancer types:")
    print(results_df.head(10))

if __name__ == "__main__":
    main()
