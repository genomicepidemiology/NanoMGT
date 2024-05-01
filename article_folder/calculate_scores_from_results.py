import pandas as pd
import sys
import os

def load_mutations_from_files(file_paths):
    """
    Loads mutations from multiple files and returns a dictionary with gene IDs as keys
    and sets of mutations as values.

    :param file_paths: List of file paths to load mutations from
    :return: Dictionary with gene IDs as keys and sets of mutations as values
    """
    mutations_dict = {}

    for file_path in file_paths:
        with open(file_path, 'r') as file:
            lines = file.readlines()
            for i in range(0, len(lines), 2):
                gene_id = lines[i].strip()
                mutations = set(lines[i + 1].strip().split(','))

                if gene_id in mutations_dict:
                    mutations_dict[gene_id].update(mutations)
                else:
                    mutations_dict[gene_id] = mutations

    return mutations_dict

# Adjusted function to load the data from the CSV file, accounting for initial spaces in column names
def load_data(filepath):
    # Added skipinitialspace=True to handle any initial spaces in column names
    return pd.read_csv(filepath, skipinitialspace=True)


# Function to find the ID with the highest percentage for a given batch ID
def find_highest_percentage_id(batch_id, df):
    # Filter the DataFrame for the given batch ID
    batch_data = df[df['Batch'] == batch_id]
    if batch_data.empty:
        return "Batch ID not found."
    # Iterate over the rows and find the ID with the highest percentage
    highest_percentage = 0
    highest_percentage_id = ''
    all = []
    for index, row in batch_data.iterrows():
        if 'Percentage3' in df.columns:
            for i in range(1, 4):
                # Extracting percentage value correctly after handling potential space issue
                percentage = int(row[f'Percentage{i}'][:-1])  # Remove the '%' and convert to int
                if percentage > highest_percentage:
                    highest_percentage = percentage
                    highest_percentage_id = row[f'ID{i}']
                all.append((row[f'ID{i}']))
        else:
            for i in range(1, 3):
                # Extracting percentage value correctly after handling potential space issue
                percentage = int(row[f'Percentage{i}'][:-1])  # Remove the '%' and convert to int
                if percentage > highest_percentage:
                    highest_percentage = percentage
                    highest_percentage_id = row[f'ID{i}']
                all.append((row[f'ID{i}']))

    minor = []

    for item in all:
        if item != highest_percentage_id:
            minor.append(item)

    return highest_percentage_id, minor


def load_mutations(filename):
    mutations_dict = {}
    with open(filename, 'r') as file:
        next(file)  # Skip the header line
        for line in file:
            # Split the line into components
            allele, position, _, mutation_base, _, _, _, _, _ = line.strip().split(',')
            gene = allele.split('_')[0]
            # Generate the mutation identifier as 'Position_MutationBase'
            mutation_identifier = f"{position}_{mutation_base}"

            # If the gene is not already in the dictionary, add it with an empty set
            if gene not in mutations_dict:
                mutations_dict[gene] = set()

            # Add the mutation identifier to the gene's set
            mutations_dict[gene].add(mutation_identifier)

    return mutations_dict


def benchmark_analysis_result(sample, results_folder):
    batch = sample.split('_')[-2]
    batch_id = int(sample.split('_')[-2][5:])
    #print(f"Batch ID: {batch_id}")
    seed = sample[4:7]
    specie = sample.split('_')[1:-2]
    batch_csv = seed + '_' + "_".join(specie) + "_" + 'batches' + ".csv"
    data = load_data('/home/people/malhal/data/new_nanomgt/simulated_batches/' + batch_csv)

    # Change this batch_id to test different batches
    # batch_id = 10
    #print(f"Highest percentage ID for batch {batch_id}: {find_highest_percentage_id(batch_id, data)}")

    top_id, minor = find_highest_percentage_id(batch_id, data)
    # Use names and batch ID to get the correct mutation map
    if 'salmonella_enterica' in sample:
        map_file_1 = '/home/people/malhal/data/new_nanomgt/majority_variants/minor_variant_maps/major_{}_minor_{}.txt'.format(
            top_id, minor[0])
        map_file_2 = '/home/people/malhal/data/new_nanomgt/majority_variants/minor_variant_maps/major_{}_minor_{}.txt'.format(
            top_id, minor[1])
        mutation_map = load_mutations_from_files([map_file_1, map_file_2])
    if 'ecoli' in sample:
        map_file_1 = '/home/people/malhal/data/new_nanomgt/majority_variants/minor_variant_maps/major_{}_minor_{}.txt'.format(
            top_id, minor[0])
        mutation_map = load_mutations_from_files([map_file_1])
    if 'staph_aureus' in sample:
        map_file_1 = '/home/people/malhal/data/new_nanomgt/majority_variants/minor_variant_maps/major_{}_minor_{}.txt'.format(
            top_id, minor[0])
        map_file_2 = '/home/people/malhal/data/new_nanomgt/majority_variants/minor_variant_maps/major_{}_minor_{}.txt'.format(
            top_id, minor[1])
        mutation_map = load_mutations_from_files([map_file_1, map_file_2])
    if 'campylobacter_jejuni' in sample:
        map_file_1 = '/home/people/malhal/data/new_nanomgt/majority_variants/minor_variant_maps/major_{}_minor_{}.txt'.format(
            top_id, minor[0])
        mutation_map = load_mutations_from_files([map_file_1])
    print ('My expected output')
    for item in mutation_map:
        print(item, mutation_map[item])

    return mutation_map


def calculate_metrics(expected_mutations, actual_mutations):
    # Initialize variables to calculate sum of metrics across all genes
    sum_precision, sum_recall, sum_f1 = 0, 0, 0
    genes_counted = 0

    # Iterate over expected mutations by gene to calculate metrics
    for gene, expected_set in expected_mutations.items():
        if gene in actual_mutations:
            actual_set = actual_mutations[gene]
            # True Positives (TP): Mutations correctly predicted
            tp = len(expected_set & actual_set)
            # False Positives (FP): Mutations incorrectly predicted (not in expected but in actual)
            fp = len(actual_set - expected_set)
            # False Negatives (FN): Mutations missed (in expected but not in actual)
            fn = len(expected_set - actual_set)

            # Calculate Precision, Recall, and F1 for this gene
            precision = tp / (tp + fp) if (tp + fp) > 0 else 0
            recall = tp / (tp + fn) if (tp + fn) > 0 else 0
            f1 = 2 * (precision * recall) / (precision + recall) if (precision + recall) > 0 else 0

            # Add to sum of metrics
            sum_precision += precision
            sum_recall += recall
            sum_f1 += f1
            genes_counted += 1

    # Calculate average metrics across all genes
    avg_precision = sum_precision / genes_counted if genes_counted > 0 else 0
    avg_recall = sum_recall / genes_counted if genes_counted > 0 else 0
    avg_f1 = sum_f1 / genes_counted if genes_counted > 0 else 0

    return avg_precision, avg_recall, avg_f1


# Main function to demonstrate functionality
if __name__ == "__main__":
    # Load the data
    samples = os.listdir('/home/people/malhal/test/run_nanomgt_on_all_for_alignments')
    for sample in samples:
        if sample.endswith('merged') and sample == 'seed101_salmonella_enterica_batch10_merged':
            results_folder = '/home/people/malhal/test/run_nanomgt_on_all_for_alignments/' + sample
            minor_mutation_expected = benchmark_analysis_result(sample, results_folder)

            minor_mutation_results = load_mutations(results_folder + '/minor_mutations.csv')
            print ('my actual output')
            for item in minor_mutation_results:
                print(item, minor_mutation_results[item])

            precision, recall, f1 = calculate_metrics(minor_mutation_expected, minor_mutation_results)

            print(f"For sample {sample}: Precision: {precision:.2f}, Recall: {recall:.2f}, F1 Score: {f1:.2f}")
