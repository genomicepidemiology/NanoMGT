import os
import pandas as pd

def load_data(filepath):
    """ Load CSV data, handling spaces in column names. """
    return pd.read_csv(filepath, skipinitialspace=True)

def load_mutations_from_files(file_paths):
    """ Load mutations from files and return a dictionary of gene IDs to mutations. """
    mutations_dict = {}
    for file_path in file_paths:
        with open(file_path, 'r') as file:
            lines = file.readlines()
            for i in range(0, len(lines), 2):
                gene_id = lines[i].strip()
                mutations = set(lines[i+1].strip().split(','))
                mutations_dict.setdefault(gene_id, set()).update(mutations)
    return mutations_dict

def calculate_metrics(expected_mutations, actual_mutations):
    """ Calculate precision, recall, and F1-score based on expected and actual mutations. """
    total_tp = total_fp = total_fn = 0
    for gene, expected_set in expected_mutations.items():
        actual_set = actual_mutations.get(gene, set())
        tp = len(expected_set & actual_set)
        fp = len(actual_set - expected_set)
        fn = len(expected_set - actual_set)
        total_tp += tp
        total_fp += fp
        total_fn += fn
    precision = total_tp / (total_tp + total_fp) if (total_tp + total_fp) > 0 else 0
    recall = total_tp / (total_tp + total_fn) if (total_tp + total_fn) > 0 else 0
    f1 = 2 * (precision * recall) / (precision + recall) if (precision + recall) > 0 else 0
    return precision, recall, f1

def find_highest_percentage_id(batch_id, df):
    """ Find the highest percentage ID from a DataFrame based on the batch ID. """
    batch_data = df[df['Batch'] == batch_id]
    if batch_data.empty:
        return "Batch ID not found.", []
    highest_percentage = 0
    highest_percentage_id = ''
    all_ids = []
    for index, row in batch_data.iterrows():
        num_percentages = 3 if 'Percentage3' in df.columns else 2
        for i in range(1, num_percentages + 1):
            percentage = int(row[f'Percentage{i}'][:-1])
            all_ids.append(row[f'ID{i}'])
            if percentage > highest_percentage:
                highest_percentage = percentage
                highest_percentage_id = row[f'ID{i}']
    minor_ids = [item for item in all_ids if item != highest_percentage_id]
    return highest_percentage_id, minor_ids

def benchmark_analysis_result(sample, results_folder):
    """ Analyze benchmark results and extract mutations. """
    parts = sample.split('_')
    batch_id = int(parts[-2][5:])
    species = '_'.join(parts[1:-5])
    sequencing_id = parts[2] if 'ecoli' in sample else parts[3]
    batch_csv = '_'.join(parts[1:-2]) + '.csv'
    data = load_data(os.path.join(results_folder, batch_csv))
    top_id, minor = find_highest_percentage_id(batch_id, data)

    # Load mutation files
    file_template = f"/home/people/malhal/data/new_nanomgt/majority_variants/minor_variant_maps/major_{top_id}_minor_{{}}.txt"
    mutation_files = [file_template.format(m) for m in minor]
    return species, sequencing_id, load_mutations_from_files(mutation_files)

def main(results_folder):
    """ Main function to process all results. """
    bf_folders = [folder for folder in os.listdir(results_folder) if folder.startswith('bf_')]

    all_results = []

    for bf_folder in bf_folders:
        maf_path = os.path.join(results_folder, bf_folder)
        experiments = os.listdir(maf_path)

        for experiment in experiments:
            exp_path = os.path.join(maf_path, experiment)
            print(f"Processing experiment {experiment} in {bf_folder}...")

            batch_path = "/home/people/malhal/data/new_nanomgt/simulated_batches"
            species, sequencing_id, expected_mutations = benchmark_analysis_result(experiment, batch_path)

            observed_file = next((os.path.join(exp_path, f) for f in os.listdir(exp_path) if 'contamination.csv' in f), None)
            observed_mutations = {}

            if observed_file:
                data = pd.read_csv(observed_file)
                for index, row in data.iterrows():
                    gene = row['Gene'].split('_')[0]
                    position = row['Position']
                    mutations = row['CongruentSNVs'] if pd.notna(row['CongruentSNVs']) else row['TotalSNVs']
                    if pd.notna(mutations):
                        mutation_details = mutations.split(':')
                        mutation = f"{position}_{mutation_details[0]}"
                        observed_mutations.setdefault(gene, set()).add(mutation)

            precision, recall, f1 = calculate_metrics(expected_mutations, observed_mutations)

            all_results.append({
                'tool': 'Confindr',
                'specie': species,
                'sequencing_id': sequencing_id,
                'maf': bf_folder.split('_')[1],
                'depth': experiment.split('_')[0].replace('depth', ''),
                'batch': int(experiment.split('batch')[-1].split('_')[0]),
                'precision': precision,
                'recall': recall,
                'f1score': f1
            })

            print(f"Completed processing for {experiment}. Metrics: Precision={precision}, Recall={recall}, F1 Score={f1}")

    # Output results directly to CSV without averaging
    df = pd.DataFrame(all_results)
    df.to_csv('confindr_final_individual_results.csv', index=False)

if __name__ == "__main__":
    results_folder_path = "/home/people/malhal/data/new_nanomgt/confindr_results"
    main(results_folder_path)
