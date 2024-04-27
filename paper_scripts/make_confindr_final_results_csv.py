import os
import sys
import pandas as pd

def load_mutations_from_files(file_paths):
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

def calculate_metrics(expected_mutations, actual_mutations):
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


def load_data(filepath):
    # Load data with handling for initial spaces in column names
    return pd.read_csv(filepath, skipinitialspace=True)

def load_mutations_from_files(file_paths):
    mutations_dict = {}
    for file_path in file_paths:
        with open(file_path, 'r') as file:
            lines = file.readlines()
            for i in range(0, len(lines), 2):
                gene_id = lines[i].strip()
                mutations = set(lines[i+1].strip().split(','))
                if gene_id in mutations_dict:
                    mutations_dict[gene_id].update(mutations)
                else:
                    mutations_dict[gene_id] = mutations
    return mutations_dict

def find_highest_percentage_id(batch_id, df):
    batch_data = df[df['Batch'] == batch_id]
    if batch_data.empty:
        return "Batch ID not found."
    highest_percentage = 0
    highest_percentage_id = ''
    all = []
    for index, row in batch_data.iterrows():
        num_percentages = 3 if 'Percentage3' in df.columns else 2
        for i in range(1, num_percentages+1):
            percentage = int(row[f'Percentage{i}'][:-1])
            if percentage > highest_percentage:
                highest_percentage = percentage
                highest_percentage_id = row[f'ID{i}']
            all.append(row[f'ID{i}'])
    minor = [item for item in all if item != highest_percentage_id]
    return highest_percentage_id, minor

def benchmark_analysis_result(sample, results_folder):
    batch_id = int(sample.split('_')[-2][5:])
    specie = '_'.join(sample.split('_')[1:-5])
    batch_csv = '_'.join(sample.split('_')[1:-2]) + '.csv'
    data = load_data(os.path.join(results_folder, batch_csv))
    top_id, minor = find_highest_percentage_id(batch_id, data)

    mutation_map = {}
    if 'salmonella_enterica' in sample:
        files = [f"/home/people/malhal/data/new_nanomgt/majority_variants/minor_variant_maps/major_{top_id}_minor_{m}.txt" for m in minor]
        mutation_map = load_mutations_from_files(files)
    elif 'ecoli' in sample:
        files = [f"/home/people/malhal/data/new_nanomgt/majority_variants/minor_variant_maps/major_{top_id}_minor_{m}.txt" for m in minor]
        mutation_map = load_mutations_from_files(files)
    elif 'staph_aureus' in sample:
        file = f"/home/people/malhal/data/new_nanomgt/majority_variants/minor_variant_maps/major_{top_id}_minor_{minor[0]}.txt"
        mutation_map = load_mutations_from_files([file])
    elif 'campylobacter_jejuni' in sample:
        file = f"/home/people/malhal/data/new_nanomgt/majority_variants/minor_variant_maps/major_{top_id}_minor_{minor[0]}.txt"
        mutation_map = load_mutations_from_files([file])
    elif 'klebsiella_pneumoniae' in sample:
        file = f"/home/people/malhal/data/new_nanomgt/majority_variants/minor_variant_maps/major_{top_id}_minor_{minor[0]}.txt"
        mutation_map = load_mutations_from_files([file])

    return mutation_map


def main(results_folder):
    bf_folders = [folder for folder in os.listdir(results_folder) if folder.startswith('bf_')]

    all_results = []

    for bf_folder in bf_folders:
        mrd_path = os.path.join(results_folder, bf_folder)
        experiments = os.listdir(mrd_path)

        for experiment in experiments:
            exp_path = os.path.join(mrd_path, experiment)
            depth = experiment.split('_')[0].replace('depth', '')  # Extract depth from experiment name
            species = experiment.split('_')[1]
            batch = int(experiment.split('batch')[-1].split('_')[0])

            print(f"Processing experiment {experiment} in {bf_folder}...")

            batch_path = "/home/people/malhal/data/new_nanomgt/simulated_batches"
            expected_mutations = benchmark_analysis_result(experiment, batch_path)

            observed_file = [os.path.join(exp_path, f) for f in os.listdir(exp_path) if 'contamination.csv' in f][0]
            observed_mutations = {}
            data = pd.read_csv(observed_file)
            for index, row in data.iterrows():
                gene = row['Gene'].split('_')[0]
                position = row['Position']
                mutations = row['CongruentSNVs'] if pd.notna(row['CongruentSNVs']) else row['TotalSNVs']
                if pd.notna(mutations):
                    mutation_details = mutations.split(':')
                    mutation = f"{position}_{mutation_details[0]}"
                    if gene in observed_mutations:
                        observed_mutations[gene].add(mutation)
                    else:
                        observed_mutations[gene] = {mutation}

            precision, recall, f1 = calculate_metrics(expected_mutations, observed_mutations)

            all_results.append({
                'tool': 'Confindr',
                'specie': species,
                'mrd': bf_folder.split('_')[1],
                'depth': depth,
                'batch': batch,
                'precision': precision,
                'recall': recall,
                'f1score': f1
            })

            print(f"Completed processing for {experiment}. Metrics: Precision={precision}, Recall={recall}, F1 Score={f1}")

    df = pd.DataFrame(all_results)
    grouped = df.groupby(['tool', 'specie', 'mrd', 'depth', 'batch'])
    mean_df = grouped[['precision', 'recall', 'f1score']].mean().reset_index()
    mean_df.to_csv('output.csv', index=False)
if __name__ == "__main__":
    results_folder_path = "/home/people/malhal/data/new_nanomgt/confindr_results"
    main(results_folder_path)
