import os
import csv
import json
from collections import defaultdict

def load_data_from_file(file_path):
    with open(file_path, 'r') as file:
        data = json.load(file)
    return data

def find_highest_percentage_id(batch_id, data):
    for entry in data:
        if entry['batch_id'] == batch_id:
            majority_id = list(entry['minority_abundance']['Majority'].keys())[0]
            minor_ids = list(entry['minority_abundance']['Minority'].keys())
            abundance = entry['minority_abundance']['Minority'][minor_ids[0]] if minor_ids else 0
            return majority_id, minor_ids, abundance
    return "Batch ID not found.", [], 0

def load_mutations_from_files(file_paths):
    mutations_dict = {}
    for file_path in file_paths:
        with open(file_path, 'r') as file:
            lines = file.readlines()
            for i in range(0, len(lines), 2):
                gene_id = lines[i].strip()
                mutations = set(lines[i+1].strip().split(','))
                mutations_dict[gene_id] = mutations_dict.get(gene_id, set()).union(mutations)
    return mutations_dict

def benchmark_analysis_result(major, minor_list, maps_path):
    map_files = [f'{maps_path}/major_{major}_minor_{minor_id}.txt' for minor_id in minor_list]
    return load_mutations_from_files(map_files)

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
    precision = total_tp / (total_tp + total_fp) if total_tp + total_fp > 0 else 0
    recall = total_tp / (total_tp + total_fn) if total_tp + total_fn > 0 else 0
    f1 = 2 * (precision * recall) / (precision + recall) if precision + recall > 0 else 0
    return precision, recall, f1, total_tp, total_fp, total_fn

def process_datasets(base_path, training_or_validation_extension_json, maps_path, output_file):
    datasets = ['120', '170', '220']
    bf_values = ['bf_0.01', 'bf_0.02', 'bf_0.03', 'bf_0.04', 'bf_0.05']
    results = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))  # {dataset: {MAF: {abundance: [metrics]}}}

    for dataset in datasets:
        json_info_path = '/home/projects/cge/people/malhal/nanomgt_json/simulated_batches_clean'

        for bf_value in bf_values:
            param_path = os.path.join(base_path, dataset, bf_value)
            for species_folder in os.listdir(param_path):
                if species_folder != '$base_folder':
                    species_path = os.path.join(param_path, species_folder)
                    batch_id = int(species_folder.split('_')[-1])
                    json_file = os.path.join(json_info_path, species_folder.split('_')[0] + '_' + species_folder.split('_')[1] + training_or_validation_extension_json)
                    data = load_data_from_file(json_file)
                    result = next((entry for entry in data if entry['batch_id'] == batch_id), None)

                    if result:
                        major, minor_list, abundance = find_highest_percentage_id(batch_id, data)
                        minor_mutation_expected = benchmark_analysis_result(major, minor_list, maps_path)
                        minor_mutation_results = defaultdict(set)
                        contamination_file = os.path.join(species_path, f"{species_folder}_contamination.csv")
                        if os.path.exists(contamination_file):
                            with open(contamination_file, 'r') as subfile:
                                reader = csv.DictReader(subfile)
                                for row in reader:
                                    gene = row['Gene'].split('_')[0]
                                    position = row['Position']
                                    mutation_base = row['TotalSNVs'].split(':')[0]
                                    mutation = f"{position}_{mutation_base}"
                                    minor_mutation_results[gene].add(mutation)

                        precision, recall, f1, tp, fp, fn = calculate_metrics(minor_mutation_expected, minor_mutation_results)
                        results[dataset][bf_value][abundance].append((precision, recall, f1))

    # Write aggregated results
    with open(output_file, 'w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(['Dataset', 'MAF', 'Abundance Level', 'Precision', 'Recall', 'F1 Score'])

        for dataset, maf_data in results.items():
            for maf, abundance_data in maf_data.items():
                for abundance, metrics in abundance_data.items():
                    avg_precision = sum(x[0] for x in metrics) / len(metrics)
                    avg_recall = sum(x[1] for x in metrics) / len(metrics)
                    avg_f1 = sum(x[2] for x in metrics) / len(metrics)
                    writer.writerow([dataset, maf, abundance, f"{avg_precision:.4f}", f"{avg_recall:.4f}", f"{avg_f1:.4f}"])

# Initialize paths and settings
base_path = '/home/projects/cge/people/malhal/nanomgt_depth_results'
training_or_validation_extension_json = '_validation.json'
maps_path = '/home/projects/cge/people/malhal/nanomgt_reads/variant_maps/'
output_file = 'average_metrics.csv'

# Process all datasets and bf values
process_datasets(base_path, training_or_validation_extension_json, maps_path, output_file)
