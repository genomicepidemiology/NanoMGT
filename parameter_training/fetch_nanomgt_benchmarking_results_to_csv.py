import os
import csv
from collections import defaultdict
import json

def load_data_from_file(file_path):
    with open(file_path, 'r') as file:
        data = json.load(file)
    return data

def find_highest_percentage_id(batch_id, data):
    for entry in data:
        if entry['batch_id'] == batch_id:
            majority_id = list(entry['minority_abundance']['Majority'].keys())[0]
            minor_ids = list(entry['minority_abundance']['Minority'].keys())
            # Fetch the abundance of the first minor entry, assume it's representative
            if minor_ids:
                abundance = entry['minority_abundance']['Minority'][minor_ids[0]]
            else:
                abundance = 0  # Default to 0 if no minority entries exist
            return majority_id, minor_ids, abundance
    return "Batch ID not found.", [], 0

def process_datasets(base_path, output_file):
    datasets = ['clean', 'contaminated', 'mixed']
    maf_values = ['1', '2', '3', '4', '5']
    training_or_validation_extension_json = '_validation.json'

    results = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))

    for dataset in datasets:
        json_info_path = f'/home/projects/cge/people/malhal/nanomgt_json/simulated_batches_{dataset}'

        for maf_value in maf_values:
            param_path = os.path.join(base_path, f'{dataset}_benchmark', f'maf_{maf_value}')
            for species_folder in os.listdir(param_path):
                batch_id = int(species_folder.split('_')[-1])
                top_results_path = os.path.join(param_path, species_folder, 'top_result.csv')
                if os.path.exists(top_results_path):
                    print (top_results_path)
                    with open(top_results_path, mode='r') as file:
                        csv_reader = csv.DictReader(file)
                        for row in csv_reader:
                            f1 = float(row['F1 Score'])
                            if f1 != 0:
                                precision = float(row['Precision'])
                                recall = float(row['Recall'])
                                json_file = os.path.join(json_info_path,
                                                         species_folder.split('_')[0] + '_' + species_folder.split('_')[
                                                             1] + training_or_validation_extension_json)
                                data = load_data_from_file(json_file)
                                result = next((entry for entry in data if entry['batch_id'] == batch_id), None)

                                if result:
                                    major, minor_list, abundance = find_highest_percentage_id(batch_id, data)
                                    print (abundance)

                                    results[dataset][maf_value][abundance].append((f1, precision, recall))

    # Write aggregated results to CSV
    with open(output_file, mode='w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(['Dataset', 'MAF', 'Abundance Level', 'Average F1 Score', 'Average Precision', 'Average Recall'])

        for dataset, maf_dict in results.items():
            for maf, abundance_dict in maf_dict.items():
                for abundance, scores in abundance_dict.items():
                    avg_f1 = sum(x[0] for x in scores) / len(scores)
                    avg_precision = sum(x[1] for x in scores) / len(scores)
                    avg_recall = sum(x[2] for x in scores) / len(scores)
                    writer.writerow([dataset, int(maf)/100, abundance, f"{avg_f1:.4f}", f"{avg_precision:.4f}", f"{avg_recall:.4f}"])

if __name__ == "__main__":
    base_path = '/home/projects/cge/people/malhal/nanomgt_benchmark'
    output_file = 'aggregated_metrics.csv'
    process_datasets(base_path, output_file)
