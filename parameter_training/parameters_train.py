import os
import sys
import json
import pandas as pd
import numpy as np
import csv
import re
import concurrent.futures
from collections import defaultdict
from scipy.interpolate import UnivariateSpline
from itertools import product

sys.path = [os.path.join(os.path.dirname(os.path.realpath(__file__)), '..')] + sys.path

from nanomgt import nanopore_variantcaller as nvc

# Modify these paths as needed.
alignment_results_path = '/home/projects/cge/people/malhal/nanomgt_new_results/clean/'
maps_path = '/home/projects/cge/people/malhal/nanomgt_reads/variant_maps/'
json_info_path = '/home/projects/cge/people/malhal/nanomgt_json/simulated_batches_clean'
training_or_validation_extension_json = '_training.json'
output_training_folder = 'clean_test_for_final_script_only_5_maf'
param_list = ['np', 'cor', 'pp', 'dp', 'ii']
maf_interval = [5]
cpus = 39
model_name = 'model_test'

def main():
    parameters_interval_search = {
        'cor_interval': [0.5],
        'ii_interval': [0.1, 0.2],
        'pp_interval': [0.2],
        'np_interval': [3],
        'dp_interval': [0.1]
    }
    folders = [f for f in os.listdir(alignment_results_path) if os.path.isdir(os.path.join(alignment_results_path, f))]
    run_parameter_search(folders, maf_interval, parameters_interval_search, output_training_folder)
    generate_spline_json(output_training_folder, maf_interval, model_name)

def generate_spline_json(output_folder, maf_intervals, model_name):
    parameter_files = [f"{output_folder}/5_round/maf_{maf}_average_best_params.json" for maf in maf_intervals]
    parameters_list = []

    for file_path in parameter_files:
        with open(file_path, 'r') as f:
            params = json.load(f)
            parameters_list.append(params)

    cor_intervals = []
    iteration_intervals = []
    pp_intervals = []
    np_intervals = []
    dp_intervals = []

    for params in parameters_list:
        cor_intervals.append(params.get('cor', 0))
        iteration_intervals.append(params.get('ii', 0))
        pp_intervals.append(params.get('pp', 0))
        np_intervals.append(params.get('np', 0))
        dp_intervals.append(params.get('dp', 0))

    x_values = np.linspace(0.01, 0.05, len(parameters_list))
    fine_x_values = np.linspace(0.01, 0.05, 500)
    spline_results = {}

    for name, data in zip(['cor', 'ii', 'pp', 'np', 'dp'],
                          [cor_intervals, iteration_intervals, pp_intervals, np_intervals, dp_intervals]):
        spline = UnivariateSpline(x_values, data, s=None)
        spline_fit_fine = spline(fine_x_values)
        spline_results[name] = {str(x): float(y) for x, y in zip(fine_x_values, spline_fit_fine)}

    output_file = f'{model_name}.json'
    with open(output_file, 'w') as f:
        json.dump(spline_results, f, indent=4)

    print(f"Spline coefficients have been saved to {output_file}.")
def train_parameters(maf, results_folder, min_n, cor, new_output_folder, maps_path, json_info_path,
                     ii, proxi, dp_window, pp, np, dp, consensus_dict, bio_validation_dict, minor_mutation_expected):
    confirmed_mutation_dict = nvc.derive_mutation_positions(consensus_dict, min_n, maf, cor)

    confirmed_mutation_dict, co_occurrence_tmp_dict, iteration_count, mutation_threshold_dict =\
        nvc.snv_convergence(results_folder, maf, cor, np, pp, dp, proxi, dp_window, ii,
                            confirmed_mutation_dict, consensus_dict, bio_validation_dict, min_n)

    parameter_string = f"maf_{maf}_cor_{cor}_pp_{pp}_np_{np}_dp_{dp}_ii_{ii}"

    nvc.format_output(new_output_folder, confirmed_mutation_dict, consensus_dict, bio_validation_dict,
                      co_occurrence_tmp_dict, mutation_threshold_dict, parameter_string)

    minor_mutation_results = convert_mutation_dict_to_object(confirmed_mutation_dict)

    precision, recall, f1, tp, fp, fn = calculate_metrics(minor_mutation_expected, minor_mutation_results)

    return f1, parameter_string, precision, recall, tp, fp, fn

def load_data(filepath):
    return pd.read_csv(filepath, skipinitialspace=True)

def find_highest_percentage_id(batch_id, data):
    for entry in data:
        if entry['batch_id'] == batch_id:
            majority_id = list(entry['minority_abundance']['Majority'].keys())[0]
            minor_ids = list(entry['minority_abundance']['Minority'].keys())
            return majority_id, minor_ids
    return "Batch ID not found.", []

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

def benchmark_analysis_result(sample, json_file_path, maps_path, training_or_validation_extension_json):
    batch_id = int(sample.split('_')[-1])
    species = sample.split('_')[0] + '_' + sample.split('_')[1]

    # Load JSON data
    with open(json_file_path + '/' + species + training_or_validation_extension_json, 'r') as file:
        data = json.load(file)

    # Get the highest percentage ID and minor IDs
    top_id, minor = find_highest_percentage_id(batch_id, data)

    if top_id == "Batch ID not found.":
        return top_id

    map_files = []
    for minor_id in minor:
        map_file = f'{maps_path}/major_{top_id}_minor_{minor_id}.txt'
        map_files.append(map_file)

    mutation_map = load_mutations_from_files(map_files)

    return mutation_map

def convert_mutation_dict_to_object(mutation_dict):
    mutation_object = {}
    for allele in mutation_dict:
        gene = allele.split('_')[0]
        if gene not in mutation_object:
            mutation_object[gene] = set()
        for mutation in mutation_dict[allele][0]:
            mutation_object[gene].add(mutation)
    return mutation_object

def calculate_metrics(expected_mutations, actual_mutations):
    total_tp = 0
    total_fp = 0
    total_fn = 0

    for gene, expected_set in expected_mutations.items():
        if gene in actual_mutations:
            actual_set = actual_mutations[gene]
            tp = len(expected_set & actual_set)
            fp = len(actual_set - expected_set)
            fn = len(expected_set - actual_set)

            total_tp += tp
            total_fp += fp
            total_fn += fn

    precision = total_tp / (total_tp + total_fp) if (total_tp + total_fp) > 0 else 0
    recall = total_tp / (total_tp + total_fn) if (total_tp + total_fn) > 0 else 0
    f1 = 2 * (precision * recall) / (precision + recall) if (precision + recall) > 0 else 0

    return precision, recall, f1, total_tp, total_fp, total_fn

def run_jobs_in_parallel(max_workers, new_output_folder, alignment_folder, maf, parameters, maps_path, json_info_path, training_or_validation_extension_json):
    min_n = 3
    proxi = 5
    dp_window = 15

    cor_interval = parameters['cor_interval']
    ii_interval = parameters['ii_interval']
    pp_interval = parameters['pp_interval']
    np_interval = parameters['np_interval']
    dp_interval = parameters['dp_interval']

    best_score = 0
    best_params = None
    top_precision = None
    top_recall = None
    top_tp = None
    top_fp = None
    top_fn = None

    all_params = list(product(cor_interval, ii_interval, pp_interval, np_interval, dp_interval))
    total_combinations = len(all_params)
    print(f"Total number of parameter combinations: {total_combinations}")

    results_filename = new_output_folder + "/all_results.csv"
    top_result_filename = new_output_folder + "/top_result.csv"

    all_results = []

    processed_combinations = 0

    consensus_dict = nvc.build_consensus_dict(os.path.join(alignment_folder, 'rmlst_alignment.res'),
                                              os.path.join(alignment_folder, 'rmlst_alignment.mat'))
    bio_validation_dict = nvc.bio_validation_mutations(consensus_dict, os.path.join(alignment_folder, 'specie.fsa'))

    sample = alignment_folder.split('/')[-1]

    minor_mutation_expected = benchmark_analysis_result(sample, json_info_path, maps_path, training_or_validation_extension_json)

    with concurrent.futures.ProcessPoolExecutor(max_workers=max_workers) as executor:
        futures_to_params = {
            executor.submit(
                train_parameters, maf, alignment_folder, min_n,
                combo[0], new_output_folder, maps_path, json_info_path, combo[1], proxi, dp_window,
                combo[2], combo[3], combo[4], consensus_dict, bio_validation_dict, minor_mutation_expected
            ): combo for combo in all_params
        }

        for future in concurrent.futures.as_completed(futures_to_params):
            params = futures_to_params[future]
            processed_combinations += 1
            print(f"Processed {processed_combinations}/{total_combinations} combinations.")
            try:
                result = future.result()
                f1, parameter_string, precision, recall, tp, fp, fn = result
                all_results.append([f1, parameter_string, precision, recall, tp, fp, fn])

                if f1 > best_score:
                    best_score = f1
                    best_params = parameter_string
                    top_precision = precision
                    top_recall = recall
                    top_tp = tp
                    top_fp = fp
                    top_fn = fn
            except Exception as exc:
                print(f"Generated an exception: {exc}")

    with open(results_filename, 'w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(['F1 Score', 'Parameters', 'Precision', 'Recall', 'TP', 'FP', 'FN'])
        writer.writerows(all_results)

    with open(top_result_filename, 'w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(['F1 Score', 'Parameters', 'Precision', 'Recall', 'TP', 'FP', 'FN'])
        writer.writerow([best_score, best_params, top_precision, top_recall, top_tp, top_fp, top_fn])

    if best_score != 0:
        for file in os.listdir(new_output_folder):
            if file.endswith('minor_mutations.csv'):
                if not file.startswith(best_params):
                    os.system('rm {}/{}'.format(new_output_folder, file))

def extract_parameters(param_string):
    param_pattern = r'([a-z_]+)_([0-9.]+)'
    return {key: float(value) for key, value in re.findall(param_pattern, param_string)}

def calculate_best_parameters(file_name):
    df = pd.read_csv(file_name)

    parameters_data = defaultdict(lambda: {'f1_scores': [], 'best_f1_score': float('-inf'), 'best_params': {}})

    for index, row in df.iterrows():
        f1_score = row['F1 Score']
        parameters = extract_parameters(row['Parameters'])

        if f1_score is None or f1_score == 0:
            continue  # Skip this iteration if F1 Score is None or 0

        key = frozenset(parameters.items())
        parameters_data[key]['f1_scores'].append(f1_score)

        if f1_score > parameters_data[key]['best_f1_score']:
            parameters_data[key]['best_f1_score'] = f1_score
            parameters_data[key]['best_params'] = parameters

    if not parameters_data:
        print("No data to analyze.")
        return None

    top_f1_score = max(data['best_f1_score'] for data in parameters_data.values())
    if top_f1_score is None or top_f1_score == 0:
        return None  # Return None if the best F1-score is None or 0

    top_parameters = [data['best_params'] for data in parameters_data.values() if data['best_f1_score'] == top_f1_score]

    if len(top_parameters) == 1:
        return top_parameters[0]

    penalty_order = ['np', 'pp', 'dp', 'cor', 'ii']

    def sort_key(params):
        return tuple(params.get(param, float('inf')) for param in penalty_order)

    sorted_top_parameters = sorted(top_parameters, key=sort_key)

    return sorted_top_parameters[0]

def load_default_parameters(file_path):
    with open(file_path, 'r') as file:
        return json.load(file)

def generate_test_values(average_value, num_increments, percentage):
    increment = average_value * percentage
    test_values = []
    for i in range(-num_increments, num_increments + 1):
        test_values.append(average_value + i * increment)
    return test_values

def create_test_object(default_params, param, test_values):
    test_objects = []
    for value in test_values:
        new_params = default_params.copy()
        new_params[param] = value
        test_objects.append(new_params)
    return test_objects

def determine_gradient_value(df, param):
    results = []
    grouped = df.groupby('MAF')
    for maf, group in grouped:
        param_values = group['Parameter Value'].values
        f1_scores = group['F1 Score'].values

        param_f1_map = defaultdict(list)
        for param_value, f1_score in zip(param_values, f1_scores):
            param_f1_map[param_value].append(f1_score)

        param_values_new = np.array(list(param_f1_map.keys()))
        f1_scores_new = np.array([np.mean(f1_scores) for f1_scores in param_f1_map.values()])

        if len(param_values_new) < 2 or len(f1_scores_new) < 2:
            print(f"Skipping {param}, maf {maf} due to insufficient data points")
            continue
        param_values_range = np.max(param_values_new) - np.min(param_values_new)
        f1_scores_range = np.max(f1_scores_new) - np.min(f1_scores_new)

        if param_values_range == 0:
            print(f"Skipping {param}, maf {maf} due to zero parameter value range")
            continue

        if f1_scores_range == 0:
            param_value_to_return = np.mean(param_values_new)
            results.append({
                'maf': maf,
                'param': param,
                'param_value_to_return': param_value_to_return,
                'first_f1_score': f1_scores_new[0],
                'last_f1_score': f1_scores_new[-1],
                'lowest_slope_angle': 0
            })
            continue

        f1_scores_first_value = f1_scores_new[0]
        f1_scores_normalized = f1_scores_new - f1_scores_first_value
        spline = UnivariateSpline(param_values_new, f1_scores_normalized, s=None)
        derivative = spline.derivative()
        param_dense = np.linspace(min(param_values_new), max(param_values_new), 50)
        f1_dense = spline(param_dense)

        derivative_values = derivative(param_dense)

        trend = f1_dense[-1] - f1_dense[0]
        print(f"Processing param: {param}, maf: {maf}, trend: {trend}")

        max_derivative_value = np.max(derivative_values)
        peak_index = np.argmax(derivative_values)
        target_value = 0.8 * max_derivative_value

        param_value_to_return = param_dense[peak_index]
        for idx in range(peak_index + 1, len(derivative_values)):
            if derivative_values[idx] <= target_value:
                param_value_to_return = param_dense[idx]
                break

        min_slope_angle = np.degrees(np.arctan(np.min(derivative_values)))
        first_f1_score = f1_dense[0]
        last_f1_score = f1_dense[-1]

        print(f"Processed param: {param}, maf: {maf}, param_value_to_return: {param_value_to_return}")

        results.append({
            'maf': maf,
            'param': param,
            'param_value_to_return': param_value_to_return,
            'first_f1_score': first_f1_score,
            'last_f1_score': last_f1_score,
            'lowest_slope_angle': min_slope_angle
        })
    return results

def process_total_parameter_results(total_parameter_results):
    result_dict = {}
    for param, maf_data in total_parameter_results.items():
        for maf, batch_data in maf_data.items():
            df_rows = []
            for batch_id, (param_values, f1_scores) in batch_data.items():
                for param_value, f1_score in zip(param_values, f1_scores):
                    df_rows.append(
                        {'MAF': maf, 'Batch ID': batch_id, 'F1 Score': f1_score, 'Parameter Value': param_value})
            df = pd.DataFrame(df_rows)

            processed_results = determine_gradient_value(df, param)

            if maf not in result_dict:
                result_dict[maf] = {}
            for result in processed_results:
                result_dict[maf][param] = result['param_value_to_return']

    return result_dict

def load_results(param_list, maf_interval, output_training_folder):
    total_parameter_results = {}
    for param in param_list:
        total_parameter_results[param] = {}
        for maf in maf_interval:
            total_parameter_results[param][maf] = {}
            for folder in os.listdir(os.path.join(output_training_folder, "maf_{}".format(maf))):
                if folder.startswith(param):
                    batch_id = int(folder.split('_')[-1])
                    if batch_id > 10:
                        batch_id = batch_id - 10
                    total_parameter_results[param][maf][batch_id] = [[], []]
                    results_file = os.path.join(output_training_folder, "maf_{}".format(maf), folder, 'all_results.csv')

                    with open(results_file, 'r') as csvfile:
                        reader = csv.DictReader(csvfile)
                        for row in reader:
                            parameter_value = extract_param_value(row['Parameters'], param)
                            f1_score = float(row['F1 Score'])

                            if parameter_value is not None:
                                total_parameter_results[param][maf][batch_id][0].append(parameter_value)
                                total_parameter_results[param][maf][batch_id][1].append(f1_score)

                    total_parameter_results[param][maf][batch_id][0].sort()
                    total_parameter_results[param][maf][batch_id][1].sort()

    average_parameter_results = {}
    for param in param_list:
        average_parameter_results[param] = {}
        for maf in maf_interval:
            average_parameter_results[param][maf] = {}
            for batch_id, values in total_parameter_results[param][maf].items():
                param_values = values[0]
                f1_scores = values[1]

                param_f1_map = defaultdict(list)
                for param_value, f1_score in zip(param_values, f1_scores):
                    param_f1_map[param_value].append(f1_score)

                unique_param_values = []
                average_f1_scores = []
                for param_value, f1_scores in param_f1_map.items():
                    unique_param_values.append(param_value)
                    average_f1_scores.append(sum(f1_scores) / len(f1_scores))

                average_parameter_results[param][maf][batch_id] = [unique_param_values, average_f1_scores]

    return average_parameter_results

def extract_param_value(parameter_string, param):
    param_list = parameter_string.split('_')
    for i, p in enumerate(param_list):
        if p == param and i + 1 < len(param_list):
            try:
                return float(param_list[i + 1])
            except ValueError:
                return None
    return None

def process_directory(directory):
    result_dict = {}

    for filename in os.listdir(directory):
        if filename.endswith(".csv"):
            parts = filename.split('_')
            param = parts[0]
            maf = parts[1]

            file_path = os.path.join(directory, filename)
            df = pd.read_csv(file_path)

            processed_results = determine_gradient_value(df, param)

            if maf not in result_dict:
                result_dict[maf] = {}

            for result in processed_results:
                result_dict[maf][param] = result['param_value_to_return']

    return result_dict

def setup_folders(maf_interval, alignment_results_path):
    folders = [f for f in os.listdir(alignment_results_path) if os.path.isdir(os.path.join(alignment_results_path, f))]
    os.makedirs(output_training_folder, exist_ok=True)
    for maf in maf_interval:
        os.makedirs(os.path.join(output_training_folder, f'maf_{maf}'), exist_ok=True)
    return folders

def setup_directory(base_path, subfolder):
    path = os.path.join(base_path, subfolder)
    os.makedirs(path, exist_ok=True)
    return path

def run_round_of_parameter_search(round_number, maf_interval, folders, output_training_folder,
                                  parameters_interval_search):
    for maf in maf_interval:
        round_folder = f"{round_number}_round_maf_{maf}"
        round_path = setup_directory(output_training_folder, round_folder)
        for folder in folders:
            batch_id = int(folder.split('_')[-1])
            if batch_id > 10:
                batch_id -= 10
            if batch_id >= maf:
                alignment_folder = os.path.join(alignment_results_path, folder)
                new_output_folder = setup_directory(round_path, folder)
                run_jobs_in_parallel(cpus, new_output_folder, alignment_folder, maf / 100,
                                     parameters_interval_search, maps_path, json_info_path,
                                     training_or_validation_extension_json)

def collect_and_store_best_params(maf_interval, folders, output_training_folder, param_list):
    all_best_params = defaultdict(list)
    for maf in maf_interval:
        average_best_params = {}
        for folder in folders:
            batch_id = int(folder.split('_')[-1])
            if batch_id > 10:
                batch_id -= 10
            if batch_id >= maf:
                results_filename = os.path.join(output_training_folder, f"maf_{maf}", folder, "all_results.csv")
                best_params = calculate_best_parameters(results_filename)
                if best_params is not None:
                    for param, value in best_params.items():
                        if param[1:] in param_list:
                            all_best_params[param[1:]].append(value)
        for param in param_list:
            if param in all_best_params:
                average_value = sum(all_best_params[param]) / len(all_best_params[param])
                average_best_params[param] = average_value
        json_path = os.path.join(output_training_folder, f"maf_{maf}_average_best_params.json")
        with open(json_path, 'w') as json_file:
            json.dump(average_best_params, json_file, indent=4)

def run_parameter_search(folders, maf_interval, parameters_interval_search, output_training_folder):
    for maf in maf_interval:
        maf_folder = f"maf_{maf}"
        maf_path = setup_directory(output_training_folder, maf_folder)
        for folder in folders:
            batch_id = int(folder.split('_')[-1])
            if batch_id > 10:
                batch_id -= 10
            if batch_id >= maf:
                alignment_folder = os.path.join(alignment_results_path, folder)
                new_output_folder = setup_directory(maf_path, folder)
                run_jobs_in_parallel(cpus, new_output_folder, alignment_folder, maf / 100,
                                     parameters_interval_search, maps_path, json_info_path,
                                     training_or_validation_extension_json)

    collect_and_store_best_params(maf_interval, folders, output_training_folder, param_list)

    num_increments = 1
    rounds = [2, 3, 4, 5]
    round_increment_dict = {2: 0.20, 3: 0.15, 4: 0.10, 5: 0.05}
    for round_number in rounds:
        print(f'Starting round {round_number}')
        for maf in maf_interval:
            if round_number == 2:
                previous_round_path = os.path.join(output_training_folder, f"maf_{maf}_average_best_params.json")
            else:
                previous_round_path = os.path.join(output_training_folder,
                                                   f"{round_number - 1}_round_maf_{maf}_average_best_params.json")

            default_params = load_default_parameters(previous_round_path)
            for param, default_value in default_params.items():
                test_values = generate_test_values(default_value, num_increments, round_increment_dict[round_number])
                parameters_interval_search[param + '_interval'] = test_values

            round_folder = f"{round_number}_round_maf_{maf}"
            round_path = setup_directory(output_training_folder, round_folder)
            run_round_of_parameter_search(round_number, [maf], folders, round_path, parameters_interval_search)

        round_path = f"{round_number}_round"
        collect_and_store_best_params(maf_interval, folders, setup_directory(output_training_folder, round_path),
                                      param_list)

if __name__ == "__main__":
    main()
