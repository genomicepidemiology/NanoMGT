import os
import sys
import gzip
import time
import pandas as pd
import numpy as np
import csv
import re
import multiprocessing
import json
from Bio import SeqIO
from itertools import combinations
import concurrent.futures
from itertools import product
import argparse
from collections import defaultdict
from scipy.interpolate import UnivariateSpline

sys.path = [os.path.join(os.path.dirname(os.path.realpath(__file__)), '..')] + sys.path

from nanomgt import nanopore_variantcaller as nvc

alignment_results_path = '/home/people/malhal/test/training_test/'
maps_path = '/home/people/malhal/test/training_test/maps/'
simulated_batches_csv_path = '/home/people/malhal/test/training_test/data/simulated_batches/'
files = os.listdir(alignment_results_path)
folders = [f for f in os.listdir(alignment_results_path)]

maf_interval = [2]

#cor_interval_search = [0.3, 0.4, 0.5, 0.6]
#dp_interval_search = [0.1, 0.3, 0.5]
#np_interval_search = [0.1, 1, 2, 3, 4]
#pp_interval_search = [0.1, 0.2, 0.3, 0.4, 0.5]
#ii_interval_search = [0.01, 0.1, 0.3, 0.5]

cor_interval_search = [0.3]
dp_interval_search = [0.35]
np_interval_search = [3]
pp_interval_search = [0.15]
ii_interval_search = [0.5]

parameters_interval_search = {
    'cor_interval': cor_interval_search,
    'ii_interval': ii_interval_search,
    'pp_interval': pp_interval_search,
    'np_interval': np_interval_search,
    'dp_interval': dp_interval_search
}

cpus = 40

def train_parameters(maf, results_folder, min_n, cor, new_output_folder, maps_path, simulated_batches_csv_path,
                    ii, proxi, dp_window, pp, np, dp):
    consensus_dict = nvc.build_consensus_dict(os.path.join(results_folder, 'rmlst_alignment.res'),
                                          os.path.join(results_folder, 'rmlst_alignment.mat'))

    confirmed_mutation_dict = nvc.derive_mutation_positions(consensus_dict, min_n, maf, cor)

    bio_validation_dict = nvc.bio_validation_mutations(consensus_dict, os.path.join(results_folder, 'specie.fsa'))

    confirmed_mutation_dict, co_occurrence_tmp_dict, iteration_count, mutation_threshold_dict =\
        nvc.snv_convergence(results_folder, maf, cor, np, pp, dp, proxi, dp_window, ii,
                            confirmed_mutation_dict, consensus_dict, bio_validation_dict)

    nvc.format_output(new_output_folder, confirmed_mutation_dict, consensus_dict, bio_validation_dict,
                  co_occurrence_tmp_dict, mutation_threshold_dict)

    sample = results_folder.split('/')[-1]

    minor_mutation_expected = benchmark_analysis_result(sample, simulated_batches_csv_path, maps_path)

    minor_mutation_results = convert_mutation_dict_to_object(confirmed_mutation_dict)

    print(len(minor_mutation_expected), len(minor_mutation_results))

    precision, recall, f1, tp, fp, fn = calculate_metrics(minor_mutation_expected, minor_mutation_results)

    parameter_string = f"maf_{maf}_cor_{cor}_pp_{pp}_np_{np}_dp_{dp}_ii_{ii}"

    return f1, parameter_string, precision, recall, tp, fp, fn

def load_data(filepath):
    return pd.read_csv(filepath, skipinitialspace=True)

def find_highest_percentage_id(batch_id, df):
    batch_data = df[df['Batch'] == batch_id]
    if batch_data.empty:
        return "Batch ID not found."
    highest_percentage = 0
    highest_percentage_id = ''
    all = []
    for index, row in batch_data.iterrows():
        if 'Percentage3' in df.columns:
            for i in range(1, 4):
                percentage = int(row[f'Percentage{i}'][:-1])
                if percentage > highest_percentage:
                    highest_percentage = percentage
                    highest_percentage_id = row[f'ID{i}']
                all.append((row[f'ID{i}']))
        else:
            for i in range(1, 3):
                percentage = int(row[f'Percentage{i}'][:-1])
                if percentage > highest_percentage:
                    highest_percentage = percentage
                    highest_percentage_id = row[f'ID{i}']
                all.append((row[f'ID{i}']))

    minor = []

    for item in all:
        if item != highest_percentage_id:
            minor.append(item)

    return highest_percentage_id, minor

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

def get_number_of_columns(dataframe):
    return dataframe.shape[1]

def benchmark_analysis_result(sample, batch_csv_path, maps_path):
    batch_id = int(sample.split('_')[-2][5:])
    batch_csv = batch_csv_path + "_".join(sample.split('_')[1:-2]) + ".csv"
    data = load_data(batch_csv)
    sample_number = int((get_number_of_columns(data) - 1) / 2)

    top_id, minor = find_highest_percentage_id(batch_id, data)

    map_files = []
    for i in range(sample_number - 1):
        map_file = f'{maps_path}/major_{top_id}_minor_{minor[i]}.txt'
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

def run_jobs_in_parallel(max_workers, new_output_folder, alignment_folder, maf, parameters, maps_path, simulated_batches_csv_path):
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

    with concurrent.futures.ProcessPoolExecutor(max_workers=max_workers) as executor:
        futures_to_params = {
            executor.submit(
                train_parameters, maf, alignment_folder, min_n,
                combo[0], new_output_folder, maps_path, simulated_batches_csv_path, combo[1], proxi, dp_window,
                combo[2], combo[3], combo[4]
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

def extract_parameters(param_string):
    param_pattern = r'([a-z_]+)_([0-9.]+)'
    return dict(re.findall(param_pattern, param_string))

def calculate_best_parameters(file_name):
    df = pd.read_csv(file_name)

    parameters_data = defaultdict(lambda: {'f1_scores': [], 'best_f1_score': float('-inf'), 'best_params': {}})

    for index, row in df.iterrows():
        f1_score = row['F1 Score']
        parameters = extract_parameters(row['Parameters'])

        for param, value in parameters.items():
            key = f"{param}_{value}"
            parameters_data[key]['f1_scores'].append(f1_score)

            if f1_score > parameters_data[key]['best_f1_score']:
                parameters_data[key]['best_f1_score'] = f1_score
                parameters_data[key]['best_params'] = {param: float(value)}

    best_params_per_file = {}
    for key, data in parameters_data.items():
        param_name, param_value = key.rsplit('_', 1)
        param_value = float(param_value)

        if param_name not in best_params_per_file:
            best_params_per_file[param_name] = (param_value, data['best_f1_score'])
        elif data['best_f1_score'] > best_params_per_file[param_name][1]:
            best_params_per_file[param_name] = (param_value, data['best_f1_score'])

    return {param: value for param, (value, _) in best_params_per_file.items()}

def load_default_parameters(file_path):
    with open(file_path, 'r') as json_file:
        params = json.load(json_file)
    if 'maf' in params:
        del params['maf']
    return params

def generate_test_values(default_value, num_values=40, increment=0.025):
    increments = np.linspace(-num_values // 2, num_values // 2, num_values) * increment
    return default_value * (1 + increments)

def create_test_object(default_params, param_to_test, test_values):
    param_mapping = {
        'cor': 'cor_interval',
        'ii': 'ii_interval',
        'pp': 'pp_interval',
        'np': 'np_interval',
        'dp': 'dp_interval'
    }

    test_object = {}
    for param, default_value in default_params.items():
        mapped_param = param_mapping.get(param)
        if mapped_param is None:
            continue

        if param == param_to_test:
            test_object[mapped_param] = test_values
        else:
            test_object[mapped_param] = [default_value]

    return test_object


def determine_gradient_value(df, param):
    results = []
    grouped = df.groupby('MAF')
    for maf, group in grouped:
        param_values = group['Parameter Value'].values
        f1_scores = group['F1 Score'].values

        # Merging duplicate parameter values
        param_f1_map = defaultdict(list)
        for param_value, f1_score in zip(param_values, f1_scores):
            param_f1_map[param_value].append(f1_score)

        param_values_new = np.array(list(param_f1_map.keys()))
        f1_scores_new = np.array([np.mean(f1_scores) for f1_scores in param_f1_map.values()])

        if len(param_values_new) < 2 or len(f1_scores_new) < 2:
            continue
        param_values_range = np.max(param_values_new) - np.min(param_values_new)
        f1_scores_range = np.max(f1_scores_new) - np.min(f1_scores_new)
        if param_values_range == 0 or f1_scores_range == 0:
            continue

        f1_scores_first_value = f1_scores_new[0]
        f1_scores_normalized = f1_scores_new - f1_scores_first_value
        spline = UnivariateSpline(param_values_new, f1_scores_normalized, s=None)
        derivative = spline.derivative()
        param_dense_normalized = np.linspace(min(param_values_new), max(param_values_new), 450)
        f1_dense_normalized = spline(param_dense_normalized)

        print(param_values_new)
        print(f1_scores_normalized)
        print('param_dense_normalized:')
        print(param_dense_normalized)
        print('f1_dense_normalized:')
        print(f1_dense_normalized)
        sys.exit()

        derivative_values_normalized = derivative(param_dense_normalized)
        print('derivative')
        print(derivative_values_normalized)
        sys.exit()
        target_slope = np.tan(np.radians(20))
        valid_param_values = []
        for idx in range(len(derivative_values_normalized)):
            if abs(derivative_values_normalized[idx] - target_slope) < 0.02 and f1_dense_normalized[idx] > f1_dense_normalized[0]:
                valid_param_value = param_values_new.min() + param_dense_normalized[idx] * (param_values_new.max() - np.min(param_values_new))
                valid_param_values.append(valid_param_value)
        min_slope_angle = np.degrees(np.arctan(np.min(derivative_values_normalized)))
        first_f1_score = f1_dense_normalized[0]
        last_f1_score = f1_dense_normalized[-1]
        if valid_param_values:
            param_value_to_return = valid_param_values[-1]
        else:
            param_value_to_return = param_values_new[0]
        results.append({
            'maf': maf,
            'param': param,
            'param_value_to_return': param_value_to_return,
            'first_f1_score': first_f1_score,
            'last_f1_score': last_f1_score,
            'lowest_slope_angle': min_slope_angle,
            'valid_param_values': valid_param_values
        })
    return results

# Usage example with a dummy DataFrame:
import pandas as pd

def process_total_parameter_results(total_parameter_results):
    result_dict = {}
    for param, maf_data in total_parameter_results.items():
        for maf, batch_data in maf_data.items():
            df_rows = []
            for batch_id, (param_values, f1_scores) in batch_data.items():
                for param_value, f1_score in zip(param_values, f1_scores):
                    df_rows.append({'MAF': maf, 'Batch ID': batch_id, 'F1 Score': f1_score, 'Parameter Value': param_value})
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
            print(os.path.join(output_training_folder, "maf_{}".format(maf)))
            for folder in os.listdir(os.path.join(output_training_folder, "maf_{}".format(maf))):
                if folder.startswith(param):
                    batch_id = int(folder.split('_')[-2][5:])
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

                    # Sort the lists
                    total_parameter_results[param][maf][batch_id][0].sort()
                    total_parameter_results[param][maf][batch_id][1].sort()

    # Calculate averages and create new object with unique values
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

def load_top_hit(file_path, param_to_fetch):
    df = pd.read_csv(file_path)

    if df.isnull().values.any() or df.empty:
        return None, None

    top_hit = df.iloc[0]
    f1_score = top_hit['F1 Score']
    parameters = top_hit['Parameters']

    if isinstance(parameters, bytes):
        parameters = parameters.decode('utf-8')
    elif not isinstance(parameters, str):
        parameters = str(parameters)

    param_pattern = r'{}_([0-9.]+)'.format(param_to_fetch)
    match = re.search(param_pattern, parameters)
    if match:
        param_value = float(match.group(1))
    else:
        raise ValueError(f"Parameter {param_to_fetch} not found in the parameters string.")

    return f1_score, param_value

output_training_folder = 'nanomgt_training_output'
os.makedirs(output_training_folder, exist_ok=True)
param_list = ['np', 'cor', 'pp', 'dp', 'ii']

"""
for maf in maf_interval:
    os.makedirs(output_training_folder + '/maf_' + str(maf), exist_ok=True)
    for folder in folders:
        if folder.startswith('depth220_SRR27755678'):
            batch_id = int(folder.split('_')[-2][5:])
            if batch_id >= maf:
                input_file_path = os.path.join(alignment_results_path, folder)
                alignment_folder = '/home/people/malhal/test/training_test/{}'.format(folder)
                new_output_folder = output_training_folder + '/' + 'maf_' + str(maf) + '/' + folder
                os.makedirs(new_output_folder, exist_ok=True)
                run_jobs_in_parallel(cpus, new_output_folder, alignment_folder, maf / 100,
                                     parameters_interval_search, maps_path, simulated_batches_csv_path)
                #train_parameters(maf / 100, alignment_folder, 3, 0.4, new_output_folder, maps_path, simulated_batches_csv_path,
                #    0.1, 5, 15, 0.44, 5, 0.15)

all_best_params = defaultdict(list)

for maf in maf_interval:
    print(f"maf_{maf}")
    average_best_params = {}
    for folder in folders:
        if folder.startswith('depth220_SRR27755678'):
            batch_id = int(folder.split('_')[-2][5:])
            if batch_id >= maf:
                new_output_folder = output_training_folder + '/' + 'maf_' + str(maf) + '/' + folder
                results_filename = new_output_folder + "/all_results.csv"
                best_params = calculate_best_parameters(results_filename)
                for param, value in best_params.items():
                    all_best_params[param[1:]].append(value)
    for param, values in all_best_params.items():
        average_value = sum(values) / len(values)
        average_best_params[param] = average_value
        print(f"Average of best {param}: {average_value:.4f}")
    output_file_path = os.path.join(output_training_folder, "{}_average_best_params.json".format('maf_' + str(maf)))
    with open(output_file_path, 'w') as json_file:
        json.dump(average_best_params, json_file, indent=4)

    print(f"Averages saved to {output_file_path}")


# Test individual parameters
for maf in maf_interval:
    output_file_path = os.path.join(output_training_folder, "{}_average_best_params.json".format('maf_' + str(maf)))
    default_params = load_default_parameters(output_file_path)
    for param, default_value in default_params.items():
        if param != 'af':  # Remove this later when maf is not saved
            test_values = generate_test_values(default_value)
            test_object = create_test_object(default_params, param, test_values)

            for folder in folders:
                if folder.startswith('depth'):
                    batch_id = int(folder.split('_')[-2][5:])
                    if batch_id >= maf:
                        new_output_folder = output_training_folder + '/' + 'maf_' + str(maf) + '/' + param + '_' + folder
                        input_file_path = os.path.join(alignment_results_path, folder)
                        alignment_folder = '/home/people/malhal/test/training_test/{}'.format(folder)
                        os.makedirs(new_output_folder, exist_ok=True)
                        if cpus > 41:  # Only training 20 values
                            cpus = 40
                        run_jobs_in_parallel(cpus, new_output_folder, alignment_folder, maf / 100,
                                             test_object, maps_path, simulated_batches_csv_path)
"""
# Eval each parameter value
total_parameter_results = load_results(param_list, maf_interval, output_training_folder)

print (total_parameter_results['np'][2][10])

for maf in maf_interval:
    for param in param_list:
        output_file_csv = os.path.join(output_training_folder, '{}_{}_results.csv'.format(param, maf))
        with open(output_file_csv, mode='w', newline='') as file:
            writer = csv.writer(file)
            writer.writerow(['MAF', 'Batch ID', 'F1 Score', 'Parameter Value'])
            for batch_id, (param_values, f1_scores) in total_parameter_results[param][maf].items():
                for param_value, f1_score in zip(param_values, f1_scores):
                    writer.writerow([maf, batch_id, f1_score, param_value])

processed_results = process_total_parameter_results(total_parameter_results)
print(processed_results)