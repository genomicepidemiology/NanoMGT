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

maf_interval = [4, 5]

cor_interval_search = [0.3, 0.5, 0.7]
dp_interval_search = [0.1, 0.3, 0.5]
np_interval_search = [1.5, 2.5, 3.5]
pp_interval_search = [0.3, 0.5, 0.7]
ii_interval_search = [0.01, 0.2, 0.40]

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
    arguments = argparse.Namespace()
    arguments.maf = maf
    arguments.output = results_folder
    arguments.min_n = min_n
    arguments.cor = cor
    arguments.new_output = new_output_folder
    arguments.ii = ii
    arguments.proxi = proxi
    arguments.dp = dp
    arguments.pp = pp
    arguments.np = np
    arguments.dp_window = dp_window

    consensus_dict = nvc.build_consensus_dict(os.path.join(arguments.output, 'rmlst_alignment.res'),
                                          os.path.join(arguments.output, 'rmlst_alignment.mat'))

    confirmed_mutation_dict = nvc.derive_mutation_positions(consensus_dict, arguments)

    bio_validation_dict = nvc.bio_validation_mutations(consensus_dict, os.path.join(results_folder, 'specie.fsa'))

    confirmed_mutation_dict, co_occurrence_tmp_dict, iteration_count =\
        nvc.co_occurrence_until_convergence(arguments, confirmed_mutation_dict,
                                        consensus_dict, {}, bio_validation_dict)

    format_output(new_output_folder, confirmed_mutation_dict, consensus_dict, bio_validation_dict,
                  co_occurrence_tmp_dict)

    sample = arguments.output.split('/')[-1]

    minor_mutation_expected = benchmark_analysis_result(sample, simulated_batches_csv_path, maps_path)

    minor_mutation_results = convert_mutation_dict_to_object(confirmed_mutation_dict)

    print (len(minor_mutation_expected), len(minor_mutation_results))

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
    for i in range(sample_number-1):
        map_file = f'{maps_path}/major_{top_id}_minor_{minor[i]}.txt'
        map_files.append(map_file)

    mutation_map = load_mutations_from_files(map_files)

    return mutation_map

def format_output(new_output_folder, confirmed_mutation_dict, consensus_dict, bio_validation_dict, co_occurrence_tmp_dict):
    with open(new_output_folder + '/minor_mutations.csv', 'w') as outfile:
        header = 'Gene,Position,MajorityBase,MutationBase,MutationDepth,TotalDepth,GeneLength,MutationComment,CoOccurrence'
        print(header, file=outfile)
        for allele in confirmed_mutation_dict:
            for mutation in zip(confirmed_mutation_dict[allele][0], confirmed_mutation_dict[allele][1]):
                position = mutation[0].split('_')[0]
                mutation_base = mutation[0].split('_')[1]
                mutation_depth = mutation[1]
                majority_base = consensus_dict[allele][1][int(position) - 1]
                total_depth = sum(consensus_dict[allele][0][int(position) - 1])
                biological_existence = nvc.check_single_mutation_existence(bio_validation_dict, allele, mutation[0])
                gene_length = len(consensus_dict[allele][1])
                if mutation[0] in co_occurrence_tmp_dict[allele]:
                    co_occurrence = 'Yes'
                else:
                    co_occurrence = 'No'

                if biological_existence:
                    print('{},{},{},{},{},{},{},{},{}'.format(allele, position, majority_base, mutation_base,
                                                              mutation_depth, total_depth, gene_length,
                                                              'Mutation seen in database', co_occurrence), file=outfile)
                else:
                    print('{},{},{},{},{},{},{},{},{}'.format(allele, position, majority_base, mutation_base,
                                                              mutation_depth, total_depth, gene_length,
                                                              'Novel mutation', co_occurrence), file=outfile)

def convert_mutation_dict_to_object(mutation_dict):
    mutation_object = {}
    for allele in mutation_dict:
        gene = allele.split('_')[0]
        mutation_object[gene] = set()
        for mutation in mutation_dict[allele][0]:
            mutation_object[gene].add(mutation)
    return mutation_object

def calculate_metrics(expected_mutations, actual_mutations):
    sum_precision, sum_recall, sum_f1 = 0, 0, 0
    genes_counted = 0

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
    # Exclude the 'maf' parameter
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
        print (param)
        mapped_param = param_mapping.get(param)
        if mapped_param is None:
            continue

        if param == param_to_test:
            test_object[mapped_param] = test_values
        else:
            test_object[mapped_param] = [default_value]

    return test_object

def process_data(df, param):
    results = []

    # Group by 'MAF' and 'Parameter Value'
    grouped = df.groupby('MAF')

    for maf, group in grouped:
        param_values = group['Parameter Value'].values
        f1_scores = group['F1 Score'].values

        print(f"Processing MAF: {maf}, Parameter: {param}")
        print(param_values, f1_scores)

        if len(param_values) < 2 or len(f1_scores) < 2:
            continue

        # Normalize param_values and f1_scores to range [0, 1]
        param_values_range = np.max(param_values) - np.min(param_values)
        f1_scores_range = np.max(f1_scores) - np.min(f1_scores)

        if param_values_range == 0 or f1_scores_range == 0:
            continue

        param_values_normalized = (param_values - np.min(param_values)) / param_values_range
        f1_scores_normalized = (f1_scores - np.min(f1_scores)) / f1_scores_range

        # Ensure param_values_normalized is in increasing order
        sorted_indices = np.argsort(param_values_normalized)
        param_values_normalized = param_values_normalized[sorted_indices]
        f1_scores_normalized = f1_scores_normalized[sorted_indices]

        # Fit a spline to the normalized data points
        spline = UnivariateSpline(param_values_normalized, f1_scores_normalized, s=None)
        derivative = spline.derivative()

        # Generate dense param values for a finer analysis in the normalized range
        param_dense_normalized = np.linspace(0, 1, 450)
        f1_dense_normalized = spline(param_dense_normalized)
        derivative_values_normalized = derivative(param_dense_normalized)

        # Calculate the target derivative for a 20-degree angle
        target_slope = np.tan(np.radians(20))

        # Collect all param values where the derivative is close to the 20-degree slope
        valid_param_values = []
        for idx in range(len(derivative_values_normalized)):
            if abs(derivative_values_normalized[idx] - target_slope) < 0.02 and f1_dense_normalized[idx] > f1_dense_normalized[0]:
                valid_param_value = param_values.min() + param_dense_normalized[idx] * (param_values.max() - param_values.min())
                valid_param_values.append(valid_param_value)

        # Calculate the lowest gradient slope angle in degrees
        min_slope_angle = np.degrees(np.arctan(np.min(derivative_values_normalized)))

        # Get the first and last F1 score from the normalized data
        first_f1_score = f1_dense_normalized[0]
        last_f1_score = f1_dense_normalized[-1]

        # Determine the parameter value to return
        if valid_param_values:
            param_value_to_return = valid_param_values[-1]  # Return the last valid value
        else:
            param_value_to_return = param_values[0]  # Return the first value if no valid values found

        # Append results
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

def process_directory(directory):
    result_dict = {}

    # List all files in the directory
    for filename in os.listdir(directory):
        if filename.endswith(".csv"):
            # Extract the parameter name and MAF value from the filename
            parts = filename.split('_')
            param = parts[0]
            maf = parts[1]

            # Read the CSV file
            file_path = os.path.join(directory, filename)
            df = pd.read_csv(file_path)

            # Process the data and get the parameter value to return
            processed_results = process_data(df, param)

            # Initialize the maf dictionary if not present
            if maf not in result_dict:
                result_dict[maf] = {}

            # Store the result
            for result in processed_results:
                result_dict[maf][param] = result['param_value_to_return']

    return result_dict

def load_top_hit(file_path, param_to_fetch):
    df = pd.read_csv(file_path)

    # Check if the DataFrame is empty or has incomplete data
    if df.isnull().values.any() or df.empty:
        return None, None

    top_hit = df.iloc[0]  # Assuming the top hit is the first row
    f1_score = top_hit['F1 Score']
    parameters = top_hit['Parameters']

    # Ensure parameters is a string
    if isinstance(parameters, bytes):
        parameters = parameters.decode('utf-8')
    elif not isinstance(parameters, str):
        parameters = str(parameters)

    # Extracting the specific parameter value
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
#Initial grid search
for maf in maf_interval:
    os.makedirs(output_training_folder + '/maf_' + str(maf), exist_ok=True)
    for folder in folders:
        if folder.startswith('depth'):
            batch_id = int(folder.split('_')[-2][5:])
            if batch_id >= maf:
                input_file_path = os.path.join(alignment_results_path, folder)
                alignment_folder = '/home/people/malhal/test/training_test/{}'.format(folder)
                new_output_folder = output_training_folder + '/' + 'maf_' + str(maf) + '/' + folder
                os.makedirs(new_output_folder, exist_ok=True)
                run_jobs_in_parallel(cpus, new_output_folder, alignment_folder, maf / 100, parameters_interval_search, maps_path, simulated_batches_csv_path)

all_best_params = defaultdict(list)

for maf in maf_interval:
    print(f"maf_{maf}")
    average_best_params = {}
    for folder in folders:
        if folder.startswith('depth'):
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

#Test individual parameters
for maf in maf_interval:
    output_file_path = os.path.join(output_training_folder, "{}_average_best_params.json".format('maf_' + str(maf)))
    default_params = load_default_parameters(output_file_path)
    for param, default_value in default_params.items():
        if param != 'af': #Remove this later when maf is not saved
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
                        if cpus > 41: #Only training 20 values
                            cpus = 40
                        run_jobs_in_parallel(cpus, new_output_folder, alignment_folder, maf / 100,
                                             test_object, maps_path, simulated_batches_csv_path)


#Eval each parameter value
"""

total_parameter_results = {}
for param in param_list:
    total_parameter_results[param] = {}
    for maf in maf_interval:
        total_parameter_results[param][maf] = {}
        for folder in os.path.join(output_training_folder, "{}".format('maf_' + str(maf))):
            if folder.startswith(param):
                batch_id = int(folder.split('_')[-2][5:])
                results_file = os.path.join(output_training_folder, "{}".format('maf_' + str(maf)), folder, 'all_results.csv')
                print (results_file)


"""


for maf in maf_interval:
    total_parameter_dict[maf] = {}
    output_file_path = os.path.join(output_training_folder, "{}_average_best_params.json".format('maf_' + str(maf)))
    default_params = load_default_parameters(output_file_path)

    for param, default_value in default_params.items():
        total_parameter_dict[maf][param] = {}
        test_values = generate_test_values(default_value)
        test_object = create_test_object(default_params, param, test_values)

        for folder in folders:
            if folder.startswith('depth'):
                batch_id = int(folder.split('_')[-2][5:])
                if batch_id >= maf:
                    new_output_folder = output_training_folder + '/' + 'maf_' + str(maf) + '/' + param + '_' + folder
                    results_filename = new_output_folder + "/top_result.csv"
                    f1_score, param_value = load_top_hit(results_filename, param)
                    total_parameter_dict[maf][param][batch_id] = [f1_score, param_value]

print ("Done with fine tuning")
for maf in total_parameter_dict:
    for param in total_parameter_dict[maf]:
        if 'maf' not in param:
            output_file_csv = os.path.join(output_training_folder, '{}_{}_results.csv'.format(param[1:], maf))

            # Open the CSV file for writing
            with open(output_file_csv, mode='w', newline='') as file:
                writer = csv.writer(file)
                # Write the header
                writer.writerow(['MAF', 'Batch ID', 'F1 Score', 'Parameter Value'])

                # Write the data
                for batch_id in total_parameter_dict[maf][param]:
                    maf_value = maf
                    batch_id_value = batch_id
                    f1_score = total_parameter_dict[maf][param][batch_id][0]
                    param_value = total_parameter_dict[maf][param][batch_id][1]
                    writer.writerow([maf_value, batch_id_value, f1_score, param_value])

processed_results = process_directory(output_training_folder)
print(processed_results)
"""