import os
import csv
from collections import defaultdict

# Path to the main directory
main_dir = "/home/people/malhal/test/new_nanomgt_results"


# Parsing parameters from the parameter string
def parse_parameters(param_str):
    """ Parse the parameter string and return a dictionary of parameter values. """
    try:
        params = param_str.split('_')
        param_dict = {
            'cor': float(params[3]),
            'pp': float(params[5]),
            'np': float(params[7]),
            'dp': float(params[9]),
            'iteration_increase': float(params[-1])
        }
        return param_dict
    except (IndexError, ValueError):
        return None


# Dictionary to store average calculations
average_scores = defaultdict(list)

# Traverse through each subfolder
for i in range(1, 6):
    subfolder = f"np_grid_parameter_output_{i}"
    subfolder_path = os.path.join(main_dir, subfolder)

    # Check each results directory in the subfolder
    for results_dir in os.listdir(subfolder_path):
        results_path = os.path.join(subfolder_path, results_dir)

        # Check for the CSV file
        csv_file = os.path.join(results_path, 'all_results.csv')
        if os.path.isfile(csv_file):
            # Extract the batch number from folder name
            batch_id = results_dir.split('_')[-2][5:]

            with open(csv_file, mode='r', newline='') as file:
                reader = csv.DictReader(file)
                for row in reader:
                    params = parse_parameters(row['Parameters'])
                    if params is not None:
                        np_value = params['np']
                        f1_score = float(row['F1 Score'])
                        # Use tuple (folder number, batch_id, np_value) as key
                        average_scores[(i, batch_id, np_value)].append(f1_score)

# Calculate average F1 scores and prepare to save them
output_data = []
for key, scores in average_scores.items():
    avg_score = sum(scores) / len(scores)
    folder_number, batch_id, np_value = key
    output_data.append({
        'folder_number': folder_number,
        'batch_id': batch_id,
        'np': np_value,
        'average_f1_score': avg_score
    })

# Write to CSV
output_filename = 'average_f1_scores_for_np.csv'  # Change to your desired path
with open(output_filename, mode='w', newline='') as file:
    fieldnames = ['folder_number', 'batch_id', 'np', 'average_f1_score']
    writer = csv.DictWriter(file, fieldnames=fieldnames)
    writer.writeheader()
    for data in output_data:
        writer.writerow(data)

print("CSV file has been created with average F1 scores.")
