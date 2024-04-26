import pandas as pd
import re


def load_data_and_calculate_averages(file_path):
    # Load the CSV data
    data = pd.read_csv(file_path)

    # Prepare a regex to extract the parameters
    param_regex = re.compile(r"cor_(\d+\.\d+)_pp_(\d+\.\d+)_np_(\d+)_dp_(\d+\.\d+)_iteration_increase_(\d+\.\d+)")

    # Lists to store the parameters
    cor_values = []
    pp_values = []
    np_values = []
    dp_values = []
    iteration_values = []

    # Iterate over each row in the DataFrame
    for params in data['Parameters']:
        # Ensure params is a string
        if isinstance(params, str):
            match = param_regex.search(params)
            if match:
                cor_values.append(float(match.group(1)))
                pp_values.append(float(match.group(2)))
                np_values.append(int(match.group(3)))
                dp_values.append(float(match.group(4)))
                iteration_values.append(float(match.group(5)))

    # Calculate and return averages safely avoiding division by zero
    avg_cor = sum(cor_values) / len(cor_values) if cor_values else 0
    avg_pp = sum(pp_values) / len(pp_values) if pp_values else 0
    avg_np = sum(np_values) / len(np_values) if np_values else 0
    avg_dp = sum(dp_values) / len(dp_values) if dp_values else 0
    avg_iteration = sum(iteration_values) / len(iteration_values) if iteration_values else 0

    # Output the averages
    print(f"Average cor: {avg_cor}")
    print(f"Average pp: {avg_pp}")
    print(f"Average np: {avg_np}")
    print(f"Average dp: {avg_dp}")
    print(f"Average iteration increase: {avg_iteration}")


def main():
    file_path = 'grid_search_5_all.csv'  # Modify this to the path of your CSV file
    load_data_and_calculate_averages(file_path)


if __name__ == "__main__":
    main()
