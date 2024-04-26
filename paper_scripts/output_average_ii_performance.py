import os
import pandas as pd
import re

def read_and_aggregate_results(base_directory):
    # Define regex to extract 'mrd', 'iteration_increase', and batch number from the folder name and Parameters
    param_regex = re.compile(r"mrd_(\d+\.\d+)_.*?iteration_increase_(\d+\.\d+)")
    batch_regex = re.compile(r"batch(\d+)_merged$")

    # DataFrame to hold all data
    all_data = pd.DataFrame()

    # List all directories that might contain results
    search_dirs = [d for d in os.listdir(base_directory) if 'ii_search_output' in d]

    # Process each directory
    for search_dir in search_dirs:
        current_dir = os.path.join(base_directory, search_dir)
        for root, dirs, files in os.walk(current_dir):
            for dir_name in dirs:
                if dir_name.endswith('merged'):
                    # Extract batch number
                    batch_match = batch_regex.search(dir_name)
                    if batch_match:
                        batch_number = int(batch_match.group(1))
                    else:
                        continue  # Skip if no batch number is found

                    file_path = os.path.join(root, dir_name, 'all_results.csv')
                    if os.path.exists(file_path):
                        # Read the CSV file
                        df = pd.read_csv(file_path)

                        # Extract 'mrd' and 'ii' values and add them as columns
                        df['mrd'] = df['Parameters'].apply(lambda x: float(param_regex.search(x).group(1)) if param_regex.search(x) else None)
                        df['ii'] = df['Parameters'].apply(lambda x: float(param_regex.search(x).group(2)) if param_regex.search(x) else None)
                        df['batch'] = batch_number

                        # Append to the all_data DataFrame
                        all_data = pd.concat([all_data, df], ignore_index=True)

    # Group by 'batch', 'mrd', and 'ii' and calculate the mean of F1 Score
    result = all_data.groupby(['batch', 'mrd', 'ii'])['F1 Score'].mean().reset_index()
    result.columns = ['batch', 'mrd', 'ii', 'average_f1']

    return result

def main():
    base_directory = "/home/people/malhal/test/new_nanomgt_results"
    results = read_and_aggregate_results(base_directory)
    # Output the results to a CSV file
    results.to_csv("ii_f1_scores.csv", index=False)
    print("Average F1 scores computed and saved.")

if __name__ == "__main__":
    main()
