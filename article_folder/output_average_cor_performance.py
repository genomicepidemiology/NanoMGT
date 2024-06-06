import os
import pandas as pd
import re

def read_and_aggregate_results(base_directory):
    # Define regex to extract 'maf', 'cor', and batch number from the folder name and Parameters
    param_regex = re.compile(r"maf_(\d+\.\d+)_.*?cor_(\d+\.\d+)_")
    batch_regex = re.compile(r"batch(\d+)_merged$")

    # DataFrame to hold all data
    all_data = pd.DataFrame()

    # List all 'cor_search_output_*' directories, assuming directory naming follows a similar pattern for 'cor'
    search_dirs = [d for d in os.listdir(base_directory) if d.startswith('cor_search_output')]

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

                        # Extract 'maf' and 'cor' values and add them as columns
                        params = df['Parameters'].apply(lambda x: param_regex.search(x))
                        df['maf'] = params.apply(lambda x: float(x.group(1)) if x else None)
                        df['cor'] = params.apply(lambda x: float(x.group(2)) if x else None)
                        df['batch'] = batch_number

                        # Append to the all_data DataFrame
                        all_data = pd.concat([all_data, df], ignore_index=True)

    # Group by 'batch', 'maf', and 'cor' and calculate the mean of F1 Score
    result = all_data.groupby(['batch', 'maf', 'cor'])['F1 Score'].mean().reset_index()
    result.columns = ['batch', 'maf', 'cor', 'average_f1']

    return result

def main():
    base_directory = "/home/people/malhal/test/new_nanomgt_results"
    results = read_and_aggregate_results(base_directory)
    # Output the results to a CSV file
    results.to_csv("cor_f1_scores.csv", index=False)
    print("Average F1 scores computed and saved.")

if __name__ == "__main__":
    main()
