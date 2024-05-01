import os
import pandas as pd

# Define the base directory
base_dir = '/home/people/malhal/test/new_nanomgt_results/'

# Define directory names
dir_names = [f'np_search_output_{i}' for i in range(1, 6)] # correct to 6

# Loop over each directory
for idx, dir_name in enumerate(dir_names):
    root_dir = os.path.join(base_dir, dir_name)
    dataframes = []

    # Walk through the directories within each parameter_output folder
    for subdir, dirs, files in os.walk(root_dir):
        for file in files:
            if file == 'top_result.csv':
                # Construct full file path
                file_path = os.path.join(subdir, file)
                # Read the csv file
                df = pd.read_csv(file_path)
                # Extract folder name from subdir and add as a new column
                df['experiment'] = os.path.basename(subdir)
                # Append dataframe to list
                dataframes.append(df)

    # Concatenate all dataframes into one
    all_data = pd.concat(dataframes, ignore_index=True)

    # Output to csv file in the current working directory
    output_filename = f'np_search_{idx + 1}_all.csv'
    all_data.to_csv(output_filename, index=False)

    print(f"CSV file {output_filename} has been created successfully.")