import os
import pandas as pd

# Define the folder containing the files
folder_path = './species_lists'

# List of files to exclude
exclude_files = ['all_ont_samples.txt']

# Initialize a list to store the data
data = []

# Loop through each file in the folder
for file_name in os.listdir(folder_path):
    if file_name not in exclude_files:
        file_path = os.path.join(folder_path, file_name)
        
        # Extract the basecalling model and species from the file name
        parts = file_name.split('_')
        model = parts[0]
        species = parts[1]

        # Read the file and store the data
        with open(file_path, 'r') as file:
            for line in file:
                sample_id = line.strip()
                data.append([sample_id, species, model])

# Convert the data to a pandas DataFrame
df = pd.DataFrame(data, columns=['Sample_ID', 'Species', 'Basecalling_Model'])

# Save the DataFrame to a CSV file
output_csv_path = os.path.join(folder_path, 'samples_data.csv')
df.to_csv(output_csv_path, index=False)

print(f"Data has been saved to {output_csv_path}")

