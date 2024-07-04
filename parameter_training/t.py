import pandas as pd

def calculate_average_mutations(file_path):
    # Read the CSV file into a pandas DataFrame
    df = pd.read_csv(file_path)

    # Calculate the average number of total mutations
    average_total_mutations = df['Total_Mutations'].mean()

    return average_total_mutations

# Specify the path to the CSV file
file_path = 'character.csv'

# Calculate and print the average number of mutations
average_mutations = calculate_average_mutations(file_path)
print(f"Average number of mutations across all isolates: {average_mutations:.2f}")

