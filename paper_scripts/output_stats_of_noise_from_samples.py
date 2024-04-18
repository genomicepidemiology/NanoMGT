import os
import pandas as pd

def count_proximity_mutations(df):
    proximity_count = 0
    mutations_in_proximity = set()  # Set to store mutations within proximity
    df.sort_values(by=['Gene', 'Position'], inplace=True)
    grouped = df.groupby('Gene')
    for _, group in grouped:
        positions = group['Position'].values
        # Construct mutation identifiers as 'Position_MajorityBase_to_MutationBase'
        mutations = [f"{pos}_{major}to{mut}" for pos, major, mut in zip(group['Position'], group['MajorityBase'], group['MutationBase'])]
        for i in range(len(positions) - 1):
            if positions[i+1] - positions[i] <= 5:
                proximity_count += 1
                mutations_in_proximity.update([mutations[i], mutations[i+1]])
    return proximity_count, mutations_in_proximity

def calculate_proximity_density(mutations):
    mutation_positions = {int(mut.split('_')[0]): mut for mut in mutations}

    def density(mutation):
        position = int(mutation.split('_')[0])
        count = 0
        for other_position in mutation_positions.keys():
            if abs(other_position - position) <= 15 and other_position != position:
                count += 1
        return count

    total_density = sum(density(mut) for mut in mutations)
    return total_density, len(mutations)

def analyze_mutations(folder_path):
    results = {}
    for root, dirs, files in os.walk(folder_path):
        for name in dirs:
            subfolder_path = os.path.join(root, name)
            csv_file_path = os.path.join(subfolder_path, 'minor_mutations.csv')
            if os.path.exists(csv_file_path):
                df = pd.read_csv(csv_file_path)
                total_mutations = len(df)
                proximity_mutations, proxi_mutations_set = count_proximity_mutations(df)
                novel_mutations = df[df['MutationComment'].str.contains('Novel mutation', na=False)].shape[0]
                if proxi_mutations_set:
                    total_proximity_density, _ = calculate_proximity_density(proxi_mutations_set)
                    average_proximity_density = total_proximity_density / proximity_mutations if proximity_mutations > 0 else 0
                else:
                    average_proximity_density = 0
                results[name] = {
                    'Total Mutations': total_mutations,
                    'Proximity Mutations': proximity_mutations,
                    'Novel Mutations': novel_mutations,
                    'Average Proximity Density': average_proximity_density
                }
    return results

# Path to the main directory containing the subfolders
folder_path = '/home/people/malhal/data/new_nanomgt/majority_variants'
mutation_stats = analyze_mutations(folder_path)
for sample, stats in mutation_stats.items():
    print(f"Sample: {sample}")
    for key, value in stats.items():
        print(f"{key}: {value}")
    print()
