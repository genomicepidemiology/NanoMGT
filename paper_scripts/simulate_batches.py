import os

# Define the directory for output
output_dir = 'simulated_batches'
os.makedirs(output_dir, exist_ok=True)

# Define species and their sequencing IDs
species_sequencing_ids = {
    'ecoli': ['SRR25689478', 'SRR26036455', 'ERR12533301'],
    'staph_aureus': ['SRR28370694', 'ERR8958843'],
    'campylobacter_jejuni': ['SRR27638397', 'SRR27710526'],
    'salmonella_enterica': ['SRR28399428', 'SRR27136088', 'SRR27755678'],
    'klebsiella_pneumoniae': ['ERR8958737', 'SRR27348733']
}

# Number of batches
num_batches = 10

def simulate_batches(species_name, ids):
    # For each ID in the list, it should act as the major ID in a file
    for major_index, major_id in enumerate(ids):
        filename = f'{output_dir}/{species_name}_{major_id}_majority_batches.csv'
        with open(filename, 'w') as f:
            if len(ids) == 2:
                header = "Batch, ID1, Percentage1, ID2, Percentage2"
            else:
                header = "Batch, ID1, Percentage1, ID2, Percentage2, ID3, Percentage3"
            print(header, file=f)

            for i in range(num_batches):
                if len(ids) == 2:
                    # Alternate majority between the two IDs
                    if major_index == 0:  # First ID is the majority
                        percentages = [100 - (i + 1), i + 1]
                    else:  # Second ID is the majority
                        percentages = [i + 1, 100 - (i + 1)]
                else:
                    # Handling for three or more IDs
                    major_percentage = 98 - i * 2
                    remaining_percentage = 100 - major_percentage
                    minor_percentages = [remaining_percentage // (len(ids) - 1)] * (len(ids) - 1)
                    minor_percentages[-1] += remaining_percentage % (len(ids) - 1)  # Ensure total is 100%
                    if major_index == 0:
                        percentages = [major_percentage] + minor_percentages
                    else:
                        percentages = minor_percentages[:major_index] + [major_percentage] + minor_percentages[major_index:]

                batch = {id: pct for id, pct in zip(ids, percentages)}
                batch_info = [f"{i + 1}"]
                for id in ids:
                    batch_info.append(f"{id}, {batch[id]}%")

                print(", ".join(batch_info), file=f)

# Simulate batches for each species
for species_name, ids in species_sequencing_ids.items():
    simulate_batches(species_name, ids)
