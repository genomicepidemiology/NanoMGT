import os
import json
import random

# Define the directory for output
output_dir = 'simulated_batches_sup'
os.makedirs(output_dir, exist_ok=True)

# Load genome sizes
with open('genome_sizes.json', 'r') as f:
    genome_sizes = json.load(f)

# Load training and validation sets
with open('training_set_sup_isolates.json', 'r') as f:
    training_set = json.load(f)

with open('validation_set_sup_isolates.json', 'r') as f:
    validation_set = json.load(f)

# Number of batches
num_batches = 10
sample_complexity = 5

def simulate_batches(ids, species, dataset_type, sample_complexity):
    all_batches = []
    batch_id_counter = 1

    for seed in range(sample_complexity):  # Simulate sample_complexity different seeds for each major ID
        random.seed(seed)
        random.shuffle(ids)

        for i in range(num_batches):
            if len(ids) == 2:
                # Alternate majority between the two IDs
                major_percentage = 98 - i * 2
                percentages = [major_percentage, 100 - major_percentage]
            else:
                # Handling for three or more IDs
                major_percentage = 98 - i * 2
                remaining_percentage = 100 - major_percentage
                minor_percentages = [remaining_percentage // (len(ids) - 1)] * (len(ids) - 1)
                minor_percentages[-1] += remaining_percentage % (len(ids) - 1)  # Ensure total is 100%

                percentages = [major_percentage] + minor_percentages

            batch = {id: pct for id, pct in zip(ids, percentages)}
            majority_id = max(batch, key=batch.get)
            minority_ids = {k: v for k, v in batch.items() if k != majority_id}

            batch_info = {
                "batch_id": batch_id_counter,
                "minority_abundance": {
                    "Majority": {majority_id: batch[majority_id]},
                    "Minority": minority_ids
                }
            }
            all_batches.append(batch_info)
            batch_id_counter += 1

    filename = f'{output_dir}/{species}_{dataset_type}.json'
    with open(filename, 'w') as f:
        json.dump(all_batches, f, indent=4)

def process_data(data, dataset_type):
    for species, ids in data.items():
        if len(ids) < 2:
            continue  # Skip if there are fewer than 2 IDs
        sample_complexity_actual = min(len(ids), sample_complexity)
        selected_ids = random.sample(ids, sample_complexity_actual)
        species = species.lower().replace(' ', '_')
        simulate_batches(selected_ids, species, dataset_type, sample_complexity)

# Simulate batches for each species in both training and validation sets
process_data(training_set, 'training')
process_data(validation_set, 'validation')

print("Batch simulation complete.")
