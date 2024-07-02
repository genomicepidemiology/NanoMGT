import os
import json
import itertools
import random

# Define the directory for output
output_dir = 'simulated_batches_clean'
os.makedirs(output_dir, exist_ok=True)

# Load training and validation sets
with open('training_set_clean_isolates.json', 'r') as f:
    training_set = json.load(f)

with open('validation_set_clean_isolates.json', 'r') as f:
    validation_set = json.load(f)

def generate_combinations(isolates, num_batches):
    combinations = []
    if len(isolates) == 2:
        for i in range(1, num_batches + 1):
            combinations.extend([
                {"Majority": {isolates[0]: 100 - i}, "Minority": {isolates[1]: i}},
                {"Majority": {isolates[1]: 100 - i}, "Minority": {isolates[0]: i}}
            ])
    elif len(isolates) >= 3:
        # Handle three isolates
        for i in range(1, num_batches + 1):
            remaining = (i * 2)  # Distribute equally among two minority isolates
            major = 100 - remaining
            for combo in itertools.permutations(isolates, 3):
                combinations.append(
                    {"Majority": {combo[0]: major}, "Minority": {combo[1]: i, combo[2]: i}}
                )
        # Sort by majority value and select 5 combinations for each unique majority value
        combinations.sort(key=lambda x: next(iter(x["Majority"].values())))
        if len(combinations) > 50:
            selected_combinations = []
            unique_major_values = sorted(set(next(iter(x["Majority"].values())) for x in combinations))
            for value in unique_major_values:
                filtered_combos = [c for c in combinations if next(iter(c["Majority"].values())) == value]
                selected_combinations.extend(random.sample(filtered_combos, min(5, len(filtered_combos))))
            combinations = selected_combinations

    return combinations

def simulate_batches(species_data, dataset_type, num_batches):
    for species, isolates in species_data.items():
        filename = f'{output_dir}/{species.lower().replace(" ", "_")}_{dataset_type}.json'
        sampled_combos = generate_combinations(isolates, num_batches)

        results = [{"batch_id": i+1, "minority_abundance": combo} for i, combo in enumerate(sampled_combos)]
        with open(filename, 'w') as f:
            json.dump(results, f, indent=4)

# Specify the number of batches per minority abundance level (1 to 10)
simulate_batches(training_set, 'training', 10)
simulate_batches(validation_set, 'validation', 10)

print("Batch simulation for both training and validation sets is complete.")
