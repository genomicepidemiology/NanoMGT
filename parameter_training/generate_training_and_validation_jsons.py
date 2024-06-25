import json
import math


def divide_into_sets(data):
    training_set = {}
    validation_set = {}

    for species, ids in data.items():
        num_ids = len(ids)
        if num_ids <= 3:
            training_set[species] = ids
            validation_set[species] = []
        elif num_ids == 4:
            half_point = 2
            training_set[species] = ids[:half_point]
            validation_set[species] = ids[half_point:]
        elif num_ids == 5:
            training_set[species] = ids[:3]
            validation_set[species] = ids[3:]
        else:
            # Ensure a minimum of 3 IDs in each set
            base_size = 3
            remaining_ids = num_ids - (2 * base_size)
            training_size = base_size + math.ceil(0.7 * remaining_ids)
            validation_size = num_ids - training_size

            training_set[species] = ids[:training_size]
            validation_set[species] = ids[training_size:training_size + validation_size]

    return training_set, validation_set


# Sample JSON input (replace this with the path to your actual input file)
input_file_path = "contaminated_isolates.json"  # Replace this with your actual input file path
with open(input_file_path, "r") as json_file:
    data = json.load(json_file)

# Divide the IDs into training and validation sets
training_set, validation_set = divide_into_sets(data)

# Save the training and validation sets to JSON files
output_training_path = "training_set_contaminated_isolates.json"
output_validation_path = "validation_set_contaminated_isolates.json"

with open(output_training_path, "w") as json_file:
    json.dump(training_set, json_file, indent=4)

with open(output_validation_path, "w") as json_file:
    json.dump(validation_set, json_file, indent=4)

print(f"Training set saved to {output_training_path}")
print(f"Validation set saved to {output_validation_path}")
