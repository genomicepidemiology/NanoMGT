import random

# List of sequencing IDs
sequencing_ids = ['SRR28399430', 'SRR26899125', 'SRR26899129']

# Number of batches
num_batches = 10

# Generate batches
for i in range(num_batches):
    # Randomly select 2 sequencing IDs
    selected_ids = random.sample(sequencing_ids, 2)

    # Calculate percentages
    percent_1 = 100 - i
    percent_2 = i + 1

    # Randomly assign percentages to the selected IDs
    if random.choice([True, False]):
        batch = {
            selected_ids[0]: percent_1,
            selected_ids[1]: percent_2
        }
    else:
        batch = {
            selected_ids[0]: percent_2,
            selected_ids[1]: percent_1
        }

    # Print the batch and its percentages
    print(f"Batch {i + 1}:")
    for id, percent in batch.items():
        print(f"{id}: {percent}%")
    print()  # Newline for better readability
