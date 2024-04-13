import random

# List of sequencing IDs

seed = 102
random.seed(seed)  # You can use any integer value as the seed


# Escherichia coli
ecoli_sequencing_ids = ['SRR25689478', 'SRR26899125']

# Staphylococcus aureus
staph_aureus_sequencing_ids = ['SRR28370694', 'ERR8958843', 'SRR27167517']

# Campylobacter jejuni
campylobacter_jejuni_sequencing_ids = ['SRR27638397', 'SRR27710526']

# Salmonella enterica
salmonella_enterica_sequencing_ids = ['SRR28399428', 'SRR27136088', 'SRR27755678']

# Number of batches
num_batches = 10


#with open('simulated_batches/salmonella_enterica_batches.csv', 'w') as f:
with open('simulated_batches/{}_ecoli_batches.csv'.format(seed), 'w') as f:
# Print CSV header
    print("Batch, ID1, Percentage1, ID2, Percentage2", file=f)

# Generate batches
    for i in range(num_batches):
        # Randomly select 2 sequencing IDs
        selected_ids = random.sample(ecoli_sequencing_ids, 2)

        # Calculate percentages
        percent_1 = 99 - i
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

        # Prepare IDs and percentages for printing
        id1, id2 = batch.keys()
        percent1, percent2 = batch.values()

        # Print in CSV format
        print(f"{i + 1}, {id1}, {percent1}%, {id2}, {percent2}%", file=f)

with open('simulated_batches/{}_campylobacter_jejuni_batches.csv'.format(seed), 'w') as f:
    # Print CSV header
    print("Batch, ID1, Percentage1, ID2, Percentage2", file=f)

    # Generate batches
    for i in range(num_batches):
        # Randomly select 2 sequencing IDs
        selected_ids = random.sample(campylobacter_jejuni_sequencing_ids, 2)

        # Calculate percentages
        percent_1 = 99 - i
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

        # Prepare IDs and percentages for printing
        id1, id2 = batch.keys()
        percent1, percent2 = batch.values()

        # Print in CSV format
        print(f"{i + 1}, {id1}, {percent1}%, {id2}, {percent2}%", file=f)


with open('simulated_batches/{}_staph_aureus_batches.csv'.format(seed), 'w') as f:
    # Print CSV header
    print("Batch, ID1, Percentage1, ID2, Percentage2, ID3, Percentage3", file=f)

    # Generate batches
    for i in range(num_batches):
        # Randomly select 2 sequencing IDs
        selected_ids = random.sample(salmonella_enterica_sequencing_ids, 3)

        # Calculate percentages
        percent_1 = 98 - i*2
        percent_2 = i + 1
        percent_3 = i + 1

        # Randomly assign percentages to the selected IDs
        if random.choice([True, False]):
            batch = {
                selected_ids[0]: percent_1,
                selected_ids[1]: percent_2,
                selected_ids[2]: percent_3
            }
        else:
            batch = {
                selected_ids[0]: percent_2,
                selected_ids[1]: percent_1,
                selected_ids[2]: percent_3
            }

        # Prepare IDs and percentages for printing
        id1, id2, id3 = batch.keys()
        percent1, percent2, percent3 = batch.values()

        # Print in CSV format
        print(f"{i + 1}, {id1}, {percent1}%, {id2}, {percent2}%, {id3}, {percent3}%", file=f)

with open('simulated_batches/{}_salmonella_enterica_batches.csv'.format(seed), 'w') as f:
    # Print CSV header
    print("Batch, ID1, Percentage1, ID2, Percentage2, ID3, Percentage3", file=f)

    # Generate batches
    for i in range(num_batches):
        # Randomly select 2 sequencing IDs
        selected_ids = random.sample(salmonella_enterica_sequencing_ids, 3)

        # Calculate percentages
        percent_1 = 98 - i*2
        percent_2 = i + 1
        percent_3 = i + 1

        # Randomly assign percentages to the selected IDs
        if random.choice([True, False]):
            batch = {
                selected_ids[0]: percent_1,
                selected_ids[1]: percent_2,
                selected_ids[2]: percent_3
            }
        else:
            batch = {
                selected_ids[0]: percent_2,
                selected_ids[1]: percent_1,
                selected_ids[2]: percent_3
            }

        # Prepare IDs and percentages for printing
        id1, id2, id3 = batch.keys()
        percent1, percent2, percent3 = batch.values()

        # Print in CSV format
        print(f"{i + 1}, {id1}, {percent1}%, {id2}, {percent2}%, {id3}, {percent3}%", file=f)