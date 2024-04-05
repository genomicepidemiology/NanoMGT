import random

# List of sequencing IDs
# Escherichia coli
ecoli_sequencing_ids = ['SRR28399430', 'SRR26899125', 'SRR26899129']

# Staphylococcus aureus
staph_aureus_sequencing_ids = ['SRR28370694', 'ERR8958843', 'SRR27167517']

# Acinetobacter baumannii
acinetobacter_baumannii_sequencing_ids = ['SRR26643493', 'SRR27710531', 'SRR23387317']

# Campylobacter jejuni
campylobacter_jejuni_sequencing_ids = ['SRR27638397', 'SRR26899121', 'SRR27710526']

# Salmonella enterica
salmonella_enterica_sequencing_ids = ['SRR28399428', 'SRR27136088', 'SRR27755678']

# Number of batches
num_batches = 10

# Print CSV header
print("Batch, ID1, Percentage1, ID2, Percentage2")

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
    print(f"{i + 1}, {id1}, {percent1}%, {id2}, {percent2}%")
