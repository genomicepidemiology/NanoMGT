import os

# Define the directory for output
output_dir = 'simulated_batches'
os.makedirs(output_dir, exist_ok=True)

# Escherichia coli
ecoli_sequencing_ids = ['SRR25689478', 'SRR26036455']

# Staphylococcus aureus
staph_aureus_sequencing_ids = ['SRR28370694', 'ERR8958843']

# Campylobacter jejuni
campylobacter_jejuni_sequencing_ids = ['SRR27638397', 'SRR27710526']

# Salmonella enterica
salmonella_enterica_sequencing_ids = ['SRR28399428', 'SRR27136088', 'SRR27755678']

# Klebsiella pneumoniae
klebsiella_pneumoniae_sequencing_ids = ['ERR8958737', 'SRR27348733']

# Number of batches
num_batches = 10


def simulate_batches(species_name, ids):
    for major_id in ids:
        filename = f'{output_dir}/{species_name}_{major_id}_majority_batches.csv'
        with open(filename, 'w') as f:
            if len(ids) == 2:
                header = "Batch, ID1, Percentage1, ID2, Percentage2"
            else:
                header = "Batch, ID1, Percentage1, ID2, Percentage2, ID3, Percentage3"
            print(header, file=f)

            for i in range(num_batches):
                if len(ids) == 2:
                    # When there are only two IDs
                    percentages = [100 - (i + 1), i + 1]
                else:
                    # When there are three or more IDs
                    major_percentage = 98 - i * 2
                    remaining_percentage = 100 - major_percentage
                    minor_percentages = [remaining_percentage // (len(ids) - 1)] * (len(ids) - 1)
                    minor_percentages[-1] += remaining_percentage % (len(ids) - 1)  # Ensure total is 100%
                    percentages = [major_percentage] + minor_percentages

                batch = {id: pct for id, pct in zip(ids, percentages)}
                batch_info = [f"{i + 1}"]
                for id in ids:
                    batch_info.append(f"{id}, {batch[id]}%")

                print(", ".join(batch_info), file=f)


simulate_batches('ecoli', ecoli_sequencing_ids)
simulate_batches('staph_aureus', staph_aureus_sequencing_ids)
simulate_batches('campylobacter_jejuni', campylobacter_jejuni_sequencing_ids)
simulate_batches('salmonella_enterica', salmonella_enterica_sequencing_ids)
simulate_batches('klebsiella_pneumoniae', klebsiella_pneumoniae_sequencing_ids)
