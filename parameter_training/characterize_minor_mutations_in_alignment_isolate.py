import os
import json
import csv
from collections import defaultdict

alignment_path = "/home/people/malhal/test/characterize_sup/"
isolate_file_path = "/home/people/malhal/test/characterize_sup/sup_isolates.json"
# Load isolate file data
with open(isolate_file_path, 'r') as f:
    isolate_data = json.load(f)

# Initialize data structures to hold results
results = []

# Process each isolate and corresponding alignment data
for organism, isolates in isolate_data.items():
    for isolate in isolates:
        isolate_path = os.path.join(alignment_path, isolate)
        minor_mutations_file = os.path.join(isolate_path, 'minor_mutations.csv')

        if not os.path.exists(minor_mutations_file):
            print(f"Minor mutations file not found for isolate {isolate}")
            continue

        # Initialize counters
        total_mutations = 0
        proximity_mutations = 0
        novel_mutations = 0
        co_occurring_mutations = 0
        proximity_density_list = []

        # Read minor mutations file
        mutations = []
        with open(minor_mutations_file, 'r') as csvfile:
            reader = csv.DictReader(csvfile)
            for row in reader:
                mutations.append(row)

                total_mutations += 1

                if row['MutationComment'] == 'Novel mutation':
                    novel_mutations += 1

                if row['CoOccurrence'] == 'Yes':
                    co_occurring_mutations += 1

        # Calculate proximity mutations and proximity density
        positions = [int(mutation['Position']) for mutation in mutations]
        for i, pos in enumerate(positions):
            nearby_mutations = sum(1 for p in positions if abs(p - pos) <= 15)
            proximity_density_list.append(nearby_mutations - 1)  # Subtract 1 to exclude the mutation itself

            if any(abs(pos - p) <= 5 for p in positions if pos != p):
                proximity_mutations += 1

        average_proximity_density = sum(proximity_density_list) / len(proximity_density_list) if proximity_density_list else 0

        # Calculate percentages
        percent_co_occurring = (co_occurring_mutations / total_mutations) * 100 if total_mutations else 0
        percent_novel = (novel_mutations / total_mutations) * 100 if total_mutations else 0
        percent_proximity = (proximity_mutations / total_mutations) * 100 if total_mutations else 0

        # Store results
        results.append({
            'Organism': organism,
            'Isolate': isolate,
            'Total_Mutations': total_mutations,
            'Proximity_Mutations': proximity_mutations,
            'Novel_Mutations': novel_mutations,
            'Average_Proximity_Density': average_proximity_density,
            'Co_Occurring_Mutations': co_occurring_mutations,
            'Percent_Co_Occurring': percent_co_occurring,
            'Percent_Novel': percent_novel,
            'Percent_Proximity': percent_proximity
        })

# Output results in CSV format
output_csv_path = "characterization.csv"
csv_columns = ['Organism', 'Isolate', 'Total_Mutations', 'Proximity_Mutations', 'Novel_Mutations', 'Average_Proximity_Density', 'Co_Occurring_Mutations', 'Percent_Co_Occurring', 'Percent_Novel', 'Percent_Proximity']

with open(output_csv_path, 'w', newline='') as csvfile:
    writer = csv.DictWriter(csvfile, fieldnames=csv_columns)
    writer.writeheader()
    for data in results:
        writer.writerow(data)

print(f"CSV file created at {output_csv_path}")
