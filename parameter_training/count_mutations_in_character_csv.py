import csv

clean_ids = [
    "SRR27638397", "SRR27710526", "SRR28370668", "SRR25689478", "SRR26036455", "SRR28789754", "SRR26036458",
    "ERR8958810", "ERR8958737", "SRR27348733", "SRR29213739", "SRR24833081", "SRR28370694", "SRR28370638",
    "SRR25865495", "SRR27755684", "SRR27755678", "SRR26899146", "SRR27755667", "SRR27755674", "SRR26899103",
    "SRR25999202", "SRR24833086", "ERR8958866"
]

contaminated_ids = [
    "SRR26899118", "SRR26353490", "SRR27710532", "SRR24837710", "SRR24837712", "SRR24834173", "SRR28789463",
    "SRR28789469", "SRR28800580", "ERR8958848", "SRR25890190", "ERR8958843", "SRR27136090", "SRR28399428", "SRR27136088"
]

def process_data(file_path, clean_ids, contaminated_ids):
    clean_data = {"files": 0, "mutations": 0, "proximity": 0, "co_occuring": 0, "novel": 0}
    contaminated_data = {"files": 0, "mutations": 0, "proximity": 0, "co_occuring": 0, "novel": 0}

    with open(file_path, 'r') as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            isolate = row["Isolate"]
            total_mutations = int(row["Total_Mutations"])
            proximity_mutations = int(row["Proximity_Mutations"])
            co_occuring_mutations = int(row["Co_Occurring_Mutations"])
            novel_mutations = int(row["Novel_Mutations"])

            if isolate in clean_ids:
                clean_data["files"] += 1
                clean_data["mutations"] += total_mutations
                clean_data["proximity"] += proximity_mutations
                clean_data["co_occuring"] += co_occuring_mutations
                clean_data["novel"] += novel_mutations
            elif isolate in contaminated_ids:
                contaminated_data["files"] += 1
                contaminated_data["mutations"] += total_mutations
                contaminated_data["proximity"] += proximity_mutations
                contaminated_data["co_occuring"] += co_occuring_mutations
                contaminated_data["novel"] += novel_mutations

    combined_data = {
        "files": clean_data["files"] + contaminated_data["files"],
        "mutations": clean_data["mutations"] + contaminated_data["mutations"],
        "proximity": clean_data["proximity"] + contaminated_data["proximity"],
        "co_occuring": clean_data["co_occuring"] + contaminated_data["co_occuring"],
        "novel": clean_data["novel"] + contaminated_data["novel"]
    }

    return clean_data, contaminated_data, combined_data

def print_results(data, label):
    print(f"{label} Data Set:")
    print(f"  Number of .txt files: {data['files']}")
    print(f"  Total mutations: {data['mutations']}")
    print(f"  Proximity mutations: {data['proximity']}")
    print(f"  Co-occurring mutations: {data['co_occuring']}")
    print(f"  Novel mutations: {data['novel']}")
    print()

# Process the data
file_path = 'character.csv'
clean_data, contaminated_data, combined_data = process_data(file_path, clean_ids, contaminated_ids)

# Print the results
print_results(clean_data, "Clean")
print_results(contaminated_data, "Contaminated")
print_results(combined_data, "Combined")
