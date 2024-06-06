import os

path = '/home/people/malhal/data/new_nanomgt'
output = '/home/people/malhal/data/new_nanomgt/genomes/'
database = '/home/people/malhal/databases/bacteria_20231113/bacteria_db'
files = os.listdir(path)

# Mapping filenames to their respective organisms and expected genome sizes
file_to_organism = {
    "SRR25689478": "Escherichia coli",
    "SRR26036455": "Escherichia coli",
    "ERR12533301": "Escherichia coli",
    "ERR8958737": "Klebsiella pneumoniae",
    "SRR27348733": "Klebsiella pneumoniae",
    "SRR27638397": "Campylobacter jejuni",
    "SRR27710526": "Campylobacter jejuni",
    "SRR28399428": "Salmonella enterica",
    "SRR27136088": "Salmonella enterica",
    "SRR27755678": "Salmonella enterica",
    "SRR28370694": "Staphylococcus aureus",
    "ERR8958843": "Staphylococcus aureus"
}

genome_sizes = {
    "Escherichia coli": 5000000,
    "Staphylococcus aureus": 3000000,
    "Campylobacter jejuni": 1700000,
    "Salmonella enterica": 5000000,
    "Klebsiella pneumoniae": 5500000
}

def find_best_template_from_spa_file(spa_file):
    """Returns the mapping results from the reference mapping"""
    if os.path.exists(spa_file):
        results = []
        with open(spa_file, 'r') as f:
            data = f.read().split("\n")
        data = data[:-1]  # Last line is empty
        for item in data:
            item = item.split("\t")
            if item[0][0] != "#":
                results.append((item[0], int(item[1]), float(item[2])))
        results.sort(key=lambda x: x[2], reverse=True)  # Sort by score descending
        return results
    return []

def get_genome_size(reference_header_text, database):
    """Returns the genome size of the given reference from the database"""
    genome_size = 0
    with open(database, 'r') as f:
        read = False
        for line in f:
            if line.startswith('>'):
                if read:
                    break
                if reference_header_text in line:
                    read = True
            elif read:
                genome_size += len(line.strip())
    return genome_size

def check_template_size(spa_file, database, expected_size):
    """Check if any template from spa_file meets the size criteria"""
    results = find_best_template_from_spa_file(spa_file)
    for reference_header_text, template_number, template_score in results:
        genome_size = get_genome_size(reference_header_text, database)
        if genome_size >= 0.7 * expected_size:
            return template_number, template_score, reference_header_text
    return None, None, None

for file in files:
    if file.endswith('.fastq'):
        base_name = file.split('.')[0]
        organism = file_to_organism.get(base_name, None)
        expected_size = genome_sizes.get(organism, None)
        if expected_size is None:
            print(f"Unknown organism for file {file}")
            continue

        os.system(f'mkdir {output}/{file[:-6]}')
        cmd = f'kma -i {path}/{file} -o {output}/{file[:-6]}/read_mapping -t_db {database} -t 8 -mem_mode -Sparse -ss c'
        os.system(cmd)

        template_number, template_score, reference_header_text = check_template_size(
            f'{output}/{file[:-6]}/read_mapping.spa',
            '/home/people/malhal/databases/bacteria_20231113/bacteria',
            expected_size
        )

        if template_number is not None:
            print(template_number, template_score, reference_header_text)
            cmd = f'kma -i {path}/{file} -o {output}/{file[:-6]}/{file[:-6]} -t_db {database} -t 8 -ont -mem_mode -Mt1 {template_number}'
            os.system(cmd)
        else:
            print(f"No suitable template found for {file}")

print("Process completed.")
