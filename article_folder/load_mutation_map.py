def load_mutations_from_files(file_paths):
    """
    Loads mutations from multiple files and returns a dictionary with gene IDs as keys
    and sets of mutations as values.

    :param file_paths: List of file paths to load mutations from
    :return: Dictionary with gene IDs as keys and sets of mutations as values
    """
    mutations_dict = {}

    for file_path in file_paths:
        with open(file_path, 'r') as file:
            lines = file.readlines()
            for i in range(0, len(lines), 2):
                gene_id = lines[i].strip()
                mutations = set(lines[i + 1].strip().split(','))

                if gene_id in mutations_dict:
                    mutations_dict[gene_id].update(mutations)
                else:
                    mutations_dict[gene_id] = mutations

    return mutations_dict


# Example usage
file_paths = ['/home/people/malhal/data/new_nanomgt/majority_variants/minor_variant_maps/major_SRR26643493_minor_SRR23387317.txt', '/home/people/malhal/data/new_nanomgt/majority_variants/minor_variant_maps/major_SRR26643493_minor_SRR27710531.txt']  # Add your file paths here
mutations_dict = load_mutations_from_files(file_paths)

# Printing the dictionary for demonstration
for gene_id, mutations in mutations_dict.items():
    print(f"{gene_id}: {mutations}")
