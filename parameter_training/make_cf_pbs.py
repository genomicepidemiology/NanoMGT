import os

# Paths to the main folders
folders = [
    "/home/projects/cge/people/malhal/120_validation_reads",
    "/home/projects/cge/people/malhal/nanomgt_validation_reads/confindr_clean_data",
    "/home/projects/cge/people/malhal/170_validation_reads"
]


def species_from_folders(folder_path):
    """ Extract species based on subfolder names and store them in a dictionary. """
    species_dict = {}
    for subdir in os.listdir(folder_path):
        if os.path.isdir(os.path.join(folder_path, subdir)):
            species_name = '_'.join(subdir.split('_')[:2])
            if species_name in species_dict:
                species_dict[species_name].append(subdir)
            else:
                species_dict[species_name] = [subdir]
    return species_dict


def create_pbs_scripts(folder_path, species_dict):
    """ Generate PBS scripts for each species with different bf rates. """
    folder_type = os.path.basename(folder_path)  # assuming 'confindr_mixed_data' structure
    for species, subfolders in species_dict.items():
        for bf_value in [0.01, 0.02, 0.03, 0.04, 0.05]:
            pbs_filename = f"{folder_type}_{species}_bf_{bf_value:.2f}.pbs"
            pbs_path = os.path.join(os.getcwd(), pbs_filename)  # Save in the current working directory
            with open(pbs_path, 'w') as f:
                f.write("#!/bin/sh\n")
                f.write("#PBS -W group_list=cge -A cge\n")
                f.write(f"#PBS -N {folder_type}_confindr_{species}_bf_{bf_value:.2f}\n")
                f.write(f"#PBS -e {os.getcwd()}/{folder_type}_{species}_bf_{bf_value:.2f}_error.txt\n")
                f.write(f"#PBS -o {os.getcwd()}/{folder_type}_{species}_bf_{bf_value:.2f}_output.log\n")
                f.write("#PBS -l nodes=1:ppn=20\n")
                f.write("#PBS -l mem=64gb\n")
                f.write("#PBS -l walltime=48:00:00\n\n")
                f.write("cd /home/projects/cge/people/malhal/nanomgt_article_results/clean\n")
                f.write("module load tools\n")
                f.write("module load miniconda3/23.5.0\n")
                f.write("conda activate mambaenv\n\n")

                result_base_path = f"/home/projects/cge/people/malhal/nanomgt_depth_results/{folder_type}"
                bf_path = os.path.join(result_base_path, f"bf_{bf_value:.2f}")
                os.makedirs(bf_path, exist_ok=True)

                # Loop through subfolders of the species
                f.write("for folder in " + " ".join([os.path.join(folder_path, sf) for sf in subfolders]) + ";\n")
                f.write("do\n")
                f.write("  base_folder=$(basename $folder)\n")  # This line correctly defines the basename
                output_path = os.path.join(bf_path, "$base_folder")
                os.makedirs(output_path, exist_ok=True)
                f.write(
                    f"  confindr -i $folder -o {output_path} -d /home/projects/cge/people/malhal/rmlst_db/ --rmlst -b 3 -bf {bf_value:.2f} -dt Nanopore -q 14 -t 8\n")
                f.write("done\n")
                print(f"Generated {pbs_path}")


# Process each folder
for folder in folders:
    species_dict = species_from_folders(folder)
    create_pbs_scripts(folder, species_dict)