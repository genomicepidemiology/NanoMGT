# Function to calculate average gene length
def calculate_average_gene_length(fasta_data):
  """
  This function calculates the average gene length from a FASTA string.

  Args:
      fasta_data: A string containing FASTA data in the standard format.

  Returns:
      The average gene length as a float.
  """

  # Split the data into individual genes
  genes = fasta_data.strip().split(">")[1:]

  # Calculate gene lengths, handling multi-line sequences
  gene_lengths = []
  current_gene = ""
  for line in genes:
    lines = line.splitlines()
    # First line is the header, skip it
    sequence = lines[1:]
    current_gene += "".join(sequence)
    gene_lengths.append(len(current_gene.strip()))
    current_gene = ""  # Reset for next gene

  # Calculate average gene length
  if gene_lengths:
    average_length = sum(gene_lengths) / len(gene_lengths)
  else:
    average_length = 0  # Handle case of empty file or no genes

  return average_length

# Load the FASTA data from the file
try:
  with open("rmlst.fsa", "r") as f:
    fasta_data = f.read()
except FileNotFoundError:
  print("Error: File 'rmlst.fsa' not found.")
  exit()

# Calculate and print the average gene length
average_length = calculate_average_gene_length(fasta_data)
print(f"Average gene length: {average_length}")
