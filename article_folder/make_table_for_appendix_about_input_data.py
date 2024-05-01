from pylatex import Document, Tabularx, Command, NoEscape, NewPage

# Create a new document
doc = Document(geometry_options={"margin": "1in"})

# Define the structure and content of the table
data = [
    ["Sample", "Specie", "Flow Cell Chemistry", "Base Calling Model", "Version"],
    ["SRR25689478", "Escherichia coli", "R10.4.1", "Super-accurate basecalling", "Guppy 6.5.7"],
    ["SRR26899125", "Escherichia coli", "R10.4.1", "Super-accurate basecalling", "Guppy 4.2.0"],
    ["SRR28370694", "Staphylococcus aureus", "R10.4.1", "Super-accurate basecalling", "Guppy 4.3.0"],
    ["ERR8958843", "Staphylococcus aureus", "R10.4", "Super-accurate basecalling", "Guppy 5.0.12"],
    ["SRR27167517", "Staphylococcus aureus", "R10.4.1", "Super-accurate basecalling", "Guppy 3.5.1"],
    ["SRR27638397", "Campylobacter jejuni", "R10.4.1", "Super-accurate basecalling", "Guppy 4.3.0"],
    ["SRR27710526", "Campylobacter jejuni", "R10.4.1", "Super-accurate basecalling", "Guppy 7.1.4"],
    ["SRR28399428", "Salmonella enterica", "R10.4.1", "Super-accurate basecalling", "Dorado 7.2.13"],
    ["SRR27136088", "Salmonella enterica", "R10.4.1", "Super-accurate basecalling", "Guppy 5.7.11"],
    ["SRR27755678", "Salmonella enterica", "R10.4.1", "Super-accurate basecalling", "Guppy 4.3.0"],
    ["ERR8958737", "Klebsiella pneumoniae", "R10.4", "Super-accurate basecalling", "Guppy 5.0.12"],
    ["SRR27348733", "Klebsiella pneumoniae", "R10.4.1", "Super-accurate basecalling", "Guppy 7.1.4"]
]

# Initialize the table with the specified format and column alignment
with doc.create(Tabularx("1\\textwidth", arguments="|Y|Y|Y|Y|Y|", pos="htbp")) as table:
    table.add_hline()
    table.add_row(data[0], mapper=bold_mapper, color="gray")
    table.add_hline(2, -1)
    for row in data[1:]:
        table.add_row(row)
        table.add_hline()
    table.add_hline()

# Add a caption and a label for the table
doc.append(Command('caption', 'Sequencing and Base Calling Information'))
doc.append(Command('label', 'tab:sequencing_info'))
doc.append(NoEscape(r'\textit{Note:} This table provides detailed sequencing and base calling information for various samples across multiple species, highlighting the flow cell chemistry and software versions used in the analysis.'))

# Generate the PDF
doc.generate_pdf("data_appendix.pdf", clean_tex=False)
