import pandas as pd

def parse_parameters(param_str):
    """ Parse the parameter string and return a dictionary of parameter values. """
    params = param_str.split('_')
    param_dict = {
        'mrd': float(params[1]),
        'cor': float(params[3]),
        'pp': float(params[5]),
        'bp': int(params[7]),
        'dp': float(params[9])
    }
    return param_dict

def main():
    # Load the data
    df = pd.read_csv('1_all.csv')

    # Prepare a dictionary to store parameter sums and counts for averaging
    param_sums = {}
    param_counts = {}

    for index, row in df.iterrows():
        # Extract seed and experiment name
        parts = row['experiment'].split('_')
        seed = parts[0]
        experiment = '_'.join(parts[1:])

        # Parse parameters
        params = parse_parameters(row['Parameters'])

        # Initialize the experiment dictionary if not present
        if experiment not in param_sums:
            param_sums[experiment] = {key: 0.0 for key in params}
            param_counts[experiment] = 0

        # Sum up the parameters
        for key in params:
            param_sums[experiment][key] += params[key]

        # Count the number of entries for averaging
        param_counts[experiment] += 1

    # Prepare data for output
    output_data = []
    for exp, sums in param_sums.items():
        averages = {key: value / param_counts[exp] for key, value in sums.items()}
        output_data.append([exp] + list(averages.values()))

    # Convert to DataFrame and output to CSV
    output_df = pd.DataFrame(output_data, columns=['experiment', 'mrd', 'cor', 'pp', 'bp', 'dp'])
    output_df.to_csv('1_all_averaged_parameters.csv', index=False)

    print("Averaged parameter CSV has been created successfully.")

if __name__ == '__main__':
    main()
