import pandas as pd

def parse_parameters(param_str):
    """ Parse the parameter string and return a dictionary of parameter values. """
    try:
        params = param_str.split('_')
        param_dict = {
            'cor': float(params[3]),
            'pp': float(params[5]),
            'np': float(params[7]),
            'dp': float(params[9]),
            'iteration_increase': float(params[-1])
        }
        return param_dict
    except (IndexError, ValueError):
        # Return None or a default dictionary if the parameter string is malformed
        return None

def main():
    #numbers = [1, 2, 3, 4, 5]
    numbers = [1, 2]
    for number in numbers:
        # Load the data
        df = pd.read_csv('first_grid_{}_all.csv'.format(number))

        # Prepare a dictionary to store parameter sums and counts for averaging
        param_sums = {}
        param_counts = {}

        for index, row in df.iterrows():
            # Handle cases where Parameters column might be NaN or not a string
            if pd.isna(row['Parameters']) or not isinstance(row['Parameters'], str):
                continue  # Skip this row if Parameters data is missing or not a string

            # Extract seed and experiment name
            parts = row['experiment'].split('_')
            seed = parts[0]
            experiment = '_'.join(parts[1:])

            # Parse parameters
            params = parse_parameters(row['Parameters'])
            if not params:
                continue  # Skip this row if parameters could not be parsed

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
        output_df = pd.DataFrame(output_data, columns=['experiment', 'cor', 'pp', 'bp', 'dp', 'iteration_increase'])
        output_df.to_csv('final_{}_averaged_parameters.csv'.format(number), index=False)

        print("Averaged parameter CSV has been created successfully.")

if __name__ == '__main__':
    main()
