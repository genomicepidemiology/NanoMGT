import json
import numpy as np
from scipy.interpolate import UnivariateSpline
import matplotlib.pyplot as plt

"""
model = 'mixed'
# Updated parameter sets based on your new values
parameters_list = [
    {'np_interval': [4.092], 'pp_interval': [0.265], 'dp_interval': [0.159], 'iteration_increase_interval': [0.180], 'cor_interval': [0.502]},
    {'np_interval': [3.732], 'pp_interval': [0.274], 'dp_interval': [0.169], 'iteration_increase_interval': [0.121], 'cor_interval': [0.483]},
    {'np_interval': [3.235], 'pp_interval': [0.245], 'dp_interval': [0.174], 'iteration_increase_interval': [0.116], 'cor_interval': [0.453]},
    {'np_interval': [2.811], 'pp_interval': [0.228], 'dp_interval': [0.131], 'iteration_increase_interval': [0.106], 'cor_interval': [0.528]},
    {'np_interval': [2.793], 'pp_interval': [0.218], 'dp_interval': [0.131], 'iteration_increase_interval': [0.103], 'cor_interval': [0.536]}
]


model = 'clean'
# Data for "clean" category
parameters_list = [
    {'np_interval': [3.689], 'pp_interval': [0.279], 'dp_interval': [0.234], 'iteration_increase_interval': [0.156], 'cor_interval': [0.388]},
    {'np_interval': [3.033], 'pp_interval': [0.246], 'dp_interval': [0.228], 'iteration_increase_interval': [0.129], 'cor_interval': [0.424]},
    {'np_interval': [2.400], 'pp_interval': [0.255], 'dp_interval': [0.161], 'iteration_increase_interval': [0.074], 'cor_interval': [0.524]},
    {'np_interval': [2.161], 'pp_interval': [0.215], 'dp_interval': [0.160], 'iteration_increase_interval': [0.055], 'cor_interval': [0.512]},
    {'np_interval': [2.009], 'pp_interval': [0.186], 'dp_interval': [0.144], 'iteration_increase_interval': [0.0497], 'cor_interval': [0.459]}
]

"""
model = 'contaminated'
# Data for "clean" category
parameters_list = [
    {'np_interval': [4.024], 'pp_interval': [0.289], 'dp_interval': [0.213], 'iteration_increase_interval': [0.179], 'cor_interval': [0.453]},
    {'np_interval': [3.780], 'pp_interval': [0.328], 'dp_interval': [0.167], 'iteration_increase_interval': [0.196], 'cor_interval': [0.462]},
    {'np_interval': [3.726], 'pp_interval': [0.280], 'dp_interval': [0.161], 'iteration_increase_interval': [0.102], 'cor_interval': [0.451]},
    {'np_interval': [3.719], 'pp_interval': [0.274], 'dp_interval': [0.151], 'iteration_increase_interval': [0.1301], 'cor_interval': [0.503]},
    {'np_interval': [3.451], 'pp_interval': [0.233], 'dp_interval': [0.145], 'iteration_increase_interval': [0.118], 'cor_interval': [0.513]}
]


# Initialize lists for each parameter type
cor_intervals = []
iteration_intervals = []
pp_intervals = []
np_intervals = []
dp_intervals = []

# Extract data from parameters list
for params in parameters_list:
    cor_intervals.append(params['cor_interval'][0])
    iteration_intervals.append(params['iteration_increase_interval'][0])
    pp_intervals.append(params['pp_interval'][0])
    np_intervals.append(params['np_interval'][0])
    dp_intervals.append(params['dp_interval'][0])

x_values = np.linspace(0.01, 0.05, 5)  # Sequence of indices 0.01 through 0.05 for existing data points
fine_x_values = np.linspace(0.01, 0.05, 500)  # Use 500 points for the desired output
spline_results = {}

for idx, (name, data) in enumerate(zip(['cor', 'ii', 'pp', 'np', 'dp'],
                                       [cor_intervals, iteration_intervals, pp_intervals, np_intervals, dp_intervals])):
    # Calculate weights (optional) - for now, we skip weighting for simplicity
    spline = UnivariateSpline(x_values, data, s=None)
    spline_fit_fine = spline(fine_x_values)  # Use fine_x_values for plotting
    spline_results[name] = {str(x): float(y) for x, y in zip(fine_x_values, spline_fit_fine)}

# Save all spline fits to a single JSON file
with open(f'{model}_parameters.json', 'w') as f:
    json.dump(spline_results, f, indent=4)

print("Spline coefficients have been saved to a single JSON file.")
