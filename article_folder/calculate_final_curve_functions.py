import json
import numpy as np
from scipy.interpolate import UnivariateSpline
import matplotlib.pyplot as plt

# Define the parameter sets
parameters_list = []

parameters_1 = {
    'cor_interval': [0.3],
    'iteration_increase_interval': [0.51],
    'pp_interval': [0.56],
    'np_interval': [3.3],
    'dp_interval': [0.4]
}
parameters_list.append(parameters_1)

parameters_2 = {
    'cor_interval': [0.34],
    'iteration_increase_interval': [0.23],
    'pp_interval': [0.51],
    'np_interval': [3.0],
    'dp_interval': [0.25]
}
parameters_list.append(parameters_2)

parameters_3 = {
    'cor_interval': [0.37],
    'iteration_increase_interval': [0.24],
    'pp_interval': [0.37],
    'np_interval': [2.8],
    'dp_interval': [0.11]
}
parameters_list.append(parameters_3)

parameters_4 = {
    'cor_interval': [0.58],
    'iteration_increase_interval': [0.03],
    'pp_interval': [0.26],
    'np_interval': [2.7],
    'dp_interval': [0.07]
}
parameters_list.append(parameters_4)

parameters_5 = {
    'cor_interval': [0.65],
    'iteration_increase_interval': [0.01],
    'pp_interval': [0.10],
    'np_interval': [2.65],
    'dp_interval': [0.07]
}
parameters_list.append(parameters_5)

# Initialize lists for each parameter type
cor_intervals = []
iteration_intervals = []
pp_intervals = []
np_intervals = []
dp_intervals = []

# Extract data from parameters list
for params in parameters_list:
    cor_intervals.extend(params['cor_interval'])
    iteration_intervals.extend(params['iteration_increase_interval'])
    pp_intervals.extend(params['pp_interval'])
    np_intervals.extend(params['np_interval'])
    dp_intervals.extend(params['dp_interval'])

# Fit splines and store results
splines = {}
x_values = np.linspace(0.01, 0.05, 5)  # MRD values from 0.01 to 0.05
fine_x_values = np.linspace(0.01, 0.05, 100)  # Use more points for a smoother plot
fig, axes = plt.subplots(nrows=1, ncols=5, figsize=(25, 5))
fig.subplots_adjust(hspace=0.4, wspace=0.4)

for idx, (name, data) in enumerate(zip(['cor', 'iteration', 'pp', 'np', 'dp'],
                                       [cor_intervals, iteration_intervals, pp_intervals, np_intervals, dp_intervals])):
    spline = UnivariateSpline(x_values, data, s=None)
    spline_fit_fine = spline(fine_x_values)  # Use fine_x_values for plotting
    spline_data = {str(x): float(y) for x, y in zip(fine_x_values, spline_fit_fine)}

    # Save each spline fit to a JSON file
    with open(f'{name}_spline_fits.json', 'w') as f:
        json.dump(spline_data, f)

    # Plotting on subplots
    ax = axes[idx]
    ax.plot(x_values, data, 'o', label='Original Data')
    ax.plot(fine_x_values, spline_fit_fine, label='Spline Fit')
    ax.set_title(f'{name} Spline Fit')
    ax.set_xlabel('MRD Values')
    ax.set_ylabel('Parameter Values')
    ax.legend()

# Show all plots at once
plt.show()

print("Spline coefficients have been saved to individual JSON files with MRD values as keys.")
