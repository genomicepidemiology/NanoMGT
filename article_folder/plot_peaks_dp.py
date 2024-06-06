import pandas as pd
import numpy as np
from scipy.interpolate import UnivariateSpline
import matplotlib.pyplot as plt

def process_data(csv_file, maf_value=0.01, threshold=5):  # threshold in percentage
    data = pd.read_csv(csv_file)
    grouped = data.groupby(['maf', 'batch'])

    # Prepare a figure to host multiple subplots
    fig, axes = plt.subplots(nrows=2, ncols=5, figsize=(20, 8))  # Adjusted for 10 plots (2 rows x 5 columns)
    fig.subplots_adjust(hspace=0.4, wspace=0.4)  # Adjust spacing to prevent label overlap

    for idx, batch in enumerate(range(1, 11)):  # Loop over batches 1 to 10
        group = grouped.get_group((maf_value, batch))
        if group is not None:
            dp_values = group['dp'].values
            f1_sdpes = group['average_f1'].values

            # Fit a spline to the data points
            spline = UnivariateSpline(dp_values, f1_sdpes, s=None)
            derivative = spline.derivative()

            # Generate dense dp values for a finer analysis
            dp_dense = np.linspace(dp_values.min(), dp_values.max(), 300)
            f1_dense = spline(dp_dense)
            derivative_values = derivative(dp_dense)

            # Calculate the percentage change in derivative values
            pct_change = np.abs(np.diff(derivative_values) / derivative_values[:-1]) * 100

            # Find the index where the percentage change drops below the threshold
            significant_change = np.where(pct_change < threshold)[0]
            diminishing_point = significant_change[0] + 1 if significant_change.size > 0 else len(dp_dense) - 1

            # Plotting on subplots
            ax = axes[idx // 5, idx % 5]  # Determine the row and column by indexing
            ax.plot(dp_values, f1_sdpes, 'bo', label='Original Data')
            ax.plot(dp_dense, f1_dense, 'r-', label='Spline Fit')
            ax.axvline(x=dp_dense[diminishing_point], color='g', label='Diminishing Returns Start')
            ax.set_title(f'Batch: {batch}')
            ax.set_xlabel('dp Values')
            ax.set_ylabel('F1 Sdpes')
            ax.legend()

    plt.show()  # Show all plots at once

def main():
    csv_file = 'dp_f1_scores.csv'
    threshold = 5  # Percentage change threshold
    maf_value = 0.04  # Set the constant maf value
    process_data(csv_file, maf_value, threshold)

if __name__ == "__main__":
    main()
