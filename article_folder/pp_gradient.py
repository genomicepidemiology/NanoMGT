import pandas as pd
import numpy as np
from scipy.interpolate import UnivariateSpline

def process_data(csv_file):
    data = pd.read_csv(csv_file)
    grouped = data.groupby(['maf', 'batch'])
    results = []

    for (maf, batch), group in grouped:
        pp_values = group['pp'].values  # Changed from cor_values to pp_values
        f1_scores = group['average_f1'].values

        # Normalize pp_values and f1_scores to range [0, 1]
        pp_values_normalized = (pp_values - np.min(pp_values)) / (np.max(pp_values) - np.min(pp_values))
        f1_scores_normalized = (f1_scores - np.min(f1_scores)) / (np.max(f1_scores) - np.min(f1_scores))

        # Fit a spline to the normalized data points
        spline = UnivariateSpline(pp_values_normalized, f1_scores_normalized, s=None)
        derivative = spline.derivative()

        # Generate dense pp values for a finer analysis in the normalized range
        pp_dense_normalized = np.linspace(0, 1, 450)
        f1_dense_normalized = spline(pp_dense_normalized)
        derivative_values_normalized = derivative(pp_dense_normalized)

        # Calculate the target derivative for a 30-degree angle
        target_slope = np.tan(np.radians(30))

        # Initialize the last valid index
        last_valid_index = 0

        # Iterate over all points to find derivatives close to the target slope
        for idx in range(len(derivative_values_normalized)):
            if target_slope > 0:  # Ensure the slope is positive
                if abs(derivative_values_normalized[idx] - target_slope) < 0.01:  # Close to 30-degree slope
                    if f1_dense_normalized[idx] > f1_dense_normalized[last_valid_index]:  # Ensure F1-score is increasing
                        last_valid_index = idx

        if last_valid_index == 0:
            if f1_dense_normalized[0] > f1_dense_normalized[-1]:
                last_valid_index = 0
            else:
                last_valid_index = -1
        # Get the pp value at the last valid index found
        target_pp_value = pp_values.min() + pp_dense_normalized[last_valid_index] * (pp_values.max() - pp_values.min())
        results.append((maf, batch, target_pp_value))

    return results

def main():
    csv_file = 'pp_f1_scores.csv'  # Adjusted to use the correct CSV file
    results = process_data(csv_file)
    for result in results:
        maf, batch, peak_value = result
        print(f"The pp value with a 30-degree slope for MAF: {maf:.2f}, Batch: {batch} is at pp: {peak_value}")

if __name__ == "__main__":
    main()
