import pandas as pd
import numpy as np
from scipy.interpolate import UnivariateSpline

def process_data(csv_file):
    data = pd.read_csv(csv_file)
    grouped = data.groupby(['mrd', 'batch'])
    results = []

    for (mrd, batch), group in grouped:
        ii_values = group['ii'].values
        f1_scores = group['average_f1'].values

        # Normalize ii_values and f1_scores to range [0, 1]
        ii_values_normalized = (ii_values - np.min(ii_values)) / (np.max(ii_values) - np.min(ii_values))
        f1_scores_normalized = (f1_scores - np.min(f1_scores)) / (np.max(f1_scores) - np.min(f1_scores))

        # Fit a spline to the normalized data points
        spline = UnivariateSpline(ii_values_normalized, f1_scores_normalized, s=None)
        derivative = spline.derivative()

        # Generate dense ii values for a finer analysis in the normalized range
        ii_dense_normalized = np.linspace(0, 1, 450)
        f1_dense_normalized = spline(ii_dense_normalized)
        derivative_values_normalized = derivative(ii_dense_normalized)

        # Calculate the target derivative for a 30-degree angle
        target_slope = np.tan(np.radians(30))

        # Collect all ii values with a derivative close to the target slope and ensure they are ascending
        valid_ii_values = []
        for idx in range(1, len(derivative_values_normalized)):
            if abs(derivative_values_normalized[idx] - target_slope) < 0.01 and f1_dense_normalized[idx] > f1_dense_normalized[idx - 1]:
                valid_ii_value = ii_values.min() + ii_dense_normalized[idx] * (ii_values.max() - ii_values.min())
                valid_ii_values.append(valid_ii_value)

        results.append((mrd, batch, valid_ii_values))

    return results

def main():
    csv_file = 'ii_f1_scores.csv'  # Ensure the CSV file is correctly named
    results = process_data(csv_file)
    for result in results:
        mrd, batch, peak_values = result
        if peak_values:
            print(f"Batch {batch} for MRD {mrd:.2f} has valid ii values at slopes of 30 degrees: {peak_values}")
        else:
            print(f"No valid ii values found for MRD {mrd:.2f}, Batch {batch}")

if __name__ == "__main__":
    main()
